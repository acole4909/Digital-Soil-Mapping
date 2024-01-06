##### AGDS II: Digital Soil Mapping Exercise
#### Exercise 5.1
library(dplyr)
library(ggplot2)
library(tidyverse)
library(terra)
library(here)
library(caret)
library(recipes)
library(plotROC)

### Load Data
df_full <- readRDS(here::here("data/df_full.rds"))

head(df_full) |>
  knitr::kable()

## Specify target: Waterlog.100
target <- "waterlog.100"

## Specify predictors_all
predictors_all <- names(df_full)[14:ncol(df_full)]

cat("The target is:", target,
    "\nThe predictors_all are:", paste0(predictors_all[1:8], sep = ", "), "...")

# Split dataset into training and testing sets
set.seed(42)
df_train <- df_full |> dplyr::filter(dataset == "calibration")
df_test  <- df_full |> dplyr::filter(dataset == "validation")

# Filter out any NAs
df_train <- df_train |> tidyr::drop_na()
df_test <- df_test   |> tidyr::drop_na()

n_tot <- nrow(df_train) + nrow(df_test)

perc_cal <- (nrow(df_train) / n_tot) |> round(2) * 100
perc_val <- (nrow(df_test)  / n_tot) |> round(2) * 100

cat("For model training, we have a calibration / validation split of: ",
    perc_cal, "/", perc_val, "%")

### Train Model
rf_basic <- ranger::ranger(
  probability = FALSE,
  classification = TRUE,
  y = df_train[, target],     # target variable
  x = df_train[, predictors_all], # Predictor variables
  seed = 42,                    # Specify the seed for randomization to reproduce the same model again
  num.threads = parallel::detectCores() - 1) # Use all but one CPU core for quick model training

# Print summary
print(rf_basic)

#Save Data
saveRDS(rf_basic,
        here::here("data/rf_basic_for_waterlog100.rds"))

saveRDS(df_train[, c(target, predictors_all)],
        here::here("data/cal_basic_for_waterlog100.rds"))

saveRDS(df_test[, c(target, predictors_all)],
        here::here("data/val_basic_for_waterlog100.rds"))

#### Model Analysis
# Load random forest model from data folder
rf_basic   <- readRDS(here::here("data/rf_basic_for_waterlog100.rds"))
df_train <- readRDS(here::here("data/cal_basic_for_waterlog100.rds"))
df_test  <- readRDS(here::here("data/val_basic_for_waterlog100.rds"))

# Load area to be predicted
set.seed(42)
raster_mask <- terra::rast(here::here("data-raw/geodata/study_area/area_to_be_mapped.tif"))
# Turn target raster into a dataframe, 1 px = 1 cell
df_mask <- as.data.frame(raster_mask, xy = TRUE)

# Filter for area of interest
df_mask <- df_mask |>
  dplyr::filter(area_to_be_mapped == 1)

# Display data frame
head(df_mask) |>
  knitr::kable()

files_covariates <- list.files(
  path = here::here("data-raw/geodata/covariates/"),
  pattern = ".tif$",
  recursive = TRUE,
  full.names = TRUE
)

# Filter for variables used in RF
preds_all <- names((df_train[, predictors_all]))
files_selected <- files_covariates[apply(sapply(X = preds_all,
                                                FUN = grepl,
                                                files_covariates),
                                         MARGIN =  1,
                                         FUN = any)]

# Load all rasters as a stack
raster_covariates <- terra::rast(files_selected)

# Get coordinates
df_locations <- df_mask |>
  dplyr::select(x, y)

# Extract data from covariate raster stack for all gridcells in the raster
df_predict <- terra::extract(
  raster_covariates,
  df_locations,
  ID = FALSE
)

df_predict <- cbind(df_locations, df_predict) |>
  tidyr::drop_na()

# Make predictions for validation sites
prediction <- predict(
  rf_basic,
  data = df_test,
  num.threads = parallel::detectCores() - 1
)

# Save predictions to validation data frame
df_test$pred <- prediction$predictions

# Classification Metrics
Y <- df_test$waterlog.100
Y <- as.factor(Y)
X <- df_test$pred
X <- as.factor(X)

#Confusion Matrix
conf_matrix_waterlog_rfbasic <- caret::confusionMatrix(data=X, reference=Y, positive="1")
conf_matrix_waterlog_rfbasic

mosaicplot(conf_matrix_waterlog_rfbasic$table,
           main = "Confusion matrix")

### Create Prediction Maps
df_predict$vszone <- trimws(df_predict$vszone) #vszone tif was reading in with missing values, this fixed the issue.
prediction <- predict(
  rf_basic,
  data = df_predict,
  num.threads = parallel::detectCores() - 1)

# Attach predictions to dataframe
df_predict$prediction <- prediction$predictions

# Extract dataframe with coordinates and predictions
df_map <- df_predict |>
  dplyr::select(x, y, prediction)

# Turn dataframe into raster
raster_pred <- terra::rast(
  df_map,                  # Table to be transformed
  crs = "+init=epsg:2056", # Swiss coordinate system
  extent = terra::ext(raster_covariates) # Prescribe same extent as predictor rasters
)

ggplot2::ggplot() +
  tidyterra::geom_spatraster(data = raster_pred) +
  ggplot2::scale_fill_viridis_c(
    breaks = seq(0,1, by=1),
    labels = c("0", "1"),
    na.value = NA,
    option = "viridis",
    name = "Waterlog.100",
  ) +
  ggplot2::theme_classic() +
  ggplot2::scale_x_continuous(expand = c(0, 0)) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::labs(title = "Predicted Waterlog")

# Save raster as .tif file
if (!dir.exists(here::here("data"))) system(paste0("mkdir ", here::here("data")))
terra::writeRaster(
  raster_pred,
  "data/ra_predicted_waterlog100.tif",
  datatype = "FLT4S",
  filetype = "GTiff",
  overwrite = TRUE
)

#### Exercise 5.2
#### Data Loading and Preparation
### Load Data
df_full <- readRDS(here::here("data/df_full.rds"))
df_train <- readRDS(here::here("data/cal_basic_for_waterlog100.rds"))
df_test  <- readRDS(here::here("data/val_basic_for_waterlog100.rds"))

## Specify target: Waterlog.100
target <- "waterlog.100"

## Specify predictors_all
predictors_all <- names(df_full)[14:ncol(df_full)]

### Variable Importance
set.seed(42)
rf_varimport <- ranger::ranger(
  probability = FALSE,
  y = df_train[, target],     # target variable
  x = df_train[, predictors_all],   # Predictor variables
  importance   = "permutation", # Pick permutation to calculate variable importance
  classification = TRUE,
  seed = 42,                    # Specify seed for randomization to reproduce the same model again
  num.threads = parallel::detectCores() - 1) # Use all but one CPU core for quick model training

print(rf_varimport)

# Extract the variable importance and create a long tibble
vi_rf_varimport <- rf_varimport$variable.importance |>
  dplyr::bind_rows() |>
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Variable")
vi_rf_varimport_top5 <- vi_rf_varimport |>
  top_n(5, value)
vi_rf_varimport_top5 <- vi_rf_varimport_top5 |>
  arrange(desc(value))
vi_rf_varimport_top5 <- vi_rf_varimport_top5 |>
  mutate(Importance = row_number())
vi_rf_varimport_top5 <- vi_rf_varimport_top5 |>
  rename(Value = value)
vi_rf_varimport_top5 <- vi_rf_varimport_top5 |>
  select(Importance, Variable, Value)

library(knitr)
kable(vi_rf_varimport_top5)

# Plot variable importance, ordered by decreasing value
gg <- vi_rf_varimport |>
  ggplot2::ggplot(ggplot2::aes(x = reorder(Variable, value), y = value)) +
  ggplot2::geom_bar(stat = "identity", fill = "grey50", width = 0.75) +
  ggplot2::labs(
    y = "Change in OOB accuracy after permutation",
    x = "",
    title = "Variable importance based on OOB") +
  ggplot2::theme(text=element_text(size=7))+
  ggplot2::coord_flip()

# Display plot
gg

### Variable Selection
set.seed(42)

# Run the algorithm
bor <- Boruta::Boruta(
  y = df_train[, target],
  x = df_train[, predictors_all],
  maxRuns = 50,
  seed = 42,
  num.threads = parallel::detectCores()-1)

df_bor <- Boruta::attStats(bor) |>
  tibble::rownames_to_column() |>
  dplyr::arrange(dplyr::desc(meanImp))

# Plot Result
ggplot2::ggplot(ggplot2::aes(x = reorder(rowname, meanImp),
                             y = meanImp,
                             fill = decision),
                data = df_bor) +
  ggplot2::geom_bar(stat = "identity", width = 0.75) +
  ggplot2::scale_fill_manual(values = c("grey30", "tomato", "grey70")) +
  ggplot2::labs(
    y = "Variable importance",
    x = "",
    title = "Variable importance based on Boruta") +
  ggplot2::theme(text=element_text(size=7)) +
  ggplot2::coord_flip()

# Save retained important variables
predictors_selected <- df_bor |>
  dplyr::filter(decision == "Confirmed") |>
  dplyr::pull(rowname)

# Re-train Random Forest model
set.seed(42)
rf_bor <- ranger::ranger(
  probability = FALSE,
  y = df_train[, target],              # target variable
  x = df_train[, predictors_selected], # Predictor variables
  classification = TRUE,
  seed = 42,                           # Specify the seed for randomization to reproduce the same model again
  num.threads = parallel::detectCores() - 1)
print(rf_bor)

# Save Data
saveRDS(rf_bor,
        here::here("data/rf_bor_for_waterlog.100.rds"))

saveRDS(df_train[, c(target, predictors_selected)],
        here::here("data/cal_bor_for_waterlog.100.rds"))

saveRDS(df_test[, c(target, predictors_selected)],
        here::here("data/val_bor_for_waterlog.100.rds"))

### Model Evaluation
# Load random forest model
rf_bor   <- readRDS(here::here("data/rf_bor_for_waterlog.100.rds"))
df_train <- readRDS(here::here("data/cal_bor_for_waterlog.100.rds"))
df_test  <- readRDS(here::here("data/val_bor_for_waterlog.100.rds"))

# Load area to be predicted
set.seed(42)
raster_mask <- terra::rast(here::here("data-raw/geodata/study_area/area_to_be_mapped.tif"))

# Turn target raster into a dataframe, 1 px = 1 cell
df_mask <- as.data.frame(raster_mask, xy = TRUE)

# Filter for area of interst
df_mask <- df_mask |>
  dplyr::filter(area_to_be_mapped == 1)

head(df_mask) |>
  knitr::kable()

files_covariates <- list.files(
  path = here::here("data-raw/geodata/covariates/"),
  pattern = ".tif$",
  recursive = TRUE,
  full.names = TRUE
)

# Filter for variables used in Random Forest model
preds_selected <- names(df_train[, predictors_selected])
files_selected <- files_covariates[apply(sapply(X = preds_selected,
                                                FUN = grepl,
                                                files_covariates),
                                         MARGIN =  1,
                                         FUN = any)]

# Load all rasters as a stack
raster_covariates <- terra::rast(files_selected)

# Get coordinates
df_locations <- df_mask |>
  dplyr::select(x, y)

# Extract data from covariate raster stack for all gridcells in the raster
df_predict <- terra::extract(
  raster_covariates,   # The raster we want to extract from
  df_locations,        # A matrix of x and y values to extract for
  ID = FALSE           # To not add a default ID column to the output
)

df_predict <- cbind(df_locations, df_predict) |>
  tidyr::drop_na()

# Make predictions for validation sites
set.seed(42)
prediction <- predict(
  rf_bor,           # RF model
  data = df_test,   # Predictor data
  num.threads = parallel::detectCores() - 1,
  seed = 42
)

# Save predictions to validation data frame
df_test$pred <- prediction$predictions

# Classification Metrics
Y <- df_test$waterlog.100
Y <- as.factor(Y)
X <- df_test$pred
X <- as.factor(X)

#Confusion Matrix
conf_matrix_waterlog_bor <- caret::confusionMatrix(data=X, reference=Y, positive="1")
conf_matrix_waterlog_bor
mosaicplot(conf_matrix_waterlog_bor$table,
           main = "Confusion matrix")

### Create Prediction Maps
prediction <- predict(
  rf_bor,              # RF model
  data = df_predict,
  num.threads = parallel::detectCores() - 1)

# Attach predictions to dataframe
df_predict$prediction <- prediction$predictions

# Extract dataframe with coordinates and predictions
df_map <- df_predict |>
  dplyr::select(x, y, prediction)

# Turn dataframe into a raster
raster_pred <- terra::rast(
  df_map,                  # Table to be transformed
  crs = "+init=epsg:2056", # Swiss coordinate system
  extent = terra::ext(raster_covariates) # Prescribe same extent as predictor rasters
)

# Visualise predictions
ggplot2::ggplot() +
  tidyterra::geom_spatraster(data = raster_pred) +
  ggplot2::scale_fill_viridis_c(
    breaks = seq(0,1, by=1),
    labels = c("0", "1"),
    na.value = NA,
    option = "viridis",
    name = "Waterlog.100"
  ) +
  ggplot2::theme_classic() +
  ggplot2::scale_x_continuous(expand = c(0, 0)) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::labs(title = "Predicted Waterlog")

# Save raster as .tif file
if (!dir.exists(here::here("data"))) system(paste0("mkdir ", here::here("data")))
terra::writeRaster(
  raster_pred,
  "data/ra_predicted_rfbor_waterlog100.tif",
  datatype = "FLT4S",  # FLT4S for floats, INT1U for integers (smaller file)
  filetype = "GTiff",  # GeoTiff format
  overwrite = TRUE     # Overwrite existing file
)

#### Exercise 5.3
# Load Boruta random forest model
rf_bor   <- readRDS(here::here("data/rf_bor_for_waterlog.100.rds"))
df_train <- readRDS(here::here("data/cal_bor_for_waterlog.100.rds"))
df_test  <- readRDS(here::here("data/val_bor_for_waterlog.100.rds"))

# Set target
target <- "waterlog.100"
# Specify predictors_selected
preds_selected <- names(df_train[, predictors_selected])

cat("The target is:", target,
    "\nThe predictors_selected are:", paste0(preds_selected[1:8], sep = ", "), "...")

#### Hyperparameter Tuning
# Recipes for Train
set.seed(42)
df_train$waterlog.100 <- as.factor(df_train$waterlog.100)
pp <- recipes::recipe(waterlog.100 ~ NegO + mrvbf25 + mt_rr_y + Se_diss2m_50c + Se_TWI2m + Se_curv2m_s60 + be_gwn25_vdist + Se_MRVBF2m + lsf +
                        Se_TWI2m_s15 + cindx10_25 + mt_td_y + Se_tpi_2m_50c + terrTextur + tsc25_40 + Se_NO2m_r500 + Se_curv2m_fmean_50c + Se_TWI2m_s60 +
                        Se_slope50m + tsc25_18 + vdcn25 + Se_alti2m_std_50c + mt_tt_y + Se_curv50m + be_gwn25_hdist + Se_rough2m_10c + Se_curvplan2m_s60 +
                        Se_diss2m_5c + Se_curvprof50m + Se_slope6m + Se_rough2m_5c + Se_curvplan50m + Se_slope2m_s7 + Se_slope2m_fmean_5c +
                        Se_slope2m_fmean_50c, data = df_train)

# Greedy Hyperparameter Tuning Approach
mtry_values <- c(2,3,4,5,6,7,8,9,10,12,14,16)
min.node.size_values <- c(5,10,20,25,30,35,40,45,50,55,60)
splitrule_values <- c("gini", "extratrees")
set.seed(42)
mod <- caret::train(
  pp,
  data = df_train %>%
    drop_na(),
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
  tuneGrid = expand.grid( .mtry = 6,
                          .min.node.size = min.node.size_values,
                          .splitrule = "gini"),
  metric = "Accuracy",
  replace = FALSE,
  sample.fraction = 0.5,
  num.trees = 100,
  seed = 42,                # for reproducibility
  num.threads = parallel::detectCores() - 1)
print(mod)
# Best min.node.size was 40 with mtry of 6 and splitrule gini

set.seed(42)
mod <- caret::train(
  pp,
  data = df_train %>%
    drop_na(),
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
  tuneGrid = expand.grid( .mtry = mtry_values,
                          .min.node.size = 40,
                          .splitrule = "gini"),
  metric = "Accuracy",
  replace = FALSE,
  sample.fraction = 0.5,
  num.trees = 100,
  seed = 42                # for reproducibility
)
print(mod)
# Best mtry value was 14 with min.node.size 5 and splitrule gini
set.seed(42)
mod <- caret::train(
  pp,
  data = df_train %>%
    drop_na(),
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
  tuneGrid = expand.grid( .mtry = 6,
                          .min.node.size = 40,
                          .splitrule = splitrule_values),
  metric = "Accuracy",
  replace = FALSE,
  sample.fraction = 0.5,
  num.trees = 100,
  seed = 42                # for reproducibility
)
print(mod)
# Best combination of hyperparameters found by greedy hyperparameter was mtry value=14, min.node.size=5, and splitrule=gini
set.seed(42)
mod_greedy <- caret::train(
  pp,
  data = df_train |>
    drop_na(),
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
  tuneGrid = expand.grid( .mtry = 6,
                          .min.node.size = 40,
                          .splitrule = "extratrees"),
  metric = "Accuracy",
  replace = FALSE,
  sample.fraction = 0.5,
  num.trees = 100,
  seed = 42                # for reproducibility
)
print(mod_greedy)

## Grid Hyperparameter Tuning
set.seed(42)
mod_grid <- caret::train(
  pp,
  data = df_train|>
    drop_na(),
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
  tuneGrid = expand.grid( .mtry = mtry_values,
                          .min.node.size = min.node.size_values,
                          .splitrule = splitrule_values),
  metric = "Accuracy",
  replace = FALSE,
  sample.fraction = 0.5,
  num.trees = 100,
  seed = 42                # for reproducibility
)
print(mod_grid)

# Save relevant data for model testing.
saveRDS(mod_greedy,
        here::here("data/rf_mod_greedy_for_waterlog.100.rds"))

saveRDS(mod_grid,
        here::here("data/rf_mod_grid_for_waterlog.100.rds"))

saveRDS(df_train[, c(target, predictors_selected)],
        here::here("data/cal_mod_for_waterlog.100.rds"))

saveRDS(df_test[, c(target, predictors_selected)],
        here::here("data/val_mod_for_waterlog.100.rds"))

#### Model Analysis
### Load Model and Data
mod_greedy   <- readRDS(here::here("data/rf_mod_greedy_for_waterlog.100.rds"))
mod_grid   <- readRDS(here::here("data/rf_mod_grid_for_waterlog.100.rds"))
df_train <- readRDS(here::here("data/cal_bor_for_waterlog.100.rds"))
df_test  <- readRDS(here::here("data/val_bor_for_waterlog.100.rds"))

### Evaluation
# Make predictions for validation sites using greedy model
set.seed(42)
prediction <- predict(
  mod_greedy,
  newdata = df_test,
  seed = 42,
  num.threads = parallel::detectCores() - 1
)

# Save predictions to validation df
df_test$predcv <- prediction

# Classification Metrics
Y <- df_test$waterlog.100
Y <- as.factor(Y)
X <- df_test$predcv
X <- as.factor(X)

conf_matrix_waterlog_rfmodcv <- caret::confusionMatrix(data=X, reference=Y, positive="1")
conf_matrix_waterlog_rfmodcv
mosaicplot(conf_matrix_waterlog_rfmodcv$table,
           main = "Confusion matrix")

# Make predictions for validation sites using grid model
set.seed(42)
prediction <- predict(
  mod_grid,
  newdata = df_test,
  seed = 42,
  num.threads = parallel::detectCores() - 1
)

# Save predictions to validation df
df_test$predcv <- prediction

# Classification Metrics
Y <- df_test$waterlog.100
Y <- as.factor(Y)
X <- df_test$predcv
X <- as.factor(X)

conf_matrix_waterlog_rfmodcv <- caret::confusionMatrix(data=X, reference=Y, positive="1")
conf_matrix_waterlog_rfmodcv
mosaicplot(conf_matrix_waterlog_rfmodcv$table,
           main = "Confusion matrix")


#### Exercise 5.4
set.seed(42)
rf_prob <- ranger::ranger(
  probability = TRUE,
  classification = TRUE,
  y = df_train[, "waterlog.100"],     # Target variable
  x = df_train[, predictors_selected], # Predictor variables
  seed = 42,
  num.threads = parallel::detectCores() - 1)

# Print summary
print(rf_prob)

# Save Data
saveRDS(rf_prob,
        here::here("data/rf_prob_for_waterlog.100.rds"))

saveRDS(df_train[, c(target, predictors_selected)],
        here::here("data/cal_prob_for_waterlog.100.rds"))

saveRDS(df_test[, c(target, predictors_selected)],
        here::here("data/val_prob_for_waterlog.100.rds"))

###Prob Model Analysis
# Load random forest model
rf_prob   <- readRDS(here::here("data/rf_prob_for_waterlog.100.rds"))
df_train <- readRDS(here::here("data/cal_prob_for_waterlog.100.rds"))
df_test  <- readRDS(here::here("data/val_prob_for_waterlog.100.rds"))

# Make predictions for validation sites
set.seed(42)
prediction <- predict(
  rf_prob,           # RF model
  data = df_test,   # Predictor data
  type="response",
  seed=42,
  num.threads = parallel::detectCores() - 1
)

# Save predictions to validation df
df_test$pred <- prediction$predictions
predictions <- prediction$predictions

#ROC Curve
basicplot <- ggplot(df_test, aes(d=waterlog.100, m=pred[,1])) + geom_roc()
basicplot

