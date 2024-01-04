##### Exercise 5: 5.3
#### Data Loading and Preparation
library(caret)
library(recipes)
library(dplyr)
library(tidyverse)

# Load Borruta random forest model
rf_bor   <- readRDS(here::here("data/rf_bor_for_waterlog.100.rds"))
df_train <- readRDS(here::here("data/cal_bor_for_waterlog.100.rds"))
df_test  <- readRDS(here::here("data/val_bor_for_waterlog.100.rds"))

# Set target
target <- "waterlog.100"
target <- as.factor(target)
# Specify predictors_selected
preds_selected <- names(df_train[, predictors_selected])

cat("The target is:", target,
    "\nThe predictors_selected are:", paste0(preds_selected[1:8], sep = ", "), "...")

# Split dataset into training and testing sets
df_train <- df_full |> dplyr::filter(dataset == "calibration")
df_test  <- df_full |> dplyr::filter(dataset == "validation")

# Filter out any NA to avoid error when running a Random Forest
df_train <- df_train |> tidyr::drop_na()
df_test <- df_test   |> tidyr::drop_na()

# A little bit of verbose output:
n_tot <- nrow(df_train) + nrow(df_test)

perc_cal <- (nrow(df_train) / n_tot) |> round(2) * 100
perc_val <- (nrow(df_test)  / n_tot) |> round(2) * 100

cat("For model training, we have a calibration / validation split of: ",
    perc_cal, "/", perc_val, "%")
set.seed(42)

# Recipes for Train
df_train$waterlog.100 <- as.factor(df_train$waterlog.100)
is.factor(df_train$waterlog.100)
pp <- recipes::recipe(waterlog.100 ~ NegO + mrvbf25 + mt_rr_y + Se_diss2m_50c + Se_TWI2m + Se_curv2m_s60 + be_gwn25_vdist + Se_MRVBF2m + lsf +
                        Se_TWI2m_s15 + cindx10_25 + mt_td_y + Se_tpi_2m_50c + terrTextur + tsc25_40 + Se_NO2m_r500 + Se_curv2m_fmean_50c + Se_TWI2m_s60 +
                        Se_slope50m + tsc25_18 + vdcn25 + Se_alti2m_std_50c + mt_tt_y + Se_curv50m + be_gwn25_hdist + Se_rough2m_10c + Se_curvplan2m_s60 +
                        Se_diss2m_5c + Se_curvprof50m + Se_slope6m + Se_rough2m_5c + Se_SCA2m + Se_curvplan50m + Se_slope2m_s7 + Se_slope2m_fmean_5c +
                        Se_slope2m_fmean_50c + Se_slope2m_s60, data = df_train)

## Hyperparameter Tuning
mtry_values <- c(2,3,4,5,6,7,8,9,10,12,14,16)
min.node.size_values <- c(2,5,10,20,25)
splitrule_values <- c("gini", "extratrees")
mod <- caret::train(
  pp,
  data = df_train %>%
    drop_na(),
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
  tuneGrid = expand.grid( .mtry = 9,
                          .min.node.size = min.node.size_values,
                          .splitrule = "gini"),
  metric = "Accuracy",
  replace = FALSE,
  sample.fraction = 0.5,
  num.trees = 50,
  seed = 42                # for reproducibility
)
print(mod)
# Best min.node.size was 10 with mtry of 9

mod <- caret::train(
  pp,
  data = df_train %>%
    drop_na(),
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
  tuneGrid = expand.grid( .mtry = mtry_values,
                          .min.node.size = 10,
                          .splitrule = "gini"),
  metric = "Accuracy",
  replace = FALSE,
  sample.fraction = 0.5,
  num.trees = 50,
  seed = 42                # for reproducibility
)
print(mod)
# Best mtry value was 12 with min.node.size 10

mod <- caret::train(
  pp,
  data = df_train %>%
    drop_na(),
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
  tuneGrid = expand.grid( .mtry = 12,
                          .min.node.size = 10,
                          .splitrule = splitrule_values),
  metric = "Accuracy",
  replace = FALSE,
  sample.fraction = 0.5,
  num.trees = 50,
  seed = 42                # for reproducibility
)
print(mod)
# Best combination of hyperparameters was mtry value=12, min.node.size=10, and splitrule=gini

mod <- caret::train(
  pp,
  data = df_train |>
    drop_na(),
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
  tuneGrid = expand.grid( .mtry = 12,
                          .min.node.size = 10,
                          .splitrule = "gini"),
  metric = "Accuracy",
  replace = FALSE,
  sample.fraction = 0.5,
  num.trees = 100,
  seed = 42                # for reproducibility
)
print(mod)

# Save relevant data for model testing.
saveRDS(mod,
        here::here("data/rf_mod_for_waterlog.100.rds"))

saveRDS(df_train[, c(target, predictors_selected)],
        here::here("data/cal_mod_for_waterlog.100.rds"))

saveRDS(df_test[, c(target, predictors_selected)],
        here::here("data/val_mod_for_waterlog.100.rds"))


# Evaluation
#### Model Analysis
rf_mod   <- readRDS(here::here("data/rf_mod_for_waterlog.100.rds"))
df_train <- readRDS(here::here("data/cal_bor_for_waterlog.100.rds"))
df_test  <- readRDS(here::here("data/val_bor_for_waterlog.100.rds"))

# Load area to be predicted
raster_mask <- terra::rast(here::here("data-raw/geodata/study_area/area_to_be_mapped.tif"))
# Turn target raster into a dataframe, 1 px = 1 cell
df_mask <- as.data.frame(raster_mask, xy = TRUE)

# Filter only for area of interest
df_mask <- df_mask |>
  dplyr::filter(area_to_be_mapped == 1)

# Display df
head(df_mask) |>
  knitr::kable()

files_covariates <- list.files(
  path = here::here("data-raw/geodata/covariates/"),
  pattern = ".tif$",
  recursive = TRUE,
  full.names = TRUE
)

preds_selected <- names(df_train[, predictors_selected])
files_selected <- files_covariates[apply(sapply(X = preds_selected,
                                                FUN = grepl,
                                                files_covariates),
                                         MARGIN =  1,
                                         FUN = any)]

# Load all rasters as a stack
raster_covariates <- terra::rast(files_selected)

# Get coordinates for which we want data
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
  rf_mod,
  newdata = df_test,
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
mosaicplot(conf_matrix_waterlog_rfbasic$table,
           main = "Confusion matrix")


### Create Prediction Maps
prediction <- predict(
  rf_mod,              # RF model
  newdata = df_predict,
  num.threads = parallel::detectCores() - 1)

# Attach predictions to dataframe and round them
df_predict$prediction <- prediction

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
    na.value = NA,
    option = "viridis",
    name = "Waterlog.100"
  ) +
  ggplot2::theme_classic() +
  ggplot2::scale_x_continuous(expand = c(0, 0)) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::labs(title = "Predicted Waterlog")

# Save raster as .tif file
terra::writeRaster(
  raster_pred,
  "../data/ra_predicted_waterlog100.tif",
  datatype = "FLT4S",  # FLT4S for floats, INT1U for integers (smaller file)
  filetype = "GTiff",  # GeoTiff format
  overwrite = TRUE     # Overwrite existing file
)
