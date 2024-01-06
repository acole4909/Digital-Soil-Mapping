##### Exercise 5: 5.3
#### Data Loading and Preparation
library(caret)
library(recipes)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(terra)
library(here)

# Load Boruta random forest model
rf_bor   <- readRDS(here::here("data/rf_bor_for_waterlog.100.rds"))
df_train <- readRDS(here::here("data/cal_bor_for_waterlog.100.rds"))
df_test  <- readRDS(here::here("data/val_bor_for_waterlog.100.rds"))
predictors_selected <- readRDS(here::here("data/predictors_selected.rds"))

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
        here::here("data/rf_mod_greedy_for_waterlog.100.rds"))

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

set.seed(42)

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
