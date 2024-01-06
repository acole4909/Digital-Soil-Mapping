#### AGDS II: Digital Soil Mapping Exercise 5.1
library(dplyr)
library(ggplot2)
library(tidyverse)
library(terra)
library(here)
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
