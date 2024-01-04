##### Exercise 5: 5.4
# Load data
rf_bor   <- readRDS(here::here("data/rf_bor_for_waterlog.100.rds"))
df_train <- readRDS(here::here("data/cal_bor_for_waterlog.100.rds"))
df_test  <- readRDS(here::here("data/val_bor_for_waterlog.100.rds"))
predictors_selected <- readRDS(here::here("data/predictors_selected.rds"))
library(plotROC)

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

random_files <- sample(files_covariates, 2)
terra::rast(random_files[1])
terra::rast(random_files[2])

# Load all rasters as a stack
raster_covariates <- terra::rast(random_files)

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
  rf_prob,           # RF model
  data = df_test,   # Predictor data
  type="response",
  num.threads = parallel::detectCores() - 1
)

# Save predictions to validation df
df_test$pred <- prediction$predictions
predictions <- prediction$predictions

#ROC Curve
basicplot <- ggplot(df_test, aes(d=waterlog.100, m=pred[,1])) + geom_roc()
basicplot

