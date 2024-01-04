##### Exercise 5: 5.4
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

df_test$pred0.5 <- ifelse(df_test$pred[,1]>0.5, 1, 0)
Y1 <- df_test$waterlog.100
Y1 <- as.factor(df_test$waterlog.100)
X1 <- df_test$pred0.5
X1 <- as.factor(df_test$pred0.5)

df_test$pred0.75 <- ifelse(df_test$pred[,1]>0.75, 1, 0)
Y2 <- df_test$waterlog.100
Y2 <- as.factor(df_test$waterlog.100)
X2 <- df_test$pred0.75
X2 <- as.factor(df_test$pred0.75)

df_test$pred0.25 <- ifelse(df_test$pred[,1]>0.25, 1, 0)
Y3 <- df_test$waterlog.100
Y3 <- as.factor(df_test$waterlog.100)
X3 <- df_test$pred0.25
X3 <- as.factor(df_test$pred0.25)

conf_matrix_waterlog_prob0.5 <- caret::confusionMatrix(data=X1, reference=Y1, positive="1")
conf_matrix_waterlog_prob0.5
# Sensitivity (TPR) for Tau=0.5 is 0.7143, FPR for Tau=0.5 is 0.1846
conf_matrix_waterlog_prob0.75 <- caret::confusionMatrix(data=X2, reference=Y2, positive="1")
conf_matrix_waterlog_prob0.75
# Sensitivity (TPR) for Tau=0.75 is 0.2571, FPR for Tau=0.75 is 0.0385
conf_matrix_waterlog_prob0.25 <- caret::confusionMatrix(data=X3, reference=Y3, positive="1")
conf_matrix_waterlog_prob0.25
# Sensitivity (TPR) for Tau=0.25 is 0.8571, FPR for Tau=0.25 is 0.4077

## Plot ROC
x1 <- 0.1846
y1 <- 0.7143
x2 <- 0.0385
y2 <- 0.2571
x3 <- 0.4077
y3 <- 0.8571
df <- data.frame(x = c(x1, x2, x3), y = c(y1, y2, y3))
ggplot(df, aes(x = x, y = y)) +
  geom_point()+
  xlim(0,1)+
  ylim(0,1)

#ROC Curve
basicplot <- ggplot(df_test, aes(d=waterlog.100, m=pred[,1])) + geom_roc()
basicplot

