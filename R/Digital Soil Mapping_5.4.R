##### 5.4: Probabalistic Predictions
### Data Loading and Preparation
# Load Data
df_full <- readRDS(here::here("data/df_full.rds"))

head(df_full) |>
  knitr::kable()

# Specify target: Waterlog.100
df_train$waterlog.100 <- as.factor(df_train$waterlog.100)
is.factor(df_train$waterlog.100)

# Specify predictors_all: Remove soil sampling and observational data
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

### Model Training
rf_prob <- ranger::ranger(
  probability = TRUE,
  classification = TRUE,
  y = df_train[, "waterlog.100"],     # target variable
  x = df_train[, predictors_selected], # Predictor variables
  seed = 42,                    # Specify the seed for randomization to reproduce the same model again
  num.threads = parallel::detectCores() - 1) # Use all but one CPU core for quick model training

# Print a summary of fitted model
print(rf_prob)

# Save relevant data for model testing in the next chapter.
saveRDS(rf_prob,
        here::here("data/rf_prob_for_waterlog.100.rds"))

saveRDS(df_train[, c(target, predictors_selected)],
        here::here("data/cal_prob_for_waterlog.100.rds"))

saveRDS(df_test[, c(target, predictors_selected)],
        here::here("data/val_prob_for_waterlog.100.rds"))

###Probability Model Analysis
# Load random forest model
rf_prob   <- readRDS(here::here("rf_prob_for_waterlog.100.rds"))
df_train <- readRDS(here::here("cal_prob_for_waterlog.100.rds"))
df_test  <- readRDS(here::here("val_prob_for_waterlog.100.rds"))

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

# Filter that list only for the variables used in the RF
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
  raster_covariates,   # The raster we want to extract from
  df_locations,        # A matrix of x and y values to extract for
  ID = FALSE           # To not add a default ID column to the output
)

df_predict <- cbind(df_locations, df_predict) |>
  tidyr::drop_na()  # Se_TWI2m has a small number of missing data

# Make predictions for validation sites
prediction <- predict(
  rf_prob,           # RF model
  data = df_test,   # Predictor data
  type="response",
  num.threads = parallel::detectCores() - 1
)

#Save predictions to validation data
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

#ROC Curve
library(plotROC)
basicplot <- ggplot(df_test, aes(d=waterlog.100, m=pred[,1])) + geom_roc()
basicplot

library(pROC)
rocobj <- roc(df_test$waterlog.100, predictions)
