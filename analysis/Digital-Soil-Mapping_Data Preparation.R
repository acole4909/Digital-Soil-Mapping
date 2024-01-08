##### AGDS II: Digital Soil Mapping
### Data Preparation
# Load Soil Data
df_obs <- readr::read_csv(
  here::here("data-raw/soildata/berne_soil_sampling_locations.csv")
)
head(df_obs) |>
  knitr::kable()

# Load Covariate Data
list_raster <- list.files(
  here::here("data-raw/geodata/covariates/"),
  full.names = TRUE
)
lapply(
  list_raster,
  function(x) sub(".*/(.*)", "\\1", x)
) |>
  unlist() |>
  head(5) |>
  print()

### Combine Data
all_rasters <- terra::rast(list_raster)
all_rasters

# Extract coordinates from sampling locations
sampling_xy <- df_obs |>
  dplyr::select(x, y)

# From all rasters, extract values for sampling coordinates
df_covars <- terra::extract(
  all_rasters,  # The raster we want to extract from
  sampling_xy,  # A matrix of x and y values to extract for
  ID = FALSE    # To not add a default ID column to the output
)

df_full <- cbind(df_obs, df_covars)
head(df_full) |>
  knitr::kable()

## Data Wrangling
vars_categorical <- df_covars |>

  # Get number of distinct values per variable
  dplyr::summarise(dplyr::across(dplyr::everything(), ~dplyr::n_distinct(.))) |>

  # Turn df into long format for easy filtering
  tidyr::pivot_longer(
    dplyr::everything(),
    names_to = "variable",
    values_to = "n"
  ) |>

  # Filter out variables with 10 or less distinct values
  dplyr::filter(n <= 10) |>

  # Extract the names of these variables
  dplyr::pull('variable')

cat("Variables with less than 10 distinct values:",
    ifelse(length(vars_categorical) == 0, "none", vars_categorical))

df_full <- df_full |>
  dplyr::mutate(dplyr::across(all_of(vars_categorical), ~as.factor(.)))

## Check for Missing Data
# Get number of rows to calculate percentages
n_rows <- nrow(df_full)

# Get number of distinct values per variable
df_full |>
  dplyr::summarise(dplyr::across(dplyr::everything(),
                                 ~ length(.) - sum(is.na(.)))) |>
  tidyr::pivot_longer(dplyr::everything(),
                      names_to = "variable",
                      values_to = "n") |>
  dplyr::mutate(perc_available = round(n / n_rows * 100)) |>
  dplyr::arrange(perc_available) |>
  head(10) |>
  knitr::kable()

# visualize missing data
df_full |>
  dplyr::select(1:20) |>   # reduce data for readability of the plot
  visdat::vis_miss()

if (!dir.exists(here::here("data"))) system(paste0("mkdir ", here::here("data")))
saveRDS(df_full,
        here::here("data/df_full.rds"))
