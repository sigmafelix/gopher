# gopher

<!-- badges: start -->
<!-- badges: end -->

`gopher` is a tidymodels-compatible **Gaussian Process model zoo** for R. It
provides a unified [parsnip](https://parsnip.tidymodels.org/) interface to
multiple GP/Kriging backends available on CRAN.

## Engines

| Engine | Package | Description |
|---|---|---|
| `"gstat"` | [gstat](https://cran.r-project.org/package=gstat) | Variogram-based kriging (ordinary, universal, spatiotemporal) |
| `"fields"` | [fields](https://cran.r-project.org/package=fields) | Krig / mKrig exact GP regression |
| `"GPvecchia"` | [GPvecchia](https://cran.r-project.org/package=GPvecchia) | Vecchia-approximated GP (scalable) |
| `"spNNGP"` | [spNNGP](https://cran.r-project.org/package=spNNGP) | Nearest-Neighbor GP (scalable, MCMC) |
| `"PrestoGP"` | [PrestoGP](https://github.com/NIEHS/PrestoGP) | Scalable penalized spatiotemporal GP with LOD-aware missing-value imputation |

## Features

- **tidymodels-native** — works with `parsnip::fit()`, `predict()`, and
  `workflows`.
- **Automatic parameter mapping** — Kriging-centred argument names
  (`range`, `nugget`, `sill`, `covariance_function`) are automatically
  translated to each engine's native names.
- **`sf` input support** — pass `sf` objects directly; coordinates are
  extracted from the geometry column automatically.
- **Covariates / Universal Kriging** — use a formula like `y ~ x1 + x2` to
  include fixed-effect trend terms (residual GP).
- **Spatiotemporal Kriging / GP** — supported via `gstat`, `GPvecchia`,
  and `PrestoGP` using a `time_col` engine argument.
- **Prediction intervals** — use `predict(fit, type = "pred_int")`.
- **Hyperparameter tuning** — `dials` parameter functions (`covariance_function()`,
  `gp_range()`, `gp_nugget()`, `gp_sill()`) integrate with `tune`.

## Quick Start

```r
library(gopher)
library(parsnip)
library(sf)
library(spacetime)

# Real-world PM10 data used in classic gstat/spacetime examples.
# `air` is delivered as legacy ST components: matrix + stations + dates.
data("air", package = "spacetime")

# Create a compact spatiotemporal table (station x multiple days)
day_ids <- head(which(colSums(!is.na(air)) > 0), 5)
air_st <- do.call(
  rbind,
  lapply(day_ids, function(i) {
    data.frame(
      station = rownames(air),
      pm10 = air[, i],
      day = dates[i],
      sp::coordinates(stations)
    )
  })
)
air_st <- air_st[complete.cases(air_st$pm10), ]
air_sf <- st_as_sf(
  air_st,
  coords = c("coords.x1", "coords.x2"),
  crs = 4326,
  remove = FALSE
)

# Hold out both stations and timestamps for a genuine spatiotemporal test set
stations_all <- sort(unique(air_sf$station))
days_all <- sort(unique(air_sf$day))

stations_test <- tail(stations_all, max(1, floor(length(stations_all) * 0.2)))
days_test <- tail(days_all, max(1, floor(length(days_all) * 0.2)))

train_sf <- air_sf[
  !(air_sf$station %in% stations_test) & !(air_sf$day %in% days_test),
]
test_sf <- air_sf[
  air_sf$station %in% stations_test & air_sf$day %in% days_test,
]

# Spatiotemporal model spec (gstat backend)
gp_spec <- gaussian_process_spatial(covariance_function = "exponential") |>
  set_engine("gstat", time_col = "day")

# Fit to sf training data (universal kriging)
gp_fit <- gp_spec |> fit(pm10 ~ coords.x1 + coords.x2, data = train_sf)

# Predict at new locations
predictions <- predict(gp_fit, new_data = test_sf)

# Prediction intervals
pred_int <- predict(gp_fit, new_data = test_sf, type = "pred_int")
```

## Spatiotemporal Prediction Examples (All ST-Capable Engines)

`gopher` currently supports spatiotemporal prediction for:
- `"gstat"`
- `"GPvecchia"`
- `"PrestoGP"`

`"fields"` and the current `spNNGP` adapter are spatial-only in `gopher`.
The shared `train_df` / `test_df` split below is disjoint in both station and
timestamp.

```r
library(gopher)
library(parsnip)
library(sf)
library(spacetime)

data("air", package = "spacetime")

# Build station x day table, then convert to sf
day_ids <- head(which(colSums(!is.na(air)) > 0), 1000)
air_st <- do.call(
  rbind,
  lapply(day_ids, function(i) {
    data.frame(
      station = rownames(air),
      pm10 = air[, i],
      day = dates[i],
      sp::coordinates(stations)
    )
  })
)
air_st <- air_st[complete.cases(air_st$pm10), ]
air_sf <- st_as_sf(
  air_st,
  coords = c("coords.x1", "coords.x2"),
  crs = 4326,
  remove = FALSE
)

# Spatiotemporal split: training and prediction use different stations and days
stations_all <- sort(unique(air_sf$station))
days_all <- sort(unique(air_sf$day))

stations_test <- tail(stations_all, max(2, floor(length(stations_all) * 0.2)))
days_test <- tail(days_all, max(5, floor(length(days_all) * 0.2)))

train_sf <- air_sf[
  !(air_sf$station %in% stations_test) & !(air_sf$day %in% days_test),
]
test_sf <- air_sf[
  air_sf$station %in% stations_test & air_sf$day %in% days_test,
]

train_df <- sf::st_drop_geometry(train_sf) |>
  dplyr::rename(x = coords.x1, y = coords.x2) |>
  dplyr::select(-station) |>
  dplyr::mutate(day = as.integer(day))
test_df <- sf::st_drop_geometry(test_sf) |>
  dplyr::rename(x = coords.x1, y = coords.x2) |>
  dplyr::select(-station) |>
  dplyr::mutate(day = as.integer(day))

train_df <- train_df[stats::complete.cases(train_df[, c("pm10", "x", "y", "day")]), , drop = FALSE]
test_df <- test_df[stats::complete.cases(test_df[, c("x", "y", "day")]), , drop = FALSE]

# For the spatial-only spNNGP example, use one day with held-out stations.
spatial_day <- min(days_all)
spatial_sf <- air_sf[air_sf$day == spatial_day, ]
spatial_train_sf <- spatial_sf[!(spatial_sf$station %in% stations_test), ]
spatial_test_sf <- spatial_sf[spatial_sf$station %in% stations_test, ]

spatial_train_df <- sf::st_drop_geometry(spatial_train_sf) |>
  dplyr::rename(x = coords.x1, y = coords.x2) |>
  dplyr::select(-station) |>
  dplyr::mutate(day = as.integer(day))
spatial_test_df <- sf::st_drop_geometry(spatial_test_sf) |>
  dplyr::rename(x = coords.x1, y = coords.x2) |>
  dplyr::select(-station) |>
  dplyr::mutate(day = as.integer(day))

```

### 1) `gstat` spatiotemporal kriging

```r
spec_gstat <- gaussian_process_spatial(covariance_function = "exponential") |>
  set_engine("gstat", time_col = "day")

fit_gstat <- fit(spec_gstat, pm10 ~ x + y, data = train_df)
pred_gstat <- predict(fit_gstat, new_data = test_df)
pi_gstat <- predict(fit_gstat, new_data = test_df, type = "pred_int")
```
> [!CAUTION]
> `gstat` version utilizes `STIDF` class objects to represent spatiotemporal data, which may consume long time to run.



### 2) `GPvecchia` scalable spatiotemporal GP

```r
spec_gpvecchia <- gaussian_process_spatial(covariance_function = "matern") |>
  set_engine(
    "GPvecchia",
    time_col = "day",   # use time as 3rd coordinate dimension
    time_scale = 1, # seconds -> days scale
    m = 50
  )

fit_gpvecchia <- fit(spec_gpvecchia, pm10 ~ x + y, data = train_df)

pred_gpvecchia <- predict(fit_gpvecchia, new_data = test_df)
pi_gpvecchia <- predict(fit_gpvecchia, new_data = test_df, type = "pred_int")
```

### 3) `spNNGP` spatial-only NNGP

```r
spec_spnngp <- gaussian_process_spatial(covariance_function = "exponential") |>
  set_engine(
    "spNNGP",
    n_neighbors = 12,
    n_samples = 1000,
    n_burnin = 500
  )

fit_spnngp <- fit(spec_spnngp, pm10 ~ x + y, data = spatial_train_df)
pred_spnngp <- predict(fit_spnngp, new_data = spatial_test_df)
pi_spnngp <- predict(fit_spnngp, new_data = spatial_test_df, type = "pred_int")
```

### 4) `PrestoGP` scalable spatiotemporal GP with LOD-aware imputation

```r
# Build data frames explicitly and ensure required columns are NA-free.
# (`impute_y` imputes missing outcomes, but predictors/coords/time must be complete.)

# Example LOD threshold: bottom decile in the training outcome
lod_upper <- as.numeric(stats::quantile(train_df$pm10, 0.10, na.rm = TRUE))

# Avoid "m >= n" warning on small training sets
n_neighbors_presto <- max(3L, min(15L, nrow(train_df) - 1L))

spec_prestogp <- gaussian_process_spatial(covariance_function = "matern") |>
  set_engine(
    "PrestoGP",
    time_col = "day",      # use time as 3rd coordinate dimension
    time_scale = 86400,    # seconds -> days scale
    n_neighbors = n_neighbors_presto,
    impute_y = TRUE,
    lod_upper = lod_upper,
    penalty = "lasso",
    quiet = TRUE
  )

fit_prestogp <- fit(spec_prestogp, pm10 ~ x + y, data = train_df)
pred_prestogp <- predict(fit_prestogp, new_data = test_df)
pi_prestogp <- predict(fit_prestogp, new_data = test_df, type = "pred_int")
```

## Parameter Mapping

| gopher arg | gstat (vgm) | fields (Krig) | GPvecchia | spNNGP | PrestoGP |
|---|---|---|---|---|---|
| `covariance_function` | `model` | `cov.function` | `covfun.name` | `cov.model` | Matérn-only (mapped) |
| `range` | `range` | `aRange` | `range` | `phi` | estimated internally |
| `nugget` | `nugget` | `lambda × sigma2` | `nugget` | `tau.sq` | estimated internally |
| `sill` | `psill` | `sigma2` | `sigma2` | `sigma.sq` | estimated internally |

## Installation

```r
# Install from GitHub
remotes::install_github("sigmafelix/gopher")

# Install engine packages as needed
install.packages(c("gstat", "fields"))          # most common
install.packages(c("GPvecchia", "spNNGP"))      # scalable engines
remotes::install_github("NIEHS/PrestoGP")       # PrestoGP engine
```
