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
| `"sdmTMB"` | [sdmTMB](https://cran.r-project.org/package=sdmTMB) | SPDE-based Matérn GP via TMB (spatial + spatiotemporal, IID / AR1 / RW fields) |

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
  `PrestoGP`, and `sdmTMB` using a `time_col` engine argument.
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
- `"sdmTMB"` (via `spatiotemporal = "iid"` / `"ar1"` / `"rw"`)

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

### 5) `sdmTMB` spatial-only SPDE-based GP

```r
spec_sdmtmb_spatial <- gaussian_process_spatial(covariance_function = "matern") |>
  set_engine(
    "sdmTMB",
    spatial    = "on",
    mesh_cutoff = 1   # minimum triangle edge length in coordinate units
  )

fit_sdmtmb_spatial <- fit(spec_sdmtmb_spatial, pm10 ~ x + y, data = spatial_train_df)
pred_sdmtmb_spatial <- predict(fit_sdmtmb_spatial, new_data = spatial_test_df)
pi_sdmtmb_spatial   <- predict(fit_sdmtmb_spatial, new_data = spatial_test_df, type = "pred_int")
```

### 6) `sdmTMB` spatiotemporal GP (AR1 temporal structure)

```r
# sdmTMB requires the time column to be an integer or factor
# (train_df / test_df already have day as integer from the shared setup above)

# extra_time: any time values in test data not present in training data must
# be declared at fit time so sdmTMB can pre-allocate the AR1 latent field.
extra_t <- setdiff(
  as.integer(unique(test_df$day)),
  as.integer(unique(train_df$day))
)

spec_sdmtmb_spt <- gaussian_process_spatial(covariance_function = "matern") |>
  set_engine(
    "sdmTMB",
    time_col       = "day",
    spatial        = "on",
    spatiotemporal = "ar1",  # "iid", "ar1", or "rw"
    mesh_cutoff    = 1,
    share_range    = FALSE,
    extra_time     = extra_t  # allow prediction at held-out time steps
  )

fit_sdmtmb_spt  <- fit(spec_sdmtmb_spt, pm10 ~ x + y, data = train_df)
pred_sdmtmb_spt <- predict(fit_sdmtmb_spt, new_data = test_df)
pi_sdmtmb_spt   <- predict(fit_sdmtmb_spt, new_data = test_df, type = "pred_int")
```

> [!NOTE]
> sdmTMB estimates all covariance parameters (range, nugget, sill) via maximum likelihood
> through TMB; `range`, `nugget`, and `sill` arguments in `gaussian_process_spatial()` are
> accepted for interface consistency but are not passed to the engine.

## Parameter Mapping

| gopher arg | gstat (vgm) | fields (Krig) | GPvecchia | spNNGP | PrestoGP | sdmTMB |
|---|---|---|---|---|---|---|
| `covariance_function` | `model` | `cov.function` | `covfun.name` | `cov.model` | Matérn-only (mapped) | Matérn-only via SPDE (mapped) |
| `range` | `range` | `aRange` | `range` | `phi` | estimated internally | estimated internally (ML) |
| `nugget` | `nugget` | `lambda × sigma2` | `nugget` | `tau.sq` | estimated internally | estimated internally (ML) |
| `sill` | `psill` | `sigma2` | `sigma2` | `sigma.sq` | estimated internally | estimated internally (ML) |

### sdmTMB-specific engine arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `coord_cols` | `character(2)` | `NULL` (auto) | Coordinate column names for non-`sf` input |
| `time_col` | `character` | `NULL` | Time column for spatiotemporal models |
| `spatial` | `"on"` / `"off"` | `"on"` | Include a spatial random field |
| `spatiotemporal` | `"off"` / `"iid"` / `"ar1"` / `"rw"` | `"off"` | Temporal structure of the spatiotemporal random field |
| `mesh_cutoff` | `numeric` | auto (bbox diagonal / 10) | Minimum triangle edge length passed to `sdmTMB::make_mesh()` |
| `n_knots` | `integer` | `NULL` | Number of k-means mesh knots (overrides `mesh_cutoff`) |
| `family` | `family` object | `gaussian()` | Response distribution passed to `sdmTMB::sdmTMB()` |
| `extra_time` | `integer` / `numeric` | `NULL` | Additional time values absent from training data that will appear in `newdata`; pre-allocates latent fields for those steps |
| `share_range` | `logical` | `FALSE` | Share spatial range between spatial and spatiotemporal fields |
| `silent` | `logical` | `TRUE` | Suppress fitting messages |

## Installation

```r
# Install from GitHub
remotes::install_github("sigmafelix/gopher")

# Install engine packages as needed
install.packages(c("gstat", "fields"))          # most common
install.packages(c("GPvecchia", "spNNGP"))      # scalable engines
remotes::install_github("NIEHS/PrestoGP")       # PrestoGP engine
install.packages("sdmTMB")                      # sdmTMB engine (TMB/SPDE)
```
