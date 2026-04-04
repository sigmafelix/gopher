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
- **Spatiotemporal Kriging** — supported via the `gstat` engine with a
  `time_col` engine argument.
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

# Convert one day to an sf table
day_id <- which.max(colSums(!is.na(air)))
air_day <- data.frame(
  station = rownames(air),
  pm10 = air[, day_id],
  day = dates[day_id],
  sp::coordinates(stations)
)
air_day <- air_day[complete.cases(air_day$pm10), ]
air_sf <- st_as_sf(
  air_day,
  coords = c("coords.x1", "coords.x2"),
  crs = 4326,
  remove = FALSE
)

n_train <- floor(0.8 * nrow(air_sf))
train_sf <- air_sf[seq_len(n_train), ]
test_sf <- air_sf[seq.int(n_train + 1L, nrow(air_sf)), ]

# Create a gopher model spec (spNNGP backend)
gp_spec <- gaussian_process_spatial(
  covariance_function = "exponential"
) |>
  set_engine("spNNGP")


# Create a gopher model spec (GPvecchia backend)
gp_spec <- gaussian_process_spatial(
  covariance_function = fields::Matern
) |>
  set_engine("fields")

# Fit to sf training data (universal kriging)
gp_fit <- gp_spec |> fit(pm10 ~ coords.x1 + coords.x2, data = train_sf)

# Predict at new locations
predictions <- predict(gp_fit, new_data = test_sf)

# Prediction intervals
pred_int <- predict(gp_fit, new_data = test_sf, type = "pred_int")
```

## Parameter Mapping

| gopher arg | gstat (vgm) | fields (Krig) | GPvecchia | spNNGP |
|---|---|---|---|---|
| `covariance_function` | `model` | `cov.function` | `covfun.name` | `cov.model` |
| `range` | `range` | `aRange` | `range` | `phi` |
| `nugget` | `nugget` | `lambda × sigma2` | `nugget` | `tau.sq` |
| `sill` | `psill` | `sigma2` | `sigma2` | `sigma.sq` |

## Installation

```r
# Install from GitHub
remotes::install_github("sigmafelix/gopher")

# Install engine packages as needed
install.packages(c("gstat", "fields"))          # most common
install.packages(c("GPvecchia", "spNNGP"))      # scalable engines
```
