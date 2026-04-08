#' Gaussian Process Model for Spatial and Spatiotemporal Data
#'
#' @description
#' `gaussian_process_spatial()` defines a Gaussian Process model for spatial
#' and spatiotemporal data. This model supports multiple backends ("engines"),
#' including **gstat**, **fields**, **GPvecchia**, **spNNGP**, and
#' **PrestoGP**.
#'
#' @param mode A single character string for the prediction outcome mode.
#'   The only possible value for this model is `"regression"`.
#' @param covariance_function The type of covariance function (variogram model)
#'   to use. One of `"exponential"`, `"spherical"`, `"gaussian"`,
#'   `"matern"`, or `"stein_matern"`. Defaults to `NULL` (uses engine
#'   default, typically `"exponential"`).
#' @param range The range (or scale) parameter of the covariance function.
#'   Controls the distance at which spatial correlation effectively vanishes.
#'   Defaults to `NULL` (estimated from data).
#' @param nugget The nugget variance representing micro-scale variation and
#'   measurement error. Defaults to `NULL` (estimated from data).
#' @param sill The partial sill, i.e., the spatially structured variance
#'   component. Defaults to `NULL` (estimated from data).
#'
#' @details
#' ## What does this model do?
#'
#' Gaussian Process (GP) models — commonly known as Kriging in geostatistics —
#' predict values at unobserved locations by leveraging spatial autocorrelation.
#' The model assumes observations are a realisation of a GP with a specified
#' covariance structure. The covariance structure is characterised by a
#' variogram model parameterised by `range`, `nugget`, and `sill`.
#'
#' ## Engines
#'
#' The following engines are available:
#'
#' * `"gstat"` — Uses the **gstat** package for variogram-based kriging
#'   (ordinary, universal, and simple kriging). Supports spatiotemporal
#'   kriging with a `time_col` engine argument.
#' * `"fields"` — Uses the **fields** package (`Krig`/`mKrig`) for spatial
#'   kriging. Supports large datasets via `mKrig`.
#' * `"GPvecchia"` — Uses the **GPvecchia** package for Vecchia-approximated
#'   GP inference, suitable for large spatial and spatiotemporal datasets via
#'   `time_col`.
#' * `"spNNGP"` — Uses the **spNNGP** package for Nearest Neighbor Gaussian
#'   Process models, scalable for large spatial datasets.
#' * `"PrestoGP"` — Uses the **PrestoGP** package for scalable penalized
#'   spatiotemporal Gaussian process models with built-in missing-value
#'   imputation and limit-of-detection handling.
#' * `"sdmTMB"` — Uses the **sdmTMB** package for spatial and spatiotemporal
#'   models via Template Model Builder (TMB) with an INLA-style SPDE mesh.
#'   Supports spatiotemporal modelling via `time_col` and `spatiotemporal`
#'   engine arguments.
#'
#' ## Parameter Mapping
#'
#' Parameters are automatically mapped between the unified gopher interface
#' and engine-specific argument names:
#'
#' | gopher                | gstat (vgm)   | fields (Krig) | GPvecchia  | spNNGP    | PrestoGP | sdmTMB     |
#' |-----------------------|---------------|---------------|------------|-----------|----------|------------|
#' | `covariance_function` | `model`       | `Covariance`  | `covFun`   | `cov.model` | Matérn-only (mapped) | Matérn via SPDE mesh |
#' | `range`               | `range`       | `aRange`      | `range`    | `phi`     | estimated internally | estimated internally |
#' | `nugget`              | `nugget`      | `sigma2`      | `nugget`   | `tau.sq`  | estimated internally | estimated internally |
#' | `sill`                | `psill`       | `sigma2`      | `sigma2`   | `sigma.sq` | estimated internally | estimated internally |
#'
#' ## Spatial Inputs
#'
#' Input data can be provided as:
#' * An `sf` object — geometry column is used for coordinates automatically.
#' * A `data.frame` with columns `x` and `y` (or `lon` and `lat`).
#'
#' ## Covariates (Universal / Residual Kriging)
#'
#' When covariates are included in the formula (e.g., `y ~ x1 + x2`), the
#' model performs **Universal Kriging** — fitting a linear trend model with
#' spatial correlation of residuals. Use `y ~ 1` for **Ordinary Kriging**
#' (no covariates, constant mean).
#'
#' ## Spatiotemporal Kriging
#'
#' Spatiotemporal Gaussian process modelling is supported through the `"gstat"`,
#' `"GPvecchia"`, `"PrestoGP"`, and `"sdmTMB"` engines by passing `time_col` as
#' an engine argument via `set_engine()`. The column specified must contain
#' date/time values (or numeric time indices). For `"sdmTMB"`, also pass
#' `spatiotemporal` (one of `"iid"`, `"ar1"`, `"rw"`) to enable the
#' spatiotemporal random field.
#'
#' @return A `gaussian_process_spatial` model specification of class
#'   `c("gaussian_process_spatial", "model_spec")`.
#'
#' @seealso [parsnip::set_engine()], [parsnip::fit.model_spec()],
#'   [parsnip::predict.model_fit()]
#'
#' @export
#' @examples
#' if (requireNamespace("spacetime", quietly = TRUE) &&
#'     requireNamespace("gstat", quietly = TRUE)) {
#'   data("air", package = "spacetime")
#'
#'   # Convert legacy ST components from `air` to a single-day `sf` table.
#'   day_id <- which.max(colSums(!is.na(air)))
#'   air_day <- data.frame(
#'     station = rownames(air),
#'     pm10 = air[, day_id],
#'     day = dates[day_id],
#'     sp::coordinates(stations)
#'   )
#'   air_day <- air_day[stats::complete.cases(air_day$pm10), ]
#'   air_sf <- sf::st_as_sf(
#'     air_day,
#'     coords = c("coords.x1", "coords.x2"),
#'     crs = 4326,
#'     remove = FALSE
#'   )
#'
#'   n_train <- floor(0.8 * nrow(air_sf))
#'   train_sf <- air_sf[seq_len(n_train), ]
#'   test_sf <- air_sf[seq.int(n_train + 1L, nrow(air_sf)), ]
#'
#'   gp_spec <- gaussian_process_spatial(
#'     covariance_function = "exponential"
#'   ) |>
#'     parsnip::set_engine("gstat")
#'
#'   gp_fit <- parsnip::fit(
#'     gp_spec,
#'     pm10 ~ coords.x1 + coords.x2,
#'     data = train_sf
#'   )
#'
#'   predict(gp_fit, new_data = test_sf, type = "pred_int")
#' }
gaussian_process_spatial <- function(
    mode                = "regression",
    covariance_function = NULL,
    range               = NULL,
    nugget              = NULL,
    sill                = NULL) {

  args <- list(
    covariance_function = rlang::enquo(covariance_function),
    range               = rlang::enquo(range),
    nugget              = rlang::enquo(nugget),
    sill                = rlang::enquo(sill)
  )

  parsnip::new_model_spec(
    "gaussian_process_spatial",
    args               = args,
    eng_args           = NULL,
    mode               = mode,
    user_specified_mode = !missing(mode),
    method             = NULL,
    engine             = NULL
  )
}

# ---- S3 methods -------------------------------------------------------

#' @export
print.gaussian_process_spatial <- function(x, ...) {
  cat("Gaussian Process Model Specification (", x$mode, ")\n\n", sep = "")
  parsnip::model_printer(x, ...)
  invisible(x)
}

#' @export
#' @rdname gaussian_process_spatial
update.gaussian_process_spatial <- function(
    object,
    parameters          = NULL,
    covariance_function = NULL,
    range               = NULL,
    nugget              = NULL,
    sill                = NULL,
    fresh               = FALSE,
    ...) {

  args <- list(
    covariance_function = rlang::enquo(covariance_function),
    range               = rlang::enquo(range),
    nugget              = rlang::enquo(nugget),
    sill                = rlang::enquo(sill)
  )

  parsnip::update_spec(object, parameters, args, fresh, ...)
}
