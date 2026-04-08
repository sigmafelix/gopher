# sdmTMB engine for gaussian_process_spatial ----------------------------
#
# Uses the sdmTMB package for spatial and spatiotemporal models via
# Template Model Builder (TMB) with an INLA-style SPDE mesh.
#
# The sdmTMB workflow is:
#   1. Build a spatial mesh:  sdmTMB::make_mesh(data, xy_cols, cutoff/n_knots)
#   2. Fit the model:         sdmTMB::sdmTMB(formula, data, mesh, ...)
#   3. Predict:               predict(fit, newdata = new_data, se_fit = ...)

#' Fit a Gaussian Process model using the sdmTMB engine
#'
#' This function is called internally by parsnip. End users should call
#' [parsnip::fit.model_spec()] on a `gaussian_process_spatial` specification
#' with `set_engine("sdmTMB")`.
#'
#' @param formula A two-sided formula. Use `y ~ 1` for a spatial-only model
#'   without fixed effects (analogous to ordinary kriging), or include
#'   covariates (e.g. `y ~ x1 + x2`) for a model with a fixed-effects trend
#'   (analogous to universal kriging).
#' @param data An `sf` object or a `data.frame` with coordinate columns.
#' @param covariance_function Canonical covariance name (see
#'   [gaussian_process_spatial()]). sdmTMB uses a Matérn covariance via the
#'   SPDE mesh; any value is accepted but mapped to `"matern"`. A warning is
#'   issued for non-Matérn names.
#' @param range    Unused by sdmTMB in the current adapter (estimated
#'   internally via maximum likelihood). Kept for unified interface
#'   compatibility.
#' @param nugget   Unused by sdmTMB in the current adapter. Kept for unified
#'   interface compatibility.
#' @param sill     Unused by sdmTMB in the current adapter. Kept for unified
#'   interface compatibility.
#' @param coord_cols Character vector of length 2 giving the coordinate column
#'   names for non-`sf` input (x/longitude first, y/latitude second). `NULL`
#'   triggers auto-detection.
#' @param time_col Character scalar naming the time column for spatiotemporal
#'   models. `NULL` (default) = purely spatial.
#' @param spatial One of `"on"` (default) or `"off"`. Controls whether a
#'   spatial random field is included.
#' @param spatiotemporal One of `"off"` (default), `"iid"`, `"ar1"`, or
#'   `"rw"`. Controls the temporal structure of the spatiotemporal random
#'   field. Ignored when `time_col` is `NULL`.
#' @param mesh_cutoff Numeric scalar passed to `sdmTMB::make_mesh()` as the
#'   `cutoff` argument (minimum triangle edge length in coordinate units).
#'   `NULL` = automatic (1/10 of the bounding box diagonal).
#' @param n_knots Integer. Number of mesh knots for the k-means mesh type.
#'   If provided, `mesh_cutoff` is ignored and `type = "kmeans"` is used.
#' @param family A `family` object (default `stats::gaussian()`). Passed
#'   directly to `sdmTMB::sdmTMB()`.
#' @param share_range Logical. Whether to share the spatial range between the
#'   spatial and spatiotemporal random fields. Default `FALSE`.
#' @param silent Logical. Suppress model-fitting messages. Default `TRUE`.
#' @param ... Additional arguments forwarded to `sdmTMB::sdmTMB()`.
#'
#' @return A list of class `"gopher_sdmTMB_fit"` containing:
#'   * `sdmtmb_fit` — the fitted `sdmTMB` model object.
#'   * `mesh`       — the `sdmTMBmesh` used for fitting.
#'   * `formula`    — the model formula.
#'   * `xy_cols`    — internal coordinate column names used in the data.
#'   * `coord_cols` — user-supplied `coord_cols` (may be `NULL`).
#'   * `time_col`   — the time column name (may be `NULL`).
#'   * `family`     — the `family` object.
#'
#' @examples
#' if (requireNamespace("sdmTMB", quietly = TRUE)) {
#'   train <- data.frame(
#'     x = runif(50, 0, 10),
#'     y = runif(50, 0, 10)
#'   )
#'   train$z <- sin(train$x) + cos(train$y) + rnorm(50, 0, 0.3)
#'   train_sf <- sf::st_as_sf(train, coords = c("x", "y"), crs = 4326)
#'
#'   fit <- sdmTMB_gp_fit(z ~ 1, data = train_sf)
#'   fit
#' }
#'
#' @export
sdmTMB_gp_fit <- function(
    formula,
    data,
    covariance_function = "matern",
    range               = NULL,
    nugget              = NULL,
    sill                = NULL,
    coord_cols          = NULL,
    time_col            = NULL,
    spatial             = "on",
    spatiotemporal      = "off",
    mesh_cutoff         = NULL,
    n_knots             = NULL,
    family              = stats::gaussian(),
    share_range         = FALSE,
    silent              = TRUE,
    ...) {

  .eval_parsnip_arg <- function(x) {
    if (rlang::is_quosure(x)) {
      return(rlang::eval_tidy(x))
    }
    if (rlang::is_formula(x)) {
      return(rlang::eval_tidy(rlang::f_rhs(x), env = rlang::f_env(x)))
    }
    x
  }

  covariance_function <- .eval_parsnip_arg(covariance_function)
  range               <- .eval_parsnip_arg(range)
  nugget              <- .eval_parsnip_arg(nugget)
  sill                <- .eval_parsnip_arg(sill)
  coord_cols          <- .eval_parsnip_arg(coord_cols)
  time_col            <- .eval_parsnip_arg(time_col)
  spatial             <- .eval_parsnip_arg(spatial)
  spatiotemporal      <- .eval_parsnip_arg(spatiotemporal)
  mesh_cutoff         <- .eval_parsnip_arg(mesh_cutoff)
  n_knots             <- .eval_parsnip_arg(n_knots)
  share_range         <- .eval_parsnip_arg(share_range)
  silent              <- .eval_parsnip_arg(silent)

  rlang::check_installed("sdmTMB", reason = "for the sdmTMB engine")

  spatial        <- rlang::arg_match(spatial,        c("on", "off"))
  spatiotemporal <- rlang::arg_match(spatiotemporal, c("off", "iid", "ar1", "rw"))

  # Warn if a non-Matern covariance was requested
  cov_name <- covariance_function %||% "matern"
  if (!tolower(cov_name) %in% c("matern", "stein_matern")) {
    cli::cli_warn(
      c(
        "sdmTMB uses Matern covariance via the SPDE mesh.",
        "i" = paste0(
          "The requested covariance '", cov_name,
          "' is not directly supported and has been mapped to Matern."
        )
      )
    )
  }

  # ---- Extract coordinates and build a plain data.frame -----------------
  coords     <- extract_coords(data, coord_cols)
  plain_data <- drop_geometry(data)

  # Use safe internal column names to avoid clashing with user columns
  xy_cols <- c(".gopher_X", ".gopher_Y")
  plain_data[[xy_cols[1]]] <- coords[, "X"]
  plain_data[[xy_cols[2]]] <- coords[, "Y"]

  # ---- Build the SPDE mesh ----------------------------------------------
  if (!is.null(n_knots)) {
    mesh_args <- list(
      data    = plain_data,
      xy_cols = xy_cols,
      type    = "kmeans",
      n_knots = as.integer(n_knots)
    )
  } else {
    cutoff <- mesh_cutoff %||% {
      xr <- diff(range(coords[, "X"]))
      yr <- diff(range(coords[, "Y"]))
      sqrt(xr^2 + yr^2) / 10
    }
    mesh_args <- list(
      data    = plain_data,
      xy_cols = xy_cols,
      cutoff  = cutoff
    )
  }

  mesh <- tryCatch(
    do.call(sdmTMB::make_mesh, mesh_args),
    error = function(e) {
      cli::cli_abort(
        c("sdmTMB::make_mesh() failed.", "x" = conditionMessage(e))
      )
    }
  )

  # ---- Fit the sdmTMB model ---------------------------------------------
  fit_args <- list(
    formula        = formula,
    data           = plain_data,
    mesh           = mesh,
    family         = family,
    spatial        = spatial,
    spatiotemporal = spatiotemporal,
    share_range    = share_range,
    silent         = silent
  )

  if (!is.null(time_col)) {
    if (!time_col %in% names(plain_data)) {
      cli::cli_abort(
        "Time column {.val {time_col}} not found in training data."
      )
    }
    fit_args$time <- time_col
  }

  extra <- list(...)
  fit_args <- c(fit_args, extra)

  fit_obj <- tryCatch(
    do.call(sdmTMB::sdmTMB, fit_args),
    error = function(e) {
      cli::cli_abort(
        c("sdmTMB::sdmTMB() failed.", "x" = conditionMessage(e))
      )
    }
  )

  structure(
    list(
      sdmtmb_fit     = fit_obj,
      mesh           = mesh,
      formula        = formula,
      xy_cols        = xy_cols,
      coord_cols     = coord_cols,
      time_col       = time_col,
      family         = family,
      spatiotemporal = spatiotemporal
    ),
    class = "gopher_sdmTMB_fit"
  )
}

#' Predict from an sdmTMB-fitted Gaussian Process model
#'
#' @param object   A `"gopher_sdmTMB_fit"` object returned by
#'   [sdmTMB_gp_fit()].
#' @param new_data An `sf` object or `data.frame` with the same coordinate
#'   structure (and covariate columns, if any) as the training data.
#' @param type     One of `"numeric"` (default, returns `.pred`) or
#'   `"pred_int"` (returns `.pred`, `.pred_lower`, `.pred_upper`).
#' @param level    Confidence level for prediction intervals (default `0.95`).
#' @param coord_cols Character vector of length 2 giving coordinate column
#'   names for non-`sf` new data. `NULL` uses the coord_cols from the fit.
#' @param ... Additional arguments forwarded to `predict.sdmTMB()`.
#'
#' @return A [tibble::tibble()] with prediction columns.
#' @examples
#' if (requireNamespace("sdmTMB", quietly = TRUE)) {
#'   set.seed(1)
#'   train <- data.frame(x = runif(50, 0, 10), y = runif(50, 0, 10))
#'   train$z <- sin(train$x) + cos(train$y) + rnorm(50, 0, 0.3)
#'   train_sf <- sf::st_as_sf(train, coords = c("x", "y"), crs = 4326)
#'
#'   newdat <- sf::st_as_sf(
#'     expand.grid(x = seq(1, 9, l = 4), y = seq(1, 9, l = 4)),
#'     coords = c("x", "y"), crs = 4326
#'   )
#'
#'   fit   <- sdmTMB_gp_fit(z ~ 1, data = train_sf)
#'   sdmTMB_gp_predict(fit, new_data = newdat)
#' }
#'
#' @export
sdmTMB_gp_predict <- function(
    object,
    new_data,
    type       = "numeric",
    level      = 0.95,
    coord_cols = NULL,
    ...) {

  .eval_parsnip_arg <- function(x) {
    if (rlang::is_quosure(x)) {
      return(rlang::eval_tidy(x))
    }
    if (rlang::is_formula(x)) {
      return(rlang::eval_tidy(rlang::f_rhs(x), env = rlang::f_env(x)))
    }
    x
  }

  type       <- .eval_parsnip_arg(type)
  level      <- .eval_parsnip_arg(level)
  coord_cols <- .eval_parsnip_arg(coord_cols)

  rlang::check_installed("sdmTMB", reason = "for the sdmTMB engine")
  type <- rlang::arg_match(type, c("numeric", "pred_int"))

  se_fit <- identical(type, "pred_int")

  # ---- Extract coordinates from new_data --------------------------------
  use_coord_cols <- coord_cols %||% object$coord_cols
  coords_new <- extract_coords(new_data, use_coord_cols)
  plain_new  <- drop_geometry(new_data)

  # Add coordinate columns using the same internal names as in training
  plain_new[[object$xy_cols[1]]] <- coords_new[, "X"]
  plain_new[[object$xy_cols[2]]] <- coords_new[, "Y"]

  # ---- Predict ----------------------------------------------------------
  pred_result <- tryCatch(
    predict(
      object$sdmtmb_fit,
      newdata = plain_new,
      se_fit  = se_fit,
      ...
    ),
    error = function(e) {
      cli::cli_abort(
        c("predict.sdmTMB() failed.", "x" = conditionMessage(e))
      )
    }
  )

  preds <- pred_result$est

  if (type == "numeric") {
    return(tibble::tibble(.pred = preds))
  }

  # pred_int: use delta-method SE from sdmTMB
  se <- pred_result$est_se
  if (is.null(se) || anyNA(se)) {
    cli::cli_abort(
      "sdmTMB did not return valid prediction standard errors for {.val type = 'pred_int'}."
    )
  }

  z <- stats::qnorm((1 + level) / 2)
  tibble::tibble(
    .pred       = preds,
    .pred_lower = preds - z * se,
    .pred_upper = preds + z * se
  )
}
