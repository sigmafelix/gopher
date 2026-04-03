# gstat engine for gaussian_process_spatial ----------------------------
#
# This file implements the fit and predict bridge functions that parsnip
# calls when the engine is "gstat".
#
# The gstat workflow for kriging is:
#   1. Compute empirical variogram:  gstat::variogram(formula, data)
#   2. Fit variogram model:          gstat::fit.variogram(emp_v, vgm(...))
#   3. Kriging prediction:           gstat::krige(formula, locations,
#                                                 newdata, model)
#
# For spatiotemporal kriging gstat uses gstat::krigeST().

#' Fit a Gaussian Process model using the gstat engine
#'
#' This function is called internally by parsnip. End users should call
#' [parsnip::fit.model_spec()] on a `gaussian_process_spatial` specification
#' with `set_engine("gstat")`.
#'
#' @param formula  A two-sided formula. Use `y ~ 1` for ordinary kriging and
#'   `y ~ x1 + x2` for universal kriging.
#' @param data     An `sf` object or a `data.frame` with coordinate columns.
#' @param covariance_function Canonical covariance name (see
#'   [gaussian_process_spatial()]). Defaults to `"exponential"`.
#' @param range    Initial/fixed range parameter. `NULL` = auto-estimate.
#' @param nugget   Initial/fixed nugget. `NULL` = auto-estimate.
#' @param sill     Initial/fixed partial sill. `NULL` = auto-estimate.
#' @param coord_cols Character(2) column names for coordinates when `data` is
#'   not `sf`. Auto-detected when `NULL`.
#' @param time_col  Character name of the date/time column for spatiotemporal
#'   kriging. `NULL` (default) = purely spatial.
#' @param fit_variogram Logical. When `TRUE` (default) the empirical
#'   variogram is computed and the model is fitted even when `range`,
#'   `nugget`, and `sill` are all provided (so initial values are used as
#'   starting values for optimisation). Set to `FALSE` to use the supplied
#'   parameters without fitting.
#' @param cutoff   Passed to `gstat::variogram()`: maximum distance for
#'   the empirical variogram. `NULL` = gstat default (one-third of bounding
#'   box diagonal).
#' @param width    Passed to `gstat::variogram()`: lag width. `NULL` =
#'   gstat default.
#' @param ... Additional arguments forwarded to `gstat::fit.variogram()`.
#'
#' @return A list of class `"gopher_gstat_fit"` containing the fitted
#'   variogram model, the original training data, and the formula.
#'
#' @export
gstat_gp_fit <- function(
    formula,
    data,
    covariance_function = "exponential",
    range               = NULL,
    nugget              = NULL,
    sill                = NULL,
    coord_cols          = NULL,
    time_col            = NULL,
    fit_variogram       = TRUE,
    cutoff              = NULL,
    width               = NULL,
    ...) {

  rlang::check_installed("gstat", reason = "for the gstat engine")

  is_spatiotemporal <- !is.null(time_col)

  # ---- Convert data to sf if not already --------------------------------
  if (!inherits(data, "sf")) {
    coord_matrix <- extract_coords(data, coord_cols)
    plain_no_xy  <- data[, setdiff(names(data), c(
      "x", "y", "X", "Y", "lon", "lat", "LON", "LAT",
      "longitude", "latitude", "LONGITUDE", "LATITUDE",
      "Longitude", "Latitude"
    )), drop = FALSE]
    data <- sf::st_as_sf(
      cbind(plain_no_xy, as.data.frame(coord_matrix)),
      coords = c("X", "Y"),
      crs    = sf::st_crs(NA)
    )
  }

  # ---- Translate covariance name ----------------------------------------
  vgm_model <- translate_covariance(
    covariance_function, "gstat", default = "Exp"
  )

  # ---- Build initial variogram model ------------------------------------
  # Get response values safely (handles simple symbol LHS only)
  response_vals <- tryCatch(
    sf::st_drop_geometry(data)[[as.character(formula[[2]])]],
    error = function(e) NULL
  )
  if (is.null(response_vals)) {
    # fallback: use model.response
    response_vals <- stats::model.response(
      stats::model.frame(formula, data = sf::st_drop_geometry(data),
                         na.action = stats::na.pass)
    )
  }

  initial_psill  <- sill   %||% stats::var(response_vals, na.rm = TRUE)
  initial_range  <- range  %||% {
    bb     <- sf::st_bbox(data)
    diag_d <- sqrt((bb["xmax"] - bb["xmin"])^2 + (bb["ymax"] - bb["ymin"])^2)
    as.numeric(diag_d) / 3
  }
  initial_nugget <- nugget %||% (initial_psill * 0.1)

  vgm_initial <- gstat::vgm(
    psill  = initial_psill,
    model  = vgm_model,
    range  = initial_range,
    nugget = initial_nugget
  )

  # ---- Empirical variogram and fitting ----------------------------------
  if (is_spatiotemporal) {
    # Spatiotemporal path (gstat::krigeST)
    if (!time_col %in% names(sf::st_drop_geometry(data))) {
      cli::cli_abort(
        "Time column {.val {time_col}} not found in training data."
      )
    }
    fitted_vgm  <- NULL   # deferred: ST variogram fitting is complex
    st_data     <- .make_stfdf(data, time_col, formula)
    fit_obj <- structure(
      list(
        variogram_fit   = vgm_initial,
        training_data   = data,
        st_data         = st_data,
        formula         = formula,
        vgm_model       = vgm_model,
        is_spatiotemporal = TRUE,
        time_col        = time_col
      ),
      class = "gopher_gstat_fit"
    )
    return(fit_obj)
  }

  vgm_call <- list(formula, data)
  if (!is.null(cutoff)) vgm_call$cutoff <- cutoff
  if (!is.null(width))  vgm_call$width  <- width

  emp_variogram <- do.call(gstat::variogram, vgm_call)

  if (fit_variogram) {
    fitted_vgm <- tryCatch(
      gstat::fit.variogram(emp_variogram, model = vgm_initial, ...),
      error = function(e) {
        cli::cli_warn(
          c(
            "Variogram fitting failed: {conditionMessage(e)}",
            "i" = "Falling back to initial variogram parameters."
          )
        )
        vgm_initial
      }
    )
  } else {
    fitted_vgm <- vgm_initial
  }

  structure(
    list(
      variogram_fit     = fitted_vgm,
      empirical_variogram = emp_variogram,
      training_data     = data,
      formula           = formula,
      vgm_model         = vgm_model,
      is_spatiotemporal = FALSE,
      time_col          = NULL
    ),
    class = "gopher_gstat_fit"
  )
}

#' Predict from a gstat-fitted Gaussian Process model
#'
#' @param object  A `"gopher_gstat_fit"` object returned by [gstat_gp_fit()].
#' @param new_data An `sf` object or `data.frame` with the same coordinate
#'   structure as the training data. For spatiotemporal kriging it must also
#'   contain the time column.
#' @param type     One of `"numeric"` (default, returns `.pred`) or
#'   `"pred_int"` (returns `.pred`, `.pred_lower`, `.pred_upper`).
#' @param level    Confidence level for prediction intervals (default `0.95`).
#' @param coord_cols Character(2) coord column names (non-sf path).
#' @param ... Forwarded to `gstat::krige()` / `gstat::krigeST()`.
#'
#' @return A [tibble::tibble()] with prediction columns.
#' @export
gstat_gp_predict <- function(
    object,
    new_data,
    type       = "numeric",
    level      = 0.95,
    coord_cols = NULL,
    ...) {

  rlang::check_installed("gstat", reason = "for the gstat engine")
  type <- rlang::arg_match(type, c("numeric", "pred_int"))

  train_data <- object$training_data
  formula    <- object$formula
  vgm_fit    <- object$variogram_fit

  # ---- Convert new_data to sf -------------------------------------------
  if (!inherits(new_data, "sf")) {
    coords       <- extract_coords(new_data, coord_cols)
    plain_no_xy  <- new_data[, setdiff(names(new_data), c(
      "x", "y", "X", "Y", "lon", "lat", "LON", "LAT",
      "longitude", "latitude", "LONGITUDE", "LATITUDE",
      "Longitude", "Latitude"
    )), drop = FALSE]
    new_data <- sf::st_as_sf(
      cbind(plain_no_xy, as.data.frame(coords)),
      coords = c("X", "Y"),
      crs    = sf::st_crs(train_data)
    )
  } else {
    if (!is.na(sf::st_crs(train_data)) &&
        sf::st_crs(new_data) != sf::st_crs(train_data)) {
      new_data <- sf::st_transform(new_data, sf::st_crs(train_data))
    }
  }

  # ---- Kriging ----------------------------------------------------------
  if (object$is_spatiotemporal) {
    # Spatiotemporal prediction via krigeST
    time_col <- object$time_col
    st_new   <- .make_stfdf(new_data, time_col, formula)
    krige_result <- gstat::krigeST(
      formula,
      data    = object$st_data,
      newdata = st_new,
      modelList = vgm_fit,
      ...
    )
    preds    <- krige_result@data$var1.pred
    variance <- krige_result@data$var1.var
  } else {
    krige_result <- gstat::krige(
      formula,
      locations = train_data,
      newdata   = new_data,
      model     = vgm_fit,
      debug.level = 0,
      ...
    )
    preds    <- krige_result$var1.pred
    variance <- krige_result$var1.var
  }

  # ---- Build output tibble ----------------------------------------------
  if (type == "numeric") {
    return(tibble::tibble(.pred = preds))
  }

  # pred_int: Normal-approximation interval
  alpha <- 1 - level
  z     <- stats::qnorm(1 - alpha / 2)
  se    <- sqrt(pmax(variance, 0))
  tibble::tibble(
    .pred       = preds,
    .pred_lower = preds - z * se,
    .pred_upper = preds + z * se
  )
}

# ---- Internal helper: build STFDF for spacetime-based gstat kriging ---

#' @keywords internal
.make_stfdf <- function(data, time_col, formula) {
  rlang::check_installed(
    c("spacetime", "sp"),
    reason = "for spatiotemporal kriging with gstat"
  )
  coords  <- extract_coords(data)
  sp_obj  <- sp::SpatialPoints(coords)
  times   <- sf::st_drop_geometry(data)[[time_col]]
  df_data <- sf::st_drop_geometry(data)
  df_data[[time_col]] <- NULL
  spacetime::STFDF(sp_obj, times, df_data)
}
