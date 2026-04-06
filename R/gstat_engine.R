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
#' @examples
#' if (requireNamespace("spacetime", quietly = TRUE)) {
#'   data("air", package = "spacetime")
#'
#'   # `air` is provided as legacy ST components (matrix + SpatialPoints + Date).
#'   # Convert one day to an `sf` object for spatial GP fitting.
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
#'   fit <- gstat_gp_fit(pm10 ~ 1, data = air_sf)
#'   fit
#' }
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

  # Save original data before sf conversion so the spatiotemporal path can
  # build STIDF with all original columns (including any coord-named columns
  # that double as formula covariates, e.g. `pm10 ~ x + y`).
  orig_data <- data

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
    rlang::check_installed(
      c("spacetime", "sp"),
      reason = "for spatiotemporal kriging with gstat"
    )
    if (!time_col %in% names(sf::st_drop_geometry(data))) {
      cli::cli_abort(
        "Time column {.val {time_col}} not found in training data."
      )
    }

    st_data <- .make_stidf(orig_data, time_col)

    # ---- Compute temporal lag structure from data -----------------------
    plain_train <- drop_geometry(orig_data)
    t_raw_train <- plain_train[[time_col]]
    t_posix_all <- if (inherits(t_raw_train, "POSIXct")) {
      t_raw_train
    } else if (inherits(t_raw_train, c("Date", "POSIXt"))) {
      as.POSIXct(as.character(t_raw_train), tz = "UTC")
    } else if (is.numeric(t_raw_train)) {
      as.POSIXct(as.Date(t_raw_train, origin = "1970-01-01"), tz = "UTC")
    } else {
      suppressWarnings(as.POSIXct(t_raw_train, tz = "UTC"))
    }
    t_uniq  <- sort(unique(t_posix_all))
    t_diffs <- as.numeric(diff(t_uniq), units = "secs")
    t_lag   <- if (length(t_diffs) > 0L) min(t_diffs[t_diffs > 0]) else 86400
    n_t_lags <- min(max(length(t_uniq) - 1L, 1L), 5L)
    tlags    <- seq(0, by = t_lag, length.out = n_t_lags + 1L)

    # ---- Spatiotemporal empirical variogram ----------------------------
    # variogramST.STIDF is broken when multiple stations share a time point;
    # try to use STFDF (requires a regular grid) for the variogram step.
    vgm_data <- .try_make_stfdf(orig_data, time_col) %||% st_data
    emp_vgm_st <- tryCatch(
      gstat::variogramST(formula, data = vgm_data, tlags = tlags,
                         cutoff = initial_range),
      error = function(e) {
        cli::cli_warn(
          c(
            "ST empirical variogram failed: {conditionMessage(e)}",
            "i" = "Falling back to initial ST variogram parameters."
          )
        )
        NULL
      }
    )

    # ---- Build and fit a separable ST variogram model ------------------
    vgm_sp   <- gstat::vgm(initial_psill * 0.9, vgm_model, initial_range,
                            nugget = initial_nugget)
    vgm_t    <- gstat::vgm(1, "Exp", t_lag * 2)
    vgmst_init <- gstat::vgmST("separable", space = vgm_sp, time = vgm_t,
                                sill = initial_psill)

    fitted_vgm_st <- if (!is.null(emp_vgm_st)) {
      tryCatch(
        gstat::fit.StVariogram(emp_vgm_st, vgmst_init),
        error = function(e) {
          cli::cli_warn(
            c(
              "ST variogram fitting failed: {conditionMessage(e)}",
              "i" = "Using initial ST variogram parameters."
            )
          )
          vgmst_init
        }
      )
    } else {
      vgmst_init
    }

    fit_obj <- structure(
      list(
        variogram_fit     = fitted_vgm_st,
        training_data     = data,
        st_data           = st_data,
        formula           = formula,
        vgm_model         = vgm_model,
        is_spatiotemporal = TRUE,
        time_col          = time_col
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
#' @examples
#' if (requireNamespace("spacetime", quietly = TRUE)) {
#'   data("air", package = "spacetime")
#'
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
#'   fit <- gstat_gp_fit(pm10 ~ coords.x1 + coords.x2, data = train_sf)
#'   gstat_gp_predict(fit, new_data = test_sf, type = "pred_int")
#' }
#'
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

  # Save original new_data before sf conversion: the spatiotemporal path
  # builds STIDF from the original (so formula covariates named after
  # coordinate columns, e.g. `pm10 ~ x + y`, remain in the data slot).
  orig_new_data <- new_data

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
    # Spatiotemporal prediction via krigeST.
    # computeVar = TRUE is needed to obtain prediction variances for pred_int.
    time_col <- object$time_col
    st_new   <- .make_stidf(orig_new_data, time_col)
    krige_result <- gstat::krigeST(
      formula,
      data        = object$st_data,
      newdata     = st_new,
      modelList   = vgm_fit,
      computeVar  = (type == "pred_int"),
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

# ---- Internal helper: build STIDF for spacetime-based gstat kriging ---

#' Build a spacetime::STIDF from a data frame or sf object
#'
#' Creates an irregular space-time data frame (STIDF) suitable for
#' `gstat::krigeST()`. One row per observation; no full regular grid required.
#' The time column is coerced to POSIXct if it is a Date or numeric value.
#'
#' @param data An `sf` object or plain `data.frame`. Coordinate columns are
#'   taken from the `sf` geometry or detected via [extract_coords()].
#' @param time_col Character scalar naming the time column in `data`.
#'
#' @return A `spacetime::STIDF` object.
#' @keywords internal
.make_stidf <- function(data, time_col) {
  rlang::check_installed(
    c("spacetime", "sp"),
    reason = "for spatiotemporal kriging with gstat"
  )
  coords  <- extract_coords(data)
  sp_obj  <- sp::SpatialPoints(coords)
  plain   <- drop_geometry(data)

  t_raw <- plain[[time_col]]
  t_posix <- if (inherits(t_raw, "POSIXct")) {
    t_raw
  } else if (inherits(t_raw, c("Date", "POSIXt"))) {
    as.POSIXct(as.character(t_raw), tz = "UTC")
  } else if (is.numeric(t_raw)) {
    as.POSIXct(as.Date(t_raw, origin = "1970-01-01"), tz = "UTC")
  } else {
    parsed <- suppressWarnings(as.POSIXct(t_raw, tz = "UTC"))
    if (all(is.na(parsed))) {
      cli::cli_abort(
        "Cannot coerce time column {.val {time_col}} to POSIXct for STIDF construction."
      )
    }
    parsed
  }

  df_data <- plain
  df_data[[time_col]] <- NULL
  spacetime::STIDF(sp_obj, t_posix, df_data)
}

# ---- Internal helper: build STFDF for variogramST -----------------------
#
# variogramST.STIDF has a subscript-out-of-bounds bug when multiple spatial
# observations share the same time point (the common case: n_stations per
# time step). variogramST works correctly on STFDF (full regular grid).
#
# This helper tries to reshape long-format data into STFDF.  It returns NULL
# if the data is not a complete regular grid (so the caller can skip the
# empirical variogram and fall back to the initial model).
#
#' @keywords internal
.try_make_stfdf <- function(data, time_col) {
  plain   <- drop_geometry(data)
  coords  <- extract_coords(data)

  t_raw <- plain[[time_col]]
  t_posix <- if (inherits(t_raw, "POSIXct")) {
    t_raw
  } else if (inherits(t_raw, c("Date", "POSIXt"))) {
    as.POSIXct(as.character(t_raw), tz = "UTC")
  } else if (is.numeric(t_raw)) {
    as.POSIXct(as.Date(t_raw, origin = "1970-01-01"), tz = "UTC")
  } else {
    suppressWarnings(as.POSIXct(t_raw, tz = "UTC"))
  }

  times_uniq <- sort(unique(t_posix))
  n_t <- length(times_uniq)

  # Round-trip coordinate matrix to a character key for deduplication
  coord_key <- paste(round(coords[, 1], 10), round(coords[, 2], 10), sep = "_")
  locs_uniq_key <- unique(coord_key)
  n_s <- length(locs_uniq_key)

  # Only proceed when the data is a complete regular (stations × times) grid
  if (nrow(plain) != n_s * n_t) {
    return(NULL)
  }

  # Determine time-index and space-index for each row
  t_idx  <- match(t_posix,   times_uniq)
  sp_idx <- match(coord_key, locs_uniq_key)

  # Sort to time-major, space-inner order (STFDF layout)
  ord  <- order(t_idx, sp_idx)
  df_data <- plain[ord, ]
  df_data[[time_col]] <- NULL

  sp_uniq  <- sp::SpatialPoints(coords[match(locs_uniq_key, coord_key), ,
                                       drop = FALSE])

  tryCatch(
    spacetime::STFDF(sp_uniq, times_uniq, df_data),
    error = function(e) NULL
  )
}
