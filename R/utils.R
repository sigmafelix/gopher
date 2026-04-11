# Utility functions shared across all GP engines --------------------------

#' Extract coordinate matrix from spatial data
#'
#' @param data An `sf` object or a `data.frame` with coordinate columns.
#' @param coord_cols Character vector of length 2 giving the names of the
#'   coordinate columns (x/longitude first, y/latitude second). Only used
#'   when `data` is not an `sf` object. If `NULL`, the function looks for
#'   columns named `c("x","y")`, `c("lon","lat")`, or `c("longitude","latitude")`.
#'
#' @return A numeric matrix with two columns (`x` and `y`).
#' @keywords internal
extract_coords <- function(data, coord_cols = NULL) {
  if (inherits(data, "sf")) {
    coords <- sf::st_coordinates(data)

    # Non-point geometries can expand to multiple coordinate rows per feature
    # (e.g. POLYGON vertices), but downstream engines expect one location per
    # observation. Use a representative point in those cases.
    if (nrow(coords) != nrow(data)) {
      cli::cli_inform(
        c(
          "Converting non-point {.pkg sf} geometries to representative points for coordinate extraction.",
          "i" = "Using {.fn sf::st_point_on_surface} so each observation contributes exactly one location."
        )
      )
      coords <- sf::st_coordinates(sf::st_point_on_surface(data))
    }

    if (nrow(coords) != nrow(data)) {
      cli::cli_abort(
        "Could not derive exactly one coordinate pair per observation from the {.pkg sf} geometry column."
      )
    }

    # st_coordinates may return X, Y, (Z), (L1, L2) — keep only X and Y
    return(coords[, c("X", "Y"), drop = FALSE])
  }

  # data.frame path: look for coordinate columns
  if (!is.null(coord_cols)) {
    if (length(coord_cols) != 2L) {
      cli::cli_abort("`coord_cols` must be a character vector of length 2.")
    }
    missing_cols <- setdiff(coord_cols, names(data))
    if (length(missing_cols) > 0L) {
      cli::cli_abort(
        "Coordinate column{?s} {.val {missing_cols}} not found in `data`."
      )
    }
    m <- as.matrix(data[, coord_cols])
    colnames(m) <- c("X", "Y")
    return(m)
  }

  # auto-detect
  candidates <- list(
    c("x", "y"),
    c("X", "Y"),
    c("lon", "lat"),
    c("LON", "LAT"),
    c("longitude", "latitude"),
    c("LONGITUDE", "LATITUDE"),
    c("Longitude", "Latitude")
  )
  for (pair in candidates) {
    if (all(pair %in% names(data))) {
      m <- as.matrix(data[, pair])
      colnames(m) <- c("X", "Y")
      return(m)
    }
  }

  cli::cli_abort(
    c(
      "Cannot determine coordinate columns from `data`.",
      "i" = "Provide an `sf` object or a `data.frame` with columns named
             {.code x}/{.code y}, {.code lon}/{.code lat}, or
             {.code longitude}/{.code latitude}.",
      "i" = "You can also pass {.arg coord_cols} to the engine via
             {.fn parsnip::set_engine}."
    )
  )
}

#' Extract coordinate matrix for spatial or spatiotemporal models
#'
#' @param data An `sf` object or a `data.frame` with coordinate columns.
#' @param coord_cols Character vector of length 2 giving spatial coordinate
#'   column names for non-`sf` inputs.
#' @param time_col Optional character scalar naming a time column to append
#'   as a third coordinate dimension (`T`).
#' @param time_scale Numeric scalar used to rescale the time coordinate.
#'
#' @return A numeric matrix with columns `X`, `Y`, and optionally `T`.
#' @keywords internal
extract_st_coords <- function(
    data,
    coord_cols = NULL,
    time_col   = NULL,
    time_scale = 1) {

  coords <- extract_coords(data, coord_cols = coord_cols)
  if (is.null(time_col)) {
    return(coords)
  }

  plain_data <- drop_geometry(data)
  if (!time_col %in% names(plain_data)) {
    cli::cli_abort("Time column {.val {time_col}} was not found in {.arg data}.")
  }

  t_raw <- plain_data[[time_col]]
  t_num <- if (inherits(t_raw, c("POSIXct", "POSIXt"))) {
    as.numeric(t_raw)
  } else if (inherits(t_raw, "Date")) {
    as.numeric(t_raw)
  } else if (is.numeric(t_raw)) {
    as.numeric(t_raw)
  } else {
    t_parsed <- as.POSIXct(t_raw, tz = "UTC")
    if (all(is.na(t_parsed))) {
      cli::cli_abort(
        "Time column {.val {time_col}} must be numeric, Date/POSIXt, or coercible to datetime."
      )
    }
    as.numeric(t_parsed)
  }

  if (length(time_scale) != 1L || !is.numeric(time_scale) ||
      is.na(time_scale) || time_scale == 0) {
    cli::cli_abort("{.arg time_scale} must be a single non-zero numeric value.")
  }

  if (anyNA(t_num)) {
    cli::cli_abort("Missing values detected in {.arg time_col}; remove or impute them first.")
  }

  st_coords <- cbind(coords, T = t_num / time_scale)
  colnames(st_coords) <- c("X", "Y", "T")
  st_coords
}

#' Drop geometry and return a plain data.frame
#'
#' @param data An `sf` object or `data.frame`.
#' @return A `data.frame` (or `tibble`) without a geometry column.
#' @keywords internal
drop_geometry <- function(data) {
  if (inherits(data, "sf")) sf::st_drop_geometry(data) else data
}

# ---- Covariance / variogram model name mapping --------------------------

# Mapping table: canonical gopher name -> engine-specific name
.cov_map <- list(
  gstat = c(
    exponential  = "Exp",
    spherical    = "Sph",
    gaussian     = "Gau",
    matern       = "Mat",
    stein_matern = "Ste",
    power        = "Pow",
    hole_effect  = "Hol",
    linear       = "Lin"
  ),
  fields = c(
    exponential  = "Exponential",
    spherical    = "stationary.cov",
    gaussian     = "Exp.cov",
    matern       = "Matern",
    stein_matern = "Matern"
  ),
  GPvecchia = c(
    exponential  = "exponential",
    spherical    = "spherical",
    gaussian     = "gaussian",
    matern       = "matern",
    stein_matern = "matern15"
  ),
  spNNGP = c(
    exponential  = "exponential",
    spherical    = "spherical",
    gaussian     = "gaussian",
    matern       = "matern",
    stein_matern = "matern"
  ),
  PrestoGP = c(
    exponential  = "matern",
    spherical    = "matern",
    gaussian     = "matern",
    matern       = "matern",
    stein_matern = "matern"
  ),
  sdmTMB = c(
    exponential  = "matern",
    spherical    = "matern",
    gaussian     = "matern",
    matern       = "matern",
    stein_matern = "matern"
  )
)

#' Translate a canonical covariance function name to an engine-specific one
#'
#' @param name A canonical gopher covariance name (e.g. `"exponential"`),
#'   or `NULL` to return the engine default.
#' @param engine One of `"gstat"`, `"fields"`, `"GPvecchia"`, `"spNNGP"`,
#'   or `"PrestoGP"`.
#' @param default The default value to return when `name` is `NULL`.
#'
#' @return A character string with the engine-specific covariance name.
#' @keywords internal
translate_covariance <- function(name, engine, default = NULL) {
  if (is.null(name)) return(default)
  map <- .cov_map[[engine]]
  if (is.null(map)) {
    cli::cli_abort("Unknown engine {.val {engine}} in {.fn translate_covariance}.")
  }
  name_lower <- tolower(name)
  # Allow direct engine-specific names to pass through
  if (name_lower %in% names(map)) {
    return(unname(map[name_lower]))
  }
  # Already an engine-specific name?
  if (name %in% unname(map)) {
    return(name)
  }
  cli::cli_warn(
    "Covariance function {.val {name}} not recognised for engine
     {.val {engine}}. Falling back to {.val {default %||% 'exponential'}}."
  )
  default %||% unname(map["exponential"])
}

# ---- Formula helpers ----------------------------------------------------

#' Parse a model formula into response, covariates, and an intercept flag
#'
#' @param formula A two-sided formula.
#' @param data    A `data.frame` (geometry already dropped).
#'
#' @return A list with elements:
#'   * `response`   – numeric vector of observed values.
#'   * `X`          – numeric model matrix (includes intercept if present).
#'   * `has_covariates` – logical, `TRUE` when predictors other than the
#'     intercept are present.
#'   * `response_name` – character name of the response variable.
#' @keywords internal
parse_formula <- function(formula, data) {
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  response      <- stats::model.response(mf)
  response_name <- as.character(formula[[2]])
  X             <- stats::model.matrix(formula, data = mf)
  # "has_covariates" means there is at least one predictor beyond the intercept
  tt            <- stats::terms(formula)
  preds         <- attr(tt, "term.labels")
  has_cov       <- length(preds) > 0L
  list(
    response      = response,
    X             = X,
    has_covariates = has_cov,
    response_name = response_name
  )
}

# ---- NULL-coalescing operator (in case rlang not attached) ---------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
