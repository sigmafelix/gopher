# fields engine for gaussian_process_spatial ---------------------------
#
# Uses fields::Krig() (exact) or fields::mKrig() (large datasets) for
# spatial Gaussian Process regression / kriging.
#
# The fields workflow is:
#   obj <- fields::Krig(x, Y, cov.function = ..., aRange = ..., sigma2 = ...,
#                       lambda = nugget/sigma2, ...)
#   predict(obj, x.pred)
#   predictSE(obj, x.pred)   # prediction standard error

#' Fit a Gaussian Process model using the fields engine
#'
#' This function is called internally by parsnip. End users should call
#' [parsnip::fit.model_spec()] on a `gaussian_process_spatial` specification
#' with `set_engine("fields")`.
#'
#' @param formula  A two-sided formula. Use `y ~ 1` for ordinary kriging;
#'   covariates (e.g. `y ~ x1 + x2`) are passed to `fields::Krig()` as the
#'   `Z` matrix (fixed effects / trend surface).
#' @param data     An `sf` object or a `data.frame` with coordinate columns.
#' @param covariance_function Canonical covariance name. Defaults to
#'   `"exponential"`. Note: `"spherical"` and `"gaussian"` are approximated
#'   via `fields::Exp.cov`. You can also pass a covariance function directly
#'   (e.g. `fields::Matern`), including through a one-sided formula used by
#'   `parsnip` (e.g. `~fields::Matern`).
#' @param range    Range (aRange) parameter for the covariance. `NULL` =
#'   estimated by `fields::Krig()`.
#' @param nugget   Nugget variance. Used to compute `lambda = nugget / sigma2`.
#' @param sill     Partial sill (sigma2). `NULL` = estimated by `Krig`.
#' @param coord_cols Character(2) coordinate column names (non-sf path).
#' @param use_mKrig Logical. Use `fields::mKrig()` instead of `fields::Krig()`
#'   for large datasets. Default `FALSE`.
#' @param ... Additional arguments forwarded to `fields::Krig()` /
#'   `fields::mKrig()`.
#'
#' @return A list of class `"gopher_fields_fit"`.
#' @examples
#' if (requireNamespace("spacetime", quietly = TRUE) &&
#'     requireNamespace("fields", quietly = TRUE)) {
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
#'   fit <- fields_gp_fit(pm10 ~ coords.x1 + coords.x2, data = air_sf)
#'   fit
#' }
#'
#' @export
fields_gp_fit <- function(
    formula,
    data,
    covariance_function = "exponential",
    range               = NULL,
    nugget              = NULL,
    sill                = NULL,
    coord_cols          = NULL,
    use_mKrig           = FALSE,
    ...) {

  rlang::check_installed("fields", reason = "for the fields engine")

  # ---- Extract coordinates and response ---------------------------------
  coords     <- extract_coords(data, coord_cols)
  plain_data <- drop_geometry(data)
  parsed     <- parse_formula(formula, plain_data)

  response <- parsed$response
  has_cov  <- parsed$has_covariates

  # Covariates Z (trend surface) — exclude intercept column
  Z <- if (has_cov) {
    X <- parsed$X
    X[, colnames(X) != "(Intercept)", drop = FALSE]
  } else {
    NULL
  }

  # ---- Resolve covariance function ------------------------------------
  cov_input <- covariance_function
  if (inherits(cov_input, "formula")) {
    cov_input <- rlang::f_rhs(cov_input)
  }

  cov_fn <- NULL
  cov_fn_name <- NULL
  cov_args_auto <- NULL

  if (is.null(cov_input) || length(cov_input) == 0L) {
    cov_fn_name <- "Exponential"
  } else if (is.function(cov_input)) {
    if (identical(cov_input, fields::Matern)) {
      cov_fn_name <- "Matern"
    } else if (identical(cov_input, fields::Exponential)) {
      cov_fn_name <- "Exponential"
    } else {
      cov_fn <- cov_input
      cov_fn_name <- "user_function"
    }
  } else if (is.language(cov_input) || is.symbol(cov_input)) {
    maybe_fn <- tryCatch(
      rlang::eval_tidy(cov_input),
      error = function(e) NULL
    )
    if (is.function(maybe_fn)) {
      if (identical(maybe_fn, fields::Matern)) {
        cov_fn_name <- "Matern"
      } else if (identical(maybe_fn, fields::Exponential)) {
        cov_fn_name <- "Exponential"
      } else {
        cov_fn <- maybe_fn
        cov_fn_name <- "user_function"
      }
    } else {
      cov_label <- paste(deparse(cov_input), collapse = "")
      cov_fn_name <- translate_covariance(
        cov_label, "fields", default = "Exponential"
      )
    }
  } else {
    cov_label <- as.character(cov_input)
    cov_fn_name <- translate_covariance(
      cov_label, "fields", default = "Exponential"
    )
  }

  if (is.null(cov_fn)) {
    # fields::Krig expects a covariance with a C argument.
    # Route Matern/Exponential through stationary.cov.
    cov_fn <- switch(
      cov_fn_name,
      Exponential  = "stationary.cov",
      Matern       = "stationary.cov",
      # fallback to Exp.cov for approximations
      "Exp.cov"
    )
    cov_args_auto <- switch(
      cov_fn_name,
      Exponential = list(Covariance = "Exponential"),
      Matern      = list(Covariance = "Matern"),
      NULL
    )
  }

  # ---- Build Krig/mKrig call arguments ----------------------------------
  krig_args <- list(
    x = coords,
    Y = response
  )

  if (!is.null(Z))     krig_args$Z          <- Z
  if (!is.null(range)) krig_args$aRange     <- range

  # lambda = nugget / sigma2 (noise-to-signal ratio)
  if (!is.null(nugget) && !is.null(sill) && sill > 0) {
    krig_args$lambda <- nugget / sill
  } else if (!is.null(nugget)) {
    krig_args$lambda <- nugget / (stats::var(response, na.rm = TRUE) * 0.9)
  }

  if (!is.null(sill)) krig_args$sigma2 <- sill

  # ---- Fit --------------------------------------------------------------
  fit_fn <- if (use_mKrig) fields::mKrig else fields::Krig

  extra_args <- list(...)
  if (!use_mKrig && !is.null(Z) && is.null(extra_args$m)) {
    # Avoid duplicated low-order trend terms between Krig's null-space basis
    # (default m = 2) and user-provided covariates in Z.
    extra_args$m <- 0
  }
  if (!is.null(cov_args_auto)) {
    if (!is.null(extra_args$cov.args)) {
      extra_args$cov.args <- utils::modifyList(cov_args_auto, extra_args$cov.args)
    } else {
      extra_args$cov.args <- cov_args_auto
    }
  }

  krig_obj <- tryCatch(
    do.call(fit_fn, c(krig_args, list(cov.function = cov_fn), extra_args)),
    error = function(e) {
      cli::cli_abort(
        c(
          "fields::{.fn {if (use_mKrig) 'mKrig' else 'Krig'}} failed.",
          "x" = conditionMessage(e)
        )
      )
    }
  )

  structure(
    list(
      krig_obj    = krig_obj,
      coords      = coords,
      training_data = plain_data,
      formula     = formula,
      has_cov     = has_cov,
      cov_fn_name = cov_fn_name,
      use_mKrig   = use_mKrig
    ),
    class = "gopher_fields_fit"
  )
}

#' Predict from a fields-fitted Gaussian Process model
#'
#' @param object   A `"gopher_fields_fit"` object returned by [fields_gp_fit()].
#' @param new_data An `sf` object or `data.frame` with coordinates.
#' @param type     `"numeric"` (default) or `"pred_int"`.
#' @param level    Confidence level for prediction intervals (default `0.95`).
#' @param coord_cols Character(2) coord column names (non-sf path).
#' @param ... Forwarded to `predict.Krig()` / `predictSE.Krig()`.
#'
#' @return A [tibble::tibble()] with prediction columns.
#' @examples
#' if (requireNamespace("spacetime", quietly = TRUE) &&
#'     requireNamespace("fields", quietly = TRUE)) {
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
#'   fit <- fields_gp_fit(pm10 ~ coords.x1 + coords.x2, data = train_sf)
#'   fields_gp_predict(fit, new_data = test_sf, type = "pred_int")
#' }
#'
#' @export
fields_gp_predict <- function(
    object,
    new_data,
    type       = "numeric",
    level      = 0.95,
    coord_cols = NULL,
    ...) {

  rlang::check_installed("fields", reason = "for the fields engine")
  type <- rlang::arg_match(type, c("numeric", "pred_int"))

  coords_new <- extract_coords(new_data, coord_cols)

  # Covariates for new data
  plain_new <- drop_geometry(new_data)
  Z_new <- if (object$has_cov) {
    tt  <- stats::terms(object$formula)
    preds <- attr(tt, "term.labels")
    # Only include columns that are covariate predictors
    available <- intersect(preds, names(plain_new))
    if (length(available) == 0L) {
      cli::cli_abort(
        "The model was fitted with covariates but {.arg new_data} does not
         contain the required predictor columns."
      )
    }
    as.matrix(plain_new[, available, drop = FALSE])
  } else {
    NULL
  }

  pred_args <- if (object$use_mKrig) {
    l <- list(object$krig_obj, xnew = coords_new)
    if (!is.null(Z_new)) l$Znew <- Z_new
    l
  } else {
    l <- list(object$krig_obj, x = coords_new)
    if (!is.null(Z_new)) l$Z <- Z_new
    l
  }

  preds <- do.call(predict, c(pred_args, list(...)))

  if (type == "numeric") {
    return(tibble::tibble(.pred = as.numeric(preds)))
  }

  # Prediction standard error — available via predictSE
  se <- tryCatch(
    {
      se_args <- if (object$use_mKrig) {
        l <- list(object$krig_obj, xnew = coords_new)
        if (!is.null(Z_new)) l$Znew <- Z_new
        l
      } else {
        l <- list(object$krig_obj, x = coords_new)
        if (!is.null(Z_new)) l$Z <- Z_new
        l
      }
      as.numeric(do.call(fields::predictSE, c(se_args, list(...))))
    },
    error = function(e) {
      cli::cli_warn(
        "Could not compute prediction SE: {conditionMessage(e)}"
      )
      rep(NA_real_, length(preds))
    }
  )

  alpha <- 1 - level
  z     <- stats::qnorm(1 - alpha / 2)
  tibble::tibble(
    .pred       = as.numeric(preds),
    .pred_lower = as.numeric(preds) - z * se,
    .pred_upper = as.numeric(preds) + z * se
  )
}
