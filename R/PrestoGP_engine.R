# PrestoGP engine for gaussian_process_spatial ----------------------------
#
# Uses the PrestoGP package's S4 API for scalable penalized spatiotemporal
# Gaussian process regression with optional missing-value imputation and
# limit-of-detection (LOD) constraints.

#' Fit a Gaussian Process model using the PrestoGP engine
#'
#' This function is called internally by parsnip. End users should call
#' [parsnip::fit.model_spec()] on a `gaussian_process_spatial` specification
#' with `set_engine("PrestoGP")`.
#'
#' @param formula A two-sided formula.
#' @param data An `sf` object or `data.frame` with coordinate columns.
#' @param covariance_function Canonical covariance name. PrestoGP currently
#'   uses Matérn covariance internally; non-Matérn values are accepted but
#'   mapped to Matérn.
#' @param range Unused by PrestoGP in the current adapter (kept for unified
#'   interface compatibility).
#' @param nugget Unused by PrestoGP in the current adapter (kept for unified
#'   interface compatibility).
#' @param sill Unused by PrestoGP in the current adapter (kept for unified
#'   interface compatibility).
#' @param n_neighbors Number of neighbors for Vecchia approximation.
#' @param model_type `"vecchia"` (default) or `"full"`.
#' @param coord_cols Character(2) coordinate column names for non-sf input.
#' @param time_col Optional character scalar specifying a time column for
#'   spatiotemporal modelling.
#' @param time_scale Numeric scalar used to rescale time when `time_col` is
#'   provided. Default `1`.
#' @param scaling Optional integer vector passed to `PrestoGP::prestogp_fit()`
#'   to group location dimensions by scale parameter.
#' @param common_scale Optional logical passed to `PrestoGP::prestogp_fit()`.
#' @param impute_y Logical or `NULL`. If `NULL`, it is automatically enabled
#'   when missing outcomes or LOD bounds are present.
#' @param lod_upper Upper LOD bound(s) passed to `lod.upper`.
#' @param lod_lower Lower LOD bound(s) passed to `lod.lower`.
#' @param n_impute Number of multiple imputations for LOD-aware missingness.
#' @param eps_impute Convergence tolerance for multiple imputation.
#' @param maxit_impute Maximum iterations for multiple imputation.
#' @param penalty Penalty type passed to PrestoGP (`"lasso"`, `"relaxed"`,
#'   `"MCP"`, `"SCAD"`).
#' @param alpha Elastic-net mixing parameter passed to PrestoGP.
#' @param family GLM family passed to PrestoGP (`"gaussian"` or `"binomial"`).
#' @param quiet Logical; suppress iterative fit messages if `TRUE`.
#' @param ... Additional arguments forwarded to `PrestoGP::prestogp_fit()`.
#'
#' @return A list of class `"gopher_PrestoGP_fit"`.
#' @export
PrestoGP_gp_fit <- function(
    formula,
    data,
    covariance_function = "matern",
    range               = NULL,
    nugget              = NULL,
    sill                = NULL,
    n_neighbors         = 15L,
    model_type          = c("vecchia", "full"),
    coord_cols          = NULL,
    time_col            = NULL,
    time_scale          = 1,
    scaling             = NULL,
    common_scale        = NULL,
    impute_y            = NULL,
    lod_upper           = NULL,
    lod_lower           = NULL,
    n_impute            = 10L,
    eps_impute          = 0.01,
    maxit_impute        = 0L,
    penalty             = c("lasso", "relaxed", "MCP", "SCAD"),
    alpha               = 1,
    family              = c("gaussian", "binomial"),
    quiet               = TRUE,
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

  .as_flag <- function(x, name) {
    if (is.null(x)) return(NULL)
    if (is.logical(x) && length(x) == 1L && !is.na(x)) return(x)
    if (is.numeric(x) && length(x) == 1L && !is.na(x)) {
      return(as.logical(as.integer(x)))
    }
    if (is.character(x) && length(x) == 1L) {
      low <- tolower(x)
      if (low %in% c("true", "t", "1")) return(TRUE)
      if (low %in% c("false", "f", "0")) return(FALSE)
    }
    cli::cli_abort("{.arg {name}} must be a single TRUE/FALSE value.")
  }

  covariance_function <- .eval_parsnip_arg(covariance_function)
  range               <- .eval_parsnip_arg(range)
  nugget              <- .eval_parsnip_arg(nugget)
  sill                <- .eval_parsnip_arg(sill)
  n_neighbors         <- .eval_parsnip_arg(n_neighbors)
  model_type          <- .eval_parsnip_arg(model_type)
  coord_cols          <- .eval_parsnip_arg(coord_cols)
  time_col            <- .eval_parsnip_arg(time_col)
  time_scale          <- .eval_parsnip_arg(time_scale)
  scaling             <- .eval_parsnip_arg(scaling)
  common_scale        <- .eval_parsnip_arg(common_scale)
  impute_y            <- .eval_parsnip_arg(impute_y)
  lod_upper           <- .eval_parsnip_arg(lod_upper)
  lod_lower           <- .eval_parsnip_arg(lod_lower)
  n_impute            <- .eval_parsnip_arg(n_impute)
  eps_impute          <- .eval_parsnip_arg(eps_impute)
  maxit_impute        <- .eval_parsnip_arg(maxit_impute)
  penalty             <- .eval_parsnip_arg(penalty)
  alpha               <- .eval_parsnip_arg(alpha)
  family              <- .eval_parsnip_arg(family)
  quiet               <- .eval_parsnip_arg(quiet)
  impute_y            <- .as_flag(impute_y, "impute_y")
  quiet               <- .as_flag(quiet, "quiet") %||% TRUE

  rlang::check_installed("PrestoGP", reason = "for the PrestoGP engine")

  model_type <- rlang::arg_match(model_type)
  penalty    <- rlang::arg_match(penalty)
  family     <- rlang::arg_match(family)

  cov_model <- translate_covariance(covariance_function, "PrestoGP", default = "matern")
  if (!identical(cov_model, "matern")) {
    cli::cli_warn(
      "PrestoGP currently supports Matérn covariance in this adapter; falling back to {.val matern}."
    )
  }

  coords <- extract_st_coords(
    data,
    coord_cols = coord_cols,
    time_col   = time_col,
    time_scale = time_scale
  )
  plain_data <- drop_geometry(data)
  parsed     <- parse_formula(formula, plain_data)

  response <- as.numeric(parsed$response)
  if (!is.numeric(response)) {
    cli::cli_abort("Outcome variable must be numeric for the PrestoGP engine.")
  }

  tt    <- stats::terms(formula)
  preds <- attr(tt, "term.labels")

  X <- if (length(preds) > 0L) {
    stats::model.matrix(
      stats::reformulate(preds, intercept = FALSE),
      data = plain_data
    )
  } else {
    matrix(1, nrow = nrow(plain_data), ncol = 1, dimnames = list(NULL, "(Intercept)"))
  }

  if (is.null(scaling) && !is.null(time_col) && ncol(coords) >= 3L) {
    scaling <- c(rep(1L, ncol(coords) - 1L), 2L)
  }

  if (is.null(impute_y)) {
    impute_y <- anyNA(response) || !is.null(lod_upper) || !is.null(lod_lower)
  }
  impute_y <- .as_flag(impute_y, "impute_y")

  na_y <- sum(is.na(response))
  na_x <- sum(is.na(X))
  na_l <- sum(is.na(coords))
  if (na_x > 0L || na_l > 0L || (na_y > 0L && !isTRUE(impute_y))) {
    cli::cli_abort(
      c(
        "PrestoGP input contains missing values that cannot be fitted.",
        "i" = "Missing counts: Y={na_y}, X={na_x}, locs={na_l}, impute_y={impute_y}.",
        "i" = "Set {.arg impute_y = TRUE} for missing outcomes and remove NAs in predictors/coordinates/time."
      )
    )
  }

  model_obj <- if (identical(model_type, "vecchia")) {
    methods::new("VecchiaModel", n_neighbors = as.integer(n_neighbors))
  } else {
    methods::new("FullModel")
  }

  fit_obj <- tryCatch(
    PrestoGP::prestogp_fit(
      model       = model_obj,
      Y           = response,
      X           = X,
      locs        = coords,
      scaling     = scaling,
      common_scale = common_scale,
      impute.y    = impute_y,
      lod.upper   = lod_upper,
      lod.lower   = lod_lower,
      n.impute    = as.integer(n_impute),
      eps.impute  = eps_impute,
      maxit.impute = as.integer(maxit_impute),
      penalty     = penalty,
      alpha       = alpha,
      family      = family,
      quiet       = quiet,
      ...
    ),
    error = function(e) {
      cli::cli_abort(
        c("PrestoGP::prestogp_fit() failed.", "x" = conditionMessage(e))
      )
    }
  )

  structure(
    list(
      prestogp_fit  = fit_obj,
      formula       = formula,
      predictor_terms = preds,
      x_colnames    = colnames(X),
      coord_cols    = coord_cols,
      time_col      = time_col,
      time_scale    = time_scale,
      n_neighbors   = as.integer(n_neighbors),
      model_type    = model_type
    ),
    class = "gopher_PrestoGP_fit"
  )
}

#' Predict from a PrestoGP-fitted Gaussian Process model
#'
#' @param object A `"gopher_PrestoGP_fit"` returned by [PrestoGP_gp_fit()].
#' @param new_data An `sf` object or `data.frame` with coordinates.
#' @param type `"numeric"` (default) or `"pred_int"`.
#' @param level Confidence level for prediction intervals (default `0.95`).
#' @param m_pred Optional number of neighbors for prediction.
#' @param coord_cols Character(2) coordinate column names for non-sf input.
#' @param time_col Optional character scalar specifying a time column for
#'   spatiotemporal modelling.
#' @param time_scale Numeric scalar used to rescale time when `time_col` is
#'   provided. Default `1`.
#' @param ... Forwarded to `PrestoGP::prestogp_predict()`.
#'
#' @return A [tibble::tibble()] with prediction columns.
#' @export
PrestoGP_gp_predict <- function(
    object,
    new_data,
    type       = "numeric",
    level      = 0.95,
    m_pred     = NULL,
    coord_cols = NULL,
    time_col   = NULL,
    time_scale = 1,
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
  m_pred     <- .eval_parsnip_arg(m_pred)
  coord_cols <- .eval_parsnip_arg(coord_cols)
  time_col   <- .eval_parsnip_arg(time_col)
  time_scale <- .eval_parsnip_arg(time_scale)

  rlang::check_installed("PrestoGP", reason = "for the PrestoGP engine")
  type <- rlang::arg_match(type, c("numeric", "pred_int"))

  if (identical(object$model_type, "full")) {
    cli::cli_abort(
      "Prediction is not supported by PrestoGP full models. Use {.val model_type = 'vecchia'}."
    )
  }

  use_coord_cols <- coord_cols %||% object$coord_cols
  use_time_col   <- time_col %||% object$time_col
  use_time_scale <- time_scale %||% object$time_scale

  coords_new <- extract_st_coords(
    new_data,
    coord_cols = use_coord_cols,
    time_col   = use_time_col,
    time_scale = use_time_scale
  )
  plain_new <- drop_geometry(new_data)

  X_new <- if (length(object$predictor_terms) > 0L) {
    stats::model.matrix(
      stats::reformulate(object$predictor_terms, intercept = FALSE),
      data = plain_new
    )
  } else {
    matrix(1, nrow = nrow(plain_new), ncol = 1, dimnames = list(NULL, "(Intercept)"))
  }

  if (!is.null(object$x_colnames)) {
    missing_cols <- setdiff(object$x_colnames, colnames(X_new))
    if (length(missing_cols) > 0L) {
      add <- matrix(
        0,
        nrow = nrow(X_new),
        ncol = length(missing_cols),
        dimnames = list(NULL, missing_cols)
      )
      X_new <- cbind(X_new, add)
    }
    X_new <- X_new[, object$x_colnames, drop = FALSE]
  }

  pred_result <- tryCatch(
    PrestoGP::prestogp_predict(
      model         = object$prestogp_fit,
      X             = X_new,
      locs          = coords_new,
      m             = as.integer(m_pred %||% object$n_neighbors),
      return.values = if (type == "pred_int") "meanvar" else "mean",
      ...
    ),
    error = function(e) {
      cli::cli_abort(
        c("PrestoGP::prestogp_predict() failed.", "x" = conditionMessage(e))
      )
    }
  )

  preds <- as.numeric(pred_result$means)
  if (type == "numeric") {
    return(tibble::tibble(.pred = preds))
  }

  sds <- as.numeric(pred_result$sds)
  if (length(sds) != length(preds) || anyNA(sds)) {
    cli::cli_abort(
      "PrestoGP did not return valid prediction standard deviations for {.val type = 'pred_int'}."
    )
  }

  z <- stats::qnorm((1 + level) / 2)
  tibble::tibble(
    .pred       = preds,
    .pred_lower = preds - z * sds,
    .pred_upper = preds + z * sds
  )
}
