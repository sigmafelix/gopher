# spNNGP engine for gaussian_process_spatial --------------------------
#
# Uses the spNNGP package for Nearest Neighbor Gaussian Process (NNGP) models.
#
# spNNGP workflow:
#   m <- spNNGP::spNNGP(formula, data = data, coords = ...,
#                       starting = ..., tuning = ..., priors = ...,
#                       cov.model = ..., n.neighbors = ...,
#                       method = "response")
#   predict(m, X.0 = ..., coords.0 = ...)


#' Fit a Gaussian Process model using the spNNGP engine
#'
#' @description
#' This function is called internally by parsnip. End users should call
#' [parsnip::fit.model_spec()] on a `gaussian_process_spatial` specification
#' with `set_engine("spNNGP")`.
#'
#' @param formula  A two-sided formula. Covariates are included as fixed
#'   effects (Universal Kriging).
#' @param data     An `sf` object or `data.frame` with coordinate columns.
#' @param covariance_function Canonical covariance name. Defaults to
#'   `"exponential"`.
#' @param range    Starting value for the range parameter (`phi`).
#'   `NULL` = uses a data-driven initial value.
#' @param nugget   Starting value for the nugget (`tau.sq`).
#'   `NULL` = 10 \% of response variance.
#' @param sill     Starting value for the partial sill (`sigma.sq`).
#'   `NULL` = 90 \% of response variance.
#' @param n_neighbors Number of nearest neighbours. Default `15`.
#' @param n_samples   Number of MCMC samples. Default `1000`.
#' @param n_burnin    MCMC burn-in. Default `500`.
#' @param method    spNNGP method: `"response"` (default) or `"latent"`.
#' @param coord_cols Character(2) coordinate column names (non-sf path).
#' @param time_col Optional character scalar specifying a time column for
#'   spatiotemporal modelling.
#' @param time_scale Numeric scalar used to rescale time when `time_col` is
#'   provided. Default `1`.
#' @param ... Additional arguments forwarded to `spNNGP::spNNGP()`.
#'
#' @return A list of class `"gopher_spNNGP_fit"`.
#' @examples
#' if (requireNamespace("spacetime", quietly = TRUE) &&
#'     requireNamespace("spNNGP", quietly = TRUE)) {
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
#'   fit <- spNNGP_gp_fit(
#'     pm10 ~ coords.x1 + coords.x2,
#'     data = air_sf,
#'     n_neighbors = 8,
#'     n_samples = 80,
#'     n_burnin = 40
#'   )
#'   fit
#' }
#'
#' @export
spNNGP_gp_fit <- function(
    formula,
    data,
    covariance_function = "exponential",
    range               = NULL,
    nugget              = NULL,
    sill                = NULL,
    n_neighbors         = 15L,
    n_samples           = 1000L,
    n_burnin            = 500L,
    method              = "response",
    coord_cols          = NULL,
    time_col            = NULL,
    time_scale          = 1,
    ...) {

  rlang::check_installed("spNNGP", reason = "for the spNNGP engine")

  coords     <- extract_st_coords(
    data,
    coord_cols = coord_cols,
    time_col   = time_col,
    time_scale = time_scale
  )
  plain_data <- drop_geometry(data)
  parsed     <- parse_formula(formula, plain_data)

  response <- parsed$response
  y_var    <- stats::var(response, na.rm = TRUE)

  cov_model <- translate_covariance(
    covariance_function, "spNNGP", default = "exponential"
  )

  # ---- Starting / prior values ------------------------------------------
  init_sigma_sq <- sill   %||% (y_var * 0.9)
  init_tau_sq   <- nugget %||% (y_var * 0.1)
  init_phi      <- range  %||% {
    dists <- as.vector(stats::dist(coords))
    3 / stats::quantile(dists, 0.75, na.rm = TRUE)  # effective range ~ 75th pct
  }

  starting <- list(
    "phi"      = init_phi,
    "sigma.sq" = init_sigma_sq,
    "tau.sq"   = init_tau_sq
  )

  tuning <- list(
    "phi"      = 0.5,
    "sigma.sq" = 0.5,
    "tau.sq"   = 0.5
  )

  priors <- list(
    "phi.Unif"      = c(0.1, 10 * init_phi),
    "sigma.sq.IG"   = c(2, init_sigma_sq),
    "tau.sq.IG"     = c(2, init_tau_sq)
  )

  fit_obj <- tryCatch(
    spNNGP::spNNGP(
      formula    = formula,
      data       = plain_data,
      coords     = coords,
      starting   = starting,
      tuning     = tuning,
      priors     = priors,
      cov.model  = cov_model,
      n.neighbors = as.integer(n_neighbors),
      n.samples  = as.integer(n_samples),
      method     = method,
      verbose    = FALSE,
      ...
    ),
    error = function(e) {
      cli::cli_abort(
        c("spNNGP::spNNGP() failed.", "x" = conditionMessage(e))
      )
    }
  )

  structure(
    list(
      spNNGP_fit    = fit_obj,
      coords        = coords,
      plain_data    = plain_data,
      formula       = formula,
      n_burnin      = as.integer(n_burnin),
      n_samples     = as.integer(n_samples),
      cov_model     = cov_model
    ),
    class = "gopher_spNNGP_fit"
  )
}

#' Predict from an spNNGP-fitted Gaussian Process model
#'
#' @param object   A `"gopher_spNNGP_fit"` returned by [spNNGP_gp_fit()].
#' @param new_data An `sf` object or `data.frame` with coordinates.
#' @param type     `"numeric"` (default) or `"pred_int"`.
#' @param level    Confidence level for prediction intervals (default `0.95`).
#' @param coord_cols Character(2) coord column names (non-sf path).
#' @param time_col Optional character scalar specifying a time column for
#'   spatiotemporal modelling.
#' @param time_scale Numeric scalar used to rescale time when `time_col` is
#'   provided. Default `1`.
#' @param ... Forwarded to `predict.spNNGP()`.
#'
#' @return A [tibble::tibble()] with prediction columns.
#' @examples
#' if (requireNamespace("spacetime", quietly = TRUE) &&
#'     requireNamespace("spNNGP", quietly = TRUE)) {
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
#'   fit <- spNNGP_gp_fit(
#'     pm10 ~ coords.x1 + coords.x2,
#'     data = train_sf,
#'     n_neighbors = 8,
#'     n_samples = 80,
#'     n_burnin = 40
#'   )
#'   spNNGP_gp_predict(fit, new_data = test_sf, type = "pred_int")
#' }
#'
#' @export
spNNGP_gp_predict <- function(
    object,
    new_data,
    type       = "numeric",
    level      = 0.95,
    coord_cols = NULL,
    time_col   = NULL,
    time_scale = 1,
    ...) {

  rlang::check_installed("spNNGP", reason = "for the spNNGP engine")
  type <- rlang::arg_match(type, c("numeric", "pred_int"))

  coords_new <- extract_st_coords(
    new_data,
    coord_cols = coord_cols,
    time_col   = time_col,
    time_scale = time_scale
  )
  plain_new  <- drop_geometry(new_data)

  # Design matrix for new locations
  tt    <- stats::terms(object$formula)
  preds <- attr(tt, "term.labels")

  X_new <- if (length(preds) > 0L) {
    available <- intersect(preds, names(plain_new))
    if (length(available) == 0L) {
      cli::cli_abort(
        "The model was fitted with covariates but {.arg new_data} does not
         contain the required predictor columns."
      )
    }
    stats::model.matrix(stats::reformulate(preds), data = plain_new)
  } else {
    matrix(1, nrow = nrow(coords_new), ncol = 1)
  }

  # Burn-in indices
  n_s      <- object$n_samples
  n_b      <- object$n_burnin
  keep_idx <- seq(n_b + 1L, n_s)

  pred_result <- tryCatch(
    predict(
      object$spNNGP_fit,
      X.0      = X_new,
      coords.0 = coords_new,
      sub.sample = list(start = n_b + 1L, end = n_s),
      verbose  = FALSE,
      ...
    ),
    error = function(e) {
      cli::cli_abort(
        c("predict.spNNGP failed.", "x" = conditionMessage(e))
      )
    }
  )

  # pred_result$p.y.0 is n_pred x n_post_samples matrix
  p_y <- pred_result$p.y.0
  preds_vec <- rowMeans(p_y)

  if (type == "numeric") {
    return(tibble::tibble(.pred = preds_vec))
  }

  alpha   <- 1 - level
  q_lo    <- alpha / 2
  q_hi    <- 1 - alpha / 2
  lo      <- apply(p_y, 1, stats::quantile, probs = q_lo, na.rm = TRUE)
  hi      <- apply(p_y, 1, stats::quantile, probs = q_hi, na.rm = TRUE)
  tibble::tibble(
    .pred       = preds_vec,
    .pred_lower = lo,
    .pred_upper = hi
  )
}
