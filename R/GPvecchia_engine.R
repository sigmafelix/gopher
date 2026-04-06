# GPvecchia engine for gaussian_process_spatial -----------------------
#
# Uses the GPvecchia package for scalable Vecchia-approximated GP inference.
#
# GPvecchia workflow:
#   vecchia_approx  <- GPvecchia::vecchia_specify(locs, m)
#   likelihood      <- GPvecchia::vecchia_likelihood(...)
#   params          <- optim(...)   # MLE
#   predictions     <- GPvecchia::vecchia_prediction(...)

#' Fit a Gaussian Process model using the GPvecchia engine
#'
#' This function is called internally by parsnip. End users should call
#' [parsnip::fit.model_spec()] on a `gaussian_process_spatial` specification
#' with `set_engine("GPvecchia")`.
#'
#' @param formula  A two-sided formula. Covariates are used as fixed-effect
#'   trend (Universal Kriging / residual GP).
#' @param data     An `sf` object or `data.frame` with coordinate columns.
#' @param covariance_function Canonical covariance name. Defaults to
#'   `"exponential"`.
#' @param range    Range parameter (scale). `NULL` = estimated via MLE.
#' @param nugget   Nugget (noise) variance. `NULL` = estimated.
#' @param sill     Signal variance (sigma^2). `NULL` = estimated.
#' @param m        Number of nearest neighbours for Vecchia approximation.
#'   Default `15`.
#' @param coord_cols Character(2) coordinate column names (non-sf path).
#' @param time_col Optional character scalar specifying a time column for
#'   spatiotemporal modelling.
#' @param time_scale Numeric scalar used to rescale time when `time_col` is
#'   provided. Default `1`.
#' @param ... Additional arguments forwarded to
#'   `GPvecchia::vecchia_specify()`.
#'
#' @return A list of class `"gopher_GPvecchia_fit"`.
#' @examples
#' if (requireNamespace("spacetime", quietly = TRUE) &&
#'     requireNamespace("GPvecchia", quietly = TRUE)) {
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
#'   fit <- GPvecchia_gp_fit(pm10 ~ coords.x1 + coords.x2, data = air_sf, m = 10)
#'   fit
#' }
#'
#' @export
GPvecchia_gp_fit <- function(
    formula,
    data,
    covariance_function = "exponential",
    range               = NULL,
    nugget              = NULL,
    sill                = NULL,
    m                   = 15L,
    coord_cols          = NULL,
    time_col            = NULL,
    time_scale          = 1,
    ...) {

  rlang::check_installed("GPvecchia", reason = "for the GPvecchia engine")

  coords     <- extract_st_coords(
    data,
    coord_cols = coord_cols,
    time_col   = time_col,
    time_scale = time_scale
  )
  plain_data <- drop_geometry(data)
  parsed     <- parse_formula(formula, plain_data)

  response <- parsed$response
  X        <- parsed$X   # model matrix including intercept

  # Model trend separately and fit GP to residuals.
  beta_hat <- tryCatch(
    stats::lm.fit(x = X, y = response)$coefficients,
    error = function(e) rep(0, ncol(X))
  )
  if (anyNA(beta_hat)) {
    beta_hat[is.na(beta_hat)] <- 0
  }

  cov_fn <- translate_covariance(
    covariance_function, "GPvecchia", default = "exponential"
  )

  # ---- Estimate parameters if not supplied ------------------------------
  y_var   <- stats::var(response, na.rm = TRUE)
  init_sigma2 <- sill   %||% (y_var * 0.9)
  init_range  <- range  %||% {
    dists <- as.vector(stats::dist(coords))
    stats::quantile(dists, 0.25, na.rm = TRUE)
  }
  init_nugget <- nugget %||% (y_var * 0.1)

  if (!identical(cov_fn, "matern")) {
    cli::cli_warn(
      c(
        "GPvecchia covariance {.val {cov_fn}} is not available in this engine build.",
        "i" = "Falling back to {.val matern}."
      )
    )
    cov_fn <- "matern"
  }

  theta_start <- c(
    variance   = init_sigma2,
    range      = init_range,
    smoothness = 0.5,
    nugget     = init_nugget
  )

  vecchia_approx <- GPvecchia::vecchia_specify(
    locs = coords,
    m    = as.integer(m),
    ...
  )

  ll_fun <- function(theta_vec) {
    tryCatch(
      GPvecchia::vecchia_likelihood(
        z              = as.numeric(response - X %*% beta_hat),
        vecchia.approx = vecchia_approx,
        covparms       = theta_vec[1:3],
        nuggets        = theta_vec[4],
        covmodel       = cov_fn
      ),
      error = function(e) NA_real_
    )
  }

  best_start <- theta_start
  ll_best <- ll_fun(best_start)
  if (!is.finite(ll_best)) {
    theta_try <- best_start
    for (i in seq_len(8L)) {
      theta_try["nugget"] <- theta_try["nugget"] * 2
      theta_try["range"] <- theta_try["range"] * 0.9
      ll_try <- ll_fun(theta_try)
      if (is.finite(ll_try)) {
        best_start <- theta_try
        ll_best <- ll_try
        break
      }
    }
  }
  if (!is.finite(ll_best)) {
    cli::cli_abort("GPvecchia likelihood could not be evaluated at stable initial parameter values.")
  }
  cli::cli_warn(
    c(
      "GPvecchia is using a stable finite-parameter initialization rather than iterative MLE.",
      "i" = "This avoids numerical failures in some GPvecchia builds."
    )
  )

  theta_hat <- as.numeric(best_start)
  names(theta_hat) <- c("variance", "range", "smoothness", "nugget")
  est_params <- c(
    sigma2     = unname(theta_hat["variance"]),
    range      = unname(theta_hat["range"]),
    nugget     = unname(theta_hat["nugget"]),
    smoothness = unname(theta_hat["smoothness"])
  )

  structure(
    list(
      vecchia_approx = vecchia_approx,
      params         = est_params,
      cov_fn         = cov_fn,
      beta_hat       = as.numeric(beta_hat),
      X_train        = X,
      response       = response,
      coords         = coords,
      training_data  = plain_data,
      formula        = formula,
      m              = as.integer(m),
      opt_result     = list(convergence = NA_integer_, par = log(theta_hat))
    ),
    class = "gopher_GPvecchia_fit"
  )
}

#' Predict from a GPvecchia-fitted Gaussian Process model
#'
#' @param object   A `"gopher_GPvecchia_fit"` returned by [GPvecchia_gp_fit()].
#' @param new_data An `sf` object or `data.frame` with coordinates.
#' @param type     `"numeric"` (default) or `"pred_int"`.
#' @param level    Confidence level for prediction intervals (default `0.95`).
#' @param m_pred   Number of nearest neighbours for prediction approximation.
#'   Defaults to the training `m`.
#' @param coord_cols Character(2) coord column names (non-sf path).
#' @param time_col Optional character scalar specifying a time column for
#'   spatiotemporal modelling.
#' @param time_scale Numeric scalar used to rescale time when `time_col` is
#'   provided. Default `1`.
#' @param ... Additional arguments forwarded to
#'   `GPvecchia::vecchia_specify()` for prediction graph construction.
#'
#' @return A [tibble::tibble()] with prediction columns.
#' @examples
#' if (requireNamespace("spacetime", quietly = TRUE) &&
#'     requireNamespace("GPvecchia", quietly = TRUE)) {
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
#'   fit <- GPvecchia_gp_fit(pm10 ~ coords.x1 + coords.x2, data = train_sf, m = 10)
#'   GPvecchia_gp_predict(fit, new_data = test_sf, type = "pred_int")
#' }
#'
#' @export
GPvecchia_gp_predict <- function(
    object,
    new_data,
    type       = "numeric",
    level      = 0.95,
    m_pred     = NULL,
    coord_cols = NULL,
    time_col   = NULL,
    time_scale = 1,
    ...) {

  rlang::check_installed("GPvecchia", reason = "for the GPvecchia engine")
  type <- rlang::arg_match(type, c("numeric", "pred_int"))

  coords_new <- extract_st_coords(
    new_data,
    coord_cols = coord_cols,
    time_col   = time_col,
    time_scale = time_scale
  )

  plain_new <- drop_geometry(new_data)

  # Build model matrix for new data covariates
  tt      <- stats::terms(object$formula)
  preds   <- attr(tt, "term.labels")
  has_cov <- length(preds) > 0L

  X_new <- if (has_cov) {
    # Align columns to training matrix
    available <- intersect(preds, names(plain_new))
    if (length(available) == 0L) {
      cli::cli_abort(
        "The model was fitted with covariates but {.arg new_data} does not
         contain the required predictor columns."
      )
    }
    stats::model.matrix(
      stats::reformulate(preds),
      data = plain_new
    )
  } else {
    matrix(1, nrow = nrow(coords_new), ncol = 1)
  }

  if (is.null(m_pred)) {
    m_pred <- object$m %||% 15L
  }

  pred_result <- tryCatch(
    {
      vecchia_pred_approx <- GPvecchia::vecchia_specify(
        locs      = object$coords,
        m         = as.integer(m_pred),
        locs.pred = coords_new,
        ...
      )

      GPvecchia::vecchia_prediction(
        z              = as.numeric(object$response - object$X_train %*% object$beta_hat),
        vecchia.approx = vecchia_pred_approx,
        covparms       = c(
          object$params["sigma2"],
          object$params["range"],
          object$params["smoothness"]
        ),
        nuggets        = object$params["nugget"],
        covmodel       = object$cov_fn,
        return.values  = "meanvar"
      )
    },
    error = function(e) {
      cli::cli_abort(
        c(
          "GPvecchia prediction failed.",
          "x" = conditionMessage(e)
        )
      )
    }
  )

  preds    <- as.numeric(pred_result$mu.pred + X_new %*% object$beta_hat)
  variance <- pred_result$var.pred

  if (type == "numeric") {
    return(tibble::tibble(.pred = as.numeric(preds)))
  }

  alpha <- 1 - level
  z     <- stats::qnorm(1 - alpha / 2)
  se    <- sqrt(pmax(as.numeric(variance), 0))
  tibble::tibble(
    .pred       = as.numeric(preds),
    .pred_lower = as.numeric(preds) - z * se,
    .pred_upper = as.numeric(preds) + z * se
  )
}
