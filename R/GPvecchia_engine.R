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
#' @param ... Additional arguments forwarded to
#'   `GPvecchia::vecchia_specify()`.
#'
#' @return A list of class `"gopher_GPvecchia_fit"`.
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
    ...) {

  rlang::check_installed("GPvecchia", reason = "for the GPvecchia engine")

  coords     <- extract_coords(data, coord_cols)
  plain_data <- drop_geometry(data)
  parsed     <- parse_formula(formula, plain_data)

  response <- parsed$response
  X        <- parsed$X   # model matrix including intercept

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

  # ---- Vecchia approximation specification ------------------------------
  vecchia_approx <- GPvecchia::vecchia_specify(
    locs = coords,
    m    = as.integer(m),
    ...
  )

  # ---- MLE using optim --------------------------------------------------
  # Log-likelihood wrapper
  neg_log_lik <- function(params) {
    sigma2_v <- exp(params[1])
    range_v  <- exp(params[2])
    nugget_v <- exp(params[3])
    tryCatch(
      -GPvecchia::vecchia_likelihood(
        y          = response,
        vecchia.approx = vecchia_approx,
        covparms   = c(sigma2_v, range_v, nugget_v),
        nugget     = nugget_v,
        covfun.name = cov_fn,
        X          = X
      ),
      error = function(e) Inf
    )
  }

  init_params <- log(c(init_sigma2, init_range, init_nugget))

  opt <- tryCatch(
    stats::optim(
      par    = init_params,
      fn     = neg_log_lik,
      method = "L-BFGS-B"
    ),
    error = function(e) {
      cli::cli_warn(
        c(
          "GPvecchia MLE optimisation failed: {conditionMessage(e)}",
          "i" = "Using initial parameter values."
        )
      )
      list(par = init_params, convergence = 1L)
    }
  )

  est_params <- exp(opt$par)
  names(est_params) <- c("sigma2", "range", "nugget")

  structure(
    list(
      vecchia_approx = vecchia_approx,
      params         = est_params,
      cov_fn         = cov_fn,
      X_train        = X,
      response       = response,
      coords         = coords,
      training_data  = plain_data,
      formula        = formula,
      opt_result     = opt
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
#' @param ... Forwarded to `GPvecchia::vecchia_prediction()`.
#'
#' @return A [tibble::tibble()] with prediction columns.
#' @export
GPvecchia_gp_predict <- function(
    object,
    new_data,
    type       = "numeric",
    level      = 0.95,
    m_pred     = NULL,
    coord_cols = NULL,
    ...) {

  rlang::check_installed("GPvecchia", reason = "for the GPvecchia engine")
  type <- rlang::arg_match(type, c("numeric", "pred_int"))

  coords_new <- extract_coords(new_data, coord_cols)

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

  params <- object$params

  pred_result <- GPvecchia::vecchia_prediction(
    y              = object$response,
    vecchia.approx = object$vecchia_approx,
    covparms       = c(params["sigma2"], params["range"], params["nugget"]),
    nugget         = params["nugget"],
    covfun.name    = object$cov_fn,
    locs.pred      = coords_new,
    X              = object$X_train,
    X.pred         = X_new,
    ...
  )

  preds    <- pred_result$mu.pred
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
