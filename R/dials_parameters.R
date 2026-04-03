# dials parameter definitions for tuning gaussian_process_spatial -----
#
# These functions return `dials::param` objects so that tidymodels' tune
# infrastructure can explore the GP hyperparameter space.

#' Covariance function type (qualitative tuning parameter)
#'
#' @param values Character vector of allowable covariance function names.
#'   Defaults to the four most common choices.
#'
#' @return A `dials` qualitative parameter object.
#' @export
#' @examples
#' covariance_function()
covariance_function <- function(
    values = c("exponential", "spherical", "gaussian", "matern")) {
  rlang::check_installed("dials", reason = "for GP hyperparameter tuning")
  dials::new_qual_param(
    type   = "character",
    values = values,
    label  = c(covariance_function = "Covariance Function"),
    finalize = NULL
  )
}

#' GP range parameter (quantitative tuning parameter)
#'
#' @param range A numeric vector of length 2 giving the lower and upper bounds
#'   of the range parameter search space. Values are on the original scale.
#' @param trans A `scales` transformation object, or `NULL` (default).
#'
#' @return A `dials` quantitative parameter object.
#' @export
#' @examples
#' gp_range()
gp_range <- function(range = c(1e-3, 1e3), trans = NULL) {
  rlang::check_installed("dials", reason = "for GP hyperparameter tuning")
  dials::new_quant_param(
    type      = "double",
    range     = range,
    inclusive = c(TRUE, TRUE),
    trans     = trans,
    label     = c(gp_range = "GP Range"),
    finalize  = NULL
  )
}

#' GP nugget variance (quantitative tuning parameter)
#'
#' @param range A numeric vector of length 2 giving lower and upper bounds.
#' @param trans A `scales` transformation, or `NULL`.
#'
#' @return A `dials` quantitative parameter object.
#' @export
#' @examples
#' gp_nugget()
gp_nugget <- function(range = c(0, 5), trans = NULL) {
  rlang::check_installed("dials", reason = "for GP hyperparameter tuning")
  dials::new_quant_param(
    type      = "double",
    range     = range,
    inclusive = c(TRUE, FALSE),
    trans     = trans,
    label     = c(gp_nugget = "GP Nugget Variance"),
    finalize  = NULL
  )
}

#' GP partial sill variance (quantitative tuning parameter)
#'
#' @param range A numeric vector of length 2 giving lower and upper bounds.
#' @param trans A `scales` transformation, or `NULL` (default).
#'
#' @return A `dials` quantitative parameter object.
#' @export
#' @examples
#' gp_sill()
gp_sill <- function(range = c(1e-6, 1e3), trans = NULL) {
  rlang::check_installed("dials", reason = "for GP hyperparameter tuning")
  dials::new_quant_param(
    type      = "double",
    range     = range,
    inclusive = c(FALSE, TRUE),
    trans     = trans,
    label     = c(gp_sill = "GP Partial Sill"),
    finalize  = NULL
  )
}
