# Model registration with parsnip ----------------------------------------
#
# This file runs when the package is loaded. It registers the
# gaussian_process_spatial model type and all its engines with parsnip's
# internal model database.

.onLoad <- function(libname, pkgname) {
  make_gaussian_process_spatial()
}

# ---- Registration function ----------------------------------------------

#' @keywords internal
make_gaussian_process_spatial <- function() {

  # ---- Model type ---------------------------------------------------------
  parsnip::set_new_model("gaussian_process_spatial")
  parsnip::set_model_mode("gaussian_process_spatial", "regression")

  # ---- Helper: register an engine with shared args ----------------------
  .register_engine <- function(eng, pkg) {

    parsnip::set_model_engine(
      "gaussian_process_spatial",
      mode = "regression",
      eng  = eng
    )

    parsnip::set_dependency(
      "gaussian_process_spatial",
      eng = eng,
      pkg = pkg
    )

    # Shared model arguments -------------------------------------------

    parsnip::set_model_arg(
      model       = "gaussian_process_spatial",
      eng         = eng,
      parsnip     = "covariance_function",
      original    = "covariance_function",
      func        = list(pkg = "gopher", fun = "covariance_function"),
      has_submodel = FALSE
    )

    parsnip::set_model_arg(
      model       = "gaussian_process_spatial",
      eng         = eng,
      parsnip     = "range",
      original    = "range",
      func        = list(pkg = "gopher", fun = "gp_range"),
      has_submodel = FALSE
    )

    parsnip::set_model_arg(
      model       = "gaussian_process_spatial",
      eng         = eng,
      parsnip     = "nugget",
      original    = "nugget",
      func        = list(pkg = "gopher", fun = "gp_nugget"),
      has_submodel = FALSE
    )

    parsnip::set_model_arg(
      model       = "gaussian_process_spatial",
      eng         = eng,
      parsnip     = "sill",
      original    = "sill",
      func        = list(pkg = "gopher", fun = "gp_sill"),
      has_submodel = FALSE
    )

    # Data encoding: keep raw coordinates, no dummy encoding -------------
    parsnip::set_encoding(
      model   = "gaussian_process_spatial",
      eng     = eng,
      mode    = "regression",
      options = list(
        predictor_indicators = "none",
        compute_intercept    = FALSE,
        remove_intercept     = FALSE,
        allow_sparse_x       = FALSE
      )
    )
  }

  # ---- Register each engine -------------------------------------------

  # gstat
  .register_engine("gstat", "gstat")

  parsnip::set_fit(
    model = "gaussian_process_spatial",
    eng   = "gstat",
    mode  = "regression",
    value = list(
      interface = "formula",
      protect   = c("formula", "data"),
      func      = c(pkg = "gopher", fun = "gstat_gp_fit"),
      defaults  = list()
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "gstat",
    mode  = "regression",
    type  = "numeric",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "gstat_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data)
      )
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "gstat",
    mode  = "regression",
    type  = "pred_int",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "gstat_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data),
        type     = "pred_int",
        level    = rlang::expr(level)
      )
    )
  )

  # fields
  .register_engine("fields", "fields")

  parsnip::set_fit(
    model = "gaussian_process_spatial",
    eng   = "fields",
    mode  = "regression",
    value = list(
      interface = "formula",
      protect   = c("formula", "data"),
      func      = c(pkg = "gopher", fun = "fields_gp_fit"),
      defaults  = list()
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "fields",
    mode  = "regression",
    type  = "numeric",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "fields_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data)
      )
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "fields",
    mode  = "regression",
    type  = "pred_int",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "fields_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data),
        type     = "pred_int",
        level    = rlang::expr(level)
      )
    )
  )

  # GPvecchia
  .register_engine("GPvecchia", "GPvecchia")

  parsnip::set_fit(
    model = "gaussian_process_spatial",
    eng   = "GPvecchia",
    mode  = "regression",
    value = list(
      interface = "formula",
      protect   = c("formula", "data"),
      func      = c(pkg = "gopher", fun = "GPvecchia_gp_fit"),
      defaults  = list()
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "GPvecchia",
    mode  = "regression",
    type  = "numeric",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "GPvecchia_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data)
      )
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "GPvecchia",
    mode  = "regression",
    type  = "pred_int",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "GPvecchia_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data),
        type     = "pred_int",
        level    = rlang::expr(level)
      )
    )
  )

  # spNNGP
  .register_engine("spNNGP", "spNNGP")

  parsnip::set_fit(
    model = "gaussian_process_spatial",
    eng   = "spNNGP",
    mode  = "regression",
    value = list(
      interface = "formula",
      protect   = c("formula", "data"),
      func      = c(pkg = "gopher", fun = "spNNGP_gp_fit"),
      defaults  = list()
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "spNNGP",
    mode  = "regression",
    type  = "numeric",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "spNNGP_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data)
      )
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "spNNGP",
    mode  = "regression",
    type  = "pred_int",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "spNNGP_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data),
        type     = "pred_int",
        level    = rlang::expr(level)
      )
    )
  )

  # PrestoGP
  .register_engine("PrestoGP", "PrestoGP")

  parsnip::set_fit(
    model = "gaussian_process_spatial",
    eng   = "PrestoGP",
    mode  = "regression",
    value = list(
      interface = "formula",
      protect   = c("formula", "data"),
      func      = c(pkg = "gopher", fun = "PrestoGP_gp_fit"),
      defaults  = list()
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "PrestoGP",
    mode  = "regression",
    type  = "numeric",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "PrestoGP_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data)
      )
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "PrestoGP",
    mode  = "regression",
    type  = "pred_int",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "PrestoGP_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data),
        type     = "pred_int",
        level    = rlang::expr(level)
      )
    )
  )

  # sdmTMB
  .register_engine("sdmTMB", "sdmTMB")

  parsnip::set_fit(
    model = "gaussian_process_spatial",
    eng   = "sdmTMB",
    mode  = "regression",
    value = list(
      interface = "formula",
      protect   = c("formula", "data"),
      func      = c(pkg = "gopher", fun = "sdmTMB_gp_fit"),
      defaults  = list()
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "sdmTMB",
    mode  = "regression",
    type  = "numeric",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "sdmTMB_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data)
      )
    )
  )

  parsnip::set_pred(
    model = "gaussian_process_spatial",
    eng   = "sdmTMB",
    mode  = "regression",
    type  = "pred_int",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "gopher", fun = "sdmTMB_gp_predict"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data),
        type     = "pred_int",
        level    = rlang::expr(level)
      )
    )
  )

  invisible(NULL)
}
