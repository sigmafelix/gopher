# Tests for the sdmTMB engine -------------------------------------------

skip_if_not_installed("sdmTMB")

# ---- sdmTMB_gp_fit ---------------------------------------------------

test_that("sdmTMB_gp_fit returns gopher_sdmTMB_fit object (sf input)", {
  train <- make_spatial_sf(40L)
  fit   <- sdmTMB_gp_fit(z ~ 1, data = train)
  expect_s3_class(fit, "gopher_sdmTMB_fit")
  expect_false(is.null(fit$sdmtmb_fit))
  expect_false(is.null(fit$mesh))
})

test_that("sdmTMB_gp_fit returns gopher_sdmTMB_fit object (data.frame input)", {
  train <- make_spatial_df(40L)
  fit   <- sdmTMB_gp_fit(z ~ 1, data = train)
  expect_s3_class(fit, "gopher_sdmTMB_fit")
})

test_that("sdmTMB_gp_fit stores formula and xy_cols", {
  train <- make_spatial_sf(40L)
  fit   <- sdmTMB_gp_fit(z ~ 1, data = train)
  expect_equal(fit$formula, z ~ 1)
  expect_length(fit$xy_cols, 2L)
})

test_that("sdmTMB_gp_fit accepts mesh_cutoff argument", {
  train <- make_spatial_sf(40L)
  fit   <- sdmTMB_gp_fit(z ~ 1, data = train, mesh_cutoff = 2)
  expect_s3_class(fit, "gopher_sdmTMB_fit")
})

test_that("sdmTMB_gp_fit accepts n_knots argument", {
  train <- make_spatial_sf(40L)
  fit   <- sdmTMB_gp_fit(z ~ 1, data = train, n_knots = 20L)
  expect_s3_class(fit, "gopher_sdmTMB_fit")
})

test_that("sdmTMB_gp_fit supports universal kriging (covariates)", {
  train <- make_spatial_sf(40L)
  fit   <- sdmTMB_gp_fit(z ~ cov1, data = train)
  expect_s3_class(fit, "gopher_sdmTMB_fit")
})

test_that("sdmTMB_gp_fit warns for non-Matern covariance", {
  train <- make_spatial_sf(40L)
  expect_warning(
    sdmTMB_gp_fit(z ~ 1, data = train, covariance_function = "exponential"),
    regexp = "Matern"
  )
})

test_that("sdmTMB_gp_fit accepts matern covariance without warning", {
  train <- make_spatial_sf(40L)
  expect_no_warning(
    sdmTMB_gp_fit(z ~ 1, data = train, covariance_function = "matern")
  )
})

# ---- sdmTMB_gp_predict -----------------------------------------------

test_that("sdmTMB_gp_predict returns tibble with .pred column", {
  train  <- make_spatial_sf(40L)
  newdat <- make_pred_sf()
  fit    <- sdmTMB_gp_fit(z ~ 1, data = train)
  preds  <- sdmTMB_gp_predict(fit, new_data = newdat)
  expect_s3_class(preds, "tbl_df")
  expect_named(preds, ".pred")
  expect_equal(nrow(preds), nrow(newdat))
  expect_true(all(is.finite(preds$.pred)))
})

test_that("sdmTMB_gp_predict with type='pred_int' returns interval columns", {
  train  <- make_spatial_sf(40L)
  newdat <- make_pred_sf()
  fit    <- sdmTMB_gp_fit(z ~ 1, data = train)
  preds  <- sdmTMB_gp_predict(fit, new_data = newdat, type = "pred_int")
  expect_named(preds, c(".pred", ".pred_lower", ".pred_upper"))
  # small tolerance for floating-point rounding at interval boundaries
  expect_true(all(preds$.pred_lower <= preds$.pred + 1e-9))
  expect_true(all(preds$.pred_upper >= preds$.pred - 1e-9))
})

test_that("sdmTMB_gp_predict works with data.frame new_data", {
  train  <- make_spatial_sf(40L)
  newdat <- make_pred_df()
  fit    <- sdmTMB_gp_fit(z ~ 1, data = train)
  preds  <- sdmTMB_gp_predict(fit, new_data = newdat)
  expect_equal(nrow(preds), nrow(newdat))
})

test_that("sdmTMB_gp_predict works with universal kriging", {
  train  <- make_spatial_sf(40L)
  newdat <- make_pred_sf()
  fit    <- sdmTMB_gp_fit(z ~ cov1, data = train)
  preds  <- sdmTMB_gp_predict(fit, new_data = newdat)
  expect_equal(nrow(preds), nrow(newdat))
})

# ---- parsnip integration --------------------------------------------

test_that("parsnip fit() + predict() pipeline works with sdmTMB", {
  train  <- make_spatial_sf(40L)
  newdat <- make_pred_sf()

  spec <- gaussian_process_spatial(covariance_function = "matern") |>
    parsnip::set_engine("sdmTMB")

  fitted <- parsnip::fit(spec, z ~ 1, data = train)
  expect_s3_class(fitted, "model_fit")

  preds <- predict(fitted, new_data = newdat)
  expect_named(preds, ".pred")
  expect_equal(nrow(preds), nrow(newdat))
})

test_that("parsnip predict(type='pred_int') pipeline works with sdmTMB", {
  train  <- make_spatial_sf(40L)
  newdat <- make_pred_sf()

  spec   <- gaussian_process_spatial() |> parsnip::set_engine("sdmTMB")
  fitted <- parsnip::fit(spec, z ~ 1, data = train)
  preds  <- predict(fitted, new_data = newdat, type = "pred_int")
  expect_true(".pred_lower" %in% names(preds))
  expect_true(".pred_upper" %in% names(preds))
})

test_that("sdmTMB engine works with non-sf data in parsnip pipeline", {
  train  <- make_spatial_df(40L)
  newdat <- make_pred_df()

  spec   <- gaussian_process_spatial() |> parsnip::set_engine("sdmTMB")
  fitted <- parsnip::fit(spec, z ~ 1, data = train)
  preds  <- predict(fitted, new_data = newdat)
  expect_equal(nrow(preds), nrow(newdat))
})
