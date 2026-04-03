# Tests for the fields engine ------------------------------------------

skip_if_not_installed("fields")

# ---- fields_gp_fit --------------------------------------------------

test_that("fields_gp_fit returns gopher_fields_fit object (sf input)", {
  train <- make_spatial_sf(30L)
  fit   <- fields_gp_fit(z ~ 1, data = train)
  expect_s3_class(fit, "gopher_fields_fit")
})

test_that("fields_gp_fit returns gopher_fields_fit object (data.frame input)", {
  train <- make_spatial_df(30L)
  fit   <- fields_gp_fit(z ~ 1, data = train)
  expect_s3_class(fit, "gopher_fields_fit")
})

test_that("fields_gp_fit stores Krig object", {
  train <- make_spatial_sf(30L)
  fit   <- fields_gp_fit(z ~ 1, data = train)
  expect_true(!is.null(fit$krig_obj))
})

test_that("fields_gp_fit accepts exponential covariance", {
  train <- make_spatial_sf(30L)
  fit   <- fields_gp_fit(z ~ 1, data = train, covariance_function = "exponential")
  expect_s3_class(fit, "gopher_fields_fit")
})

test_that("fields_gp_fit accepts matern covariance", {
  train <- make_spatial_sf(30L)
  fit   <- fields_gp_fit(z ~ 1, data = train, covariance_function = "matern")
  expect_s3_class(fit, "gopher_fields_fit")
})

test_that("fields_gp_fit supports universal kriging (covariates)", {
  train <- make_spatial_sf(30L)
  fit   <- fields_gp_fit(z ~ cov1, data = train)
  expect_s3_class(fit, "gopher_fields_fit")
  expect_true(fit$has_cov)
})

test_that("fields_gp_fit stores has_cov = FALSE for ordinary kriging", {
  train <- make_spatial_sf(30L)
  fit   <- fields_gp_fit(z ~ 1, data = train)
  expect_false(fit$has_cov)
})

# ---- fields_gp_predict ----------------------------------------------

test_that("fields_gp_predict returns tibble with .pred column", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_sf()
  fit    <- fields_gp_fit(z ~ 1, data = train)
  preds  <- fields_gp_predict(fit, new_data = newdat)
  expect_s3_class(preds, "tbl_df")
  expect_named(preds, ".pred")
  expect_equal(nrow(preds), nrow(newdat))
  expect_true(all(is.finite(preds$.pred)))
})

test_that("fields_gp_predict with type='pred_int' returns interval columns", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_sf()
  fit    <- fields_gp_fit(z ~ 1, data = train)
  preds  <- fields_gp_predict(fit, new_data = newdat, type = "pred_int")
  expect_named(preds, c(".pred", ".pred_lower", ".pred_upper"))
})

test_that("fields_gp_predict works with data.frame new_data", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_df()
  fit    <- fields_gp_fit(z ~ 1, data = train)
  preds  <- fields_gp_predict(fit, new_data = newdat)
  expect_equal(nrow(preds), nrow(newdat))
})

test_that("fields_gp_predict works for universal kriging", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_sf()
  fit    <- fields_gp_fit(z ~ cov1, data = train)
  preds  <- fields_gp_predict(fit, new_data = newdat)
  expect_equal(nrow(preds), nrow(newdat))
})

# ---- parsnip integration --------------------------------------------

test_that("parsnip fit() + predict() pipeline works with fields", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_sf()

  spec <- gaussian_process_spatial(covariance_function = "exponential") |>
    parsnip::set_engine("fields")

  fitted <- parsnip::fit(spec, z ~ 1, data = train)
  expect_s3_class(fitted, "model_fit")

  preds <- predict(fitted, new_data = newdat)
  expect_named(preds, ".pred")
  expect_equal(nrow(preds), nrow(newdat))
})

test_that("parsnip predict(type='pred_int') pipeline works with fields", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_sf()

  spec   <- gaussian_process_spatial() |> parsnip::set_engine("fields")
  fitted <- parsnip::fit(spec, z ~ 1, data = train)
  preds  <- predict(fitted, new_data = newdat, type = "pred_int")
  expect_true(".pred_lower" %in% names(preds))
  expect_true(".pred_upper" %in% names(preds))
})

test_that("fields engine works with data.frame in parsnip pipeline", {
  train  <- make_spatial_df(30L)
  newdat <- make_pred_df()

  spec   <- gaussian_process_spatial() |> parsnip::set_engine("fields")
  fitted <- parsnip::fit(spec, z ~ 1, data = train)
  preds  <- predict(fitted, new_data = newdat)
  expect_equal(nrow(preds), nrow(newdat))
})
