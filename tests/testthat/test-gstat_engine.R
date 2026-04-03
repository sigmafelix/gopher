# Tests for the gstat engine -------------------------------------------

skip_if_not_installed("gstat")

# ---- gstat_gp_fit ---------------------------------------------------

test_that("gstat_gp_fit returns gopher_gstat_fit object (sf input)", {
  train <- make_spatial_sf(30L)
  fit   <- gstat_gp_fit(z ~ 1, data = train)
  expect_s3_class(fit, "gopher_gstat_fit")
  expect_false(fit$is_spatiotemporal)
})

test_that("gstat_gp_fit returns gopher_gstat_fit object (data.frame input)", {
  train <- make_spatial_df(30L)
  fit   <- gstat_gp_fit(z ~ 1, data = train)
  expect_s3_class(fit, "gopher_gstat_fit")
})

test_that("gstat_gp_fit stores empirical variogram", {
  train <- make_spatial_sf(30L)
  fit   <- gstat_gp_fit(z ~ 1, data = train)
  expect_true(!is.null(fit$empirical_variogram))
})

test_that("gstat_gp_fit stores fitted variogram model", {
  train <- make_spatial_sf(30L)
  fit   <- gstat_gp_fit(z ~ 1, data = train)
  expect_s3_class(fit$variogram_fit, "variogramModel")
})

test_that("gstat_gp_fit accepts all covariance functions", {
  train   <- make_spatial_sf(30L)
  cov_fns <- c("exponential", "spherical", "gaussian", "matern")
  for (cf in cov_fns) {
    fit <- gstat_gp_fit(z ~ 1, data = train, covariance_function = cf)
    expect_s3_class(fit, "gopher_gstat_fit",
                    info = paste("covariance_function =", cf))
  }
})

test_that("gstat_gp_fit with manual parameters skips fitting when fit_variogram=FALSE", {
  train <- make_spatial_sf(30L)
  fit   <- gstat_gp_fit(
    z ~ 1, data = train,
    range = 3, nugget = 0.05, sill = 0.8,
    fit_variogram = FALSE
  )
  # Variogram should use our supplied values
  expect_equal(fit$variogram_fit$range[2], 3, tolerance = 1e-6)
})

test_that("gstat_gp_fit supports universal kriging formula", {
  train <- make_spatial_sf(30L)
  fit   <- gstat_gp_fit(z ~ cov1, data = train)
  expect_s3_class(fit, "gopher_gstat_fit")
})

# ---- gstat_gp_predict -----------------------------------------------

test_that("gstat_gp_predict returns tibble with .pred column", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_sf()
  fit    <- gstat_gp_fit(z ~ 1, data = train)
  preds  <- gstat_gp_predict(fit, new_data = newdat)
  expect_s3_class(preds, "tbl_df")
  expect_named(preds, ".pred")
  expect_equal(nrow(preds), nrow(newdat))
  expect_true(all(is.finite(preds$.pred)))
})

test_that("gstat_gp_predict with type='pred_int' returns interval columns", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_sf()
  fit    <- gstat_gp_fit(z ~ 1, data = train)
  preds  <- gstat_gp_predict(fit, new_data = newdat, type = "pred_int")
  expect_named(preds, c(".pred", ".pred_lower", ".pred_upper"))
  expect_true(all(preds$.pred_lower <= preds$.pred + 1e-9))
  expect_true(all(preds$.pred_upper >= preds$.pred - 1e-9))
})

test_that("gstat_gp_predict works with data.frame new_data", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_df()
  fit    <- gstat_gp_fit(z ~ 1, data = train)
  preds  <- gstat_gp_predict(fit, new_data = newdat)
  expect_equal(nrow(preds), nrow(newdat))
})

test_that("gstat_gp_predict works with universal kriging", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_sf()
  fit    <- gstat_gp_fit(z ~ cov1, data = train)
  preds  <- gstat_gp_predict(fit, new_data = newdat)
  expect_equal(nrow(preds), nrow(newdat))
})

# ---- parsnip integration --------------------------------------------

test_that("parsnip fit() + predict() pipeline works with gstat", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_sf()

  spec <- gaussian_process_spatial(covariance_function = "exponential") |>
    parsnip::set_engine("gstat")

  fitted <- parsnip::fit(spec, z ~ 1, data = train)
  expect_s3_class(fitted, "model_fit")

  preds <- predict(fitted, new_data = newdat)
  expect_named(preds, ".pred")
  expect_equal(nrow(preds), nrow(newdat))
})

test_that("parsnip predict(type='pred_int') pipeline works with gstat", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_sf()

  spec <- gaussian_process_spatial() |> parsnip::set_engine("gstat")
  fitted <- parsnip::fit(spec, z ~ 1, data = train)
  preds  <- predict(fitted, new_data = newdat, type = "pred_int")
  expect_true(".pred_lower" %in% names(preds))
  expect_true(".pred_upper" %in% names(preds))
})

test_that("parsnip universal kriging pipeline works with gstat", {
  train  <- make_spatial_sf(30L)
  newdat <- make_pred_sf()

  spec   <- gaussian_process_spatial() |> parsnip::set_engine("gstat")
  fitted <- parsnip::fit(spec, z ~ cov1, data = train)
  preds  <- predict(fitted, new_data = newdat)
  expect_equal(nrow(preds), nrow(newdat))
})

test_that("gstat engine works with non-sf data in parsnip pipeline", {
  train  <- make_spatial_df(30L)
  newdat <- make_pred_df()

  spec   <- gaussian_process_spatial() |> parsnip::set_engine("gstat")
  fitted <- parsnip::fit(spec, z ~ 1, data = train)
  preds  <- predict(fitted, new_data = newdat)
  expect_equal(nrow(preds), nrow(newdat))
})
