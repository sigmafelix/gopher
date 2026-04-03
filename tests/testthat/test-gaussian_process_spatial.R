# Tests for gaussian_process_spatial() model specification ---------------

test_that("gaussian_process_spatial() returns correct class", {
  spec <- gaussian_process_spatial()
  expect_s3_class(spec, "gaussian_process_spatial")
  expect_s3_class(spec, "model_spec")
})

test_that("gaussian_process_spatial() default mode is regression", {
  spec <- gaussian_process_spatial()
  expect_equal(spec$mode, "regression")
})

test_that("gaussian_process_spatial() stores args as quosures", {
  spec <- gaussian_process_spatial(
    covariance_function = "exponential",
    range  = 10,
    nugget = 0.1,
    sill   = 1.5
  )
  expect_true(rlang::is_quosure(spec$args$covariance_function))
  expect_true(rlang::is_quosure(spec$args$range))
  expect_true(rlang::is_quosure(spec$args$nugget))
  expect_true(rlang::is_quosure(spec$args$sill))
})

test_that("gaussian_process_spatial() args evaluate correctly", {
  spec <- gaussian_process_spatial(
    covariance_function = "spherical",
    range  = 50,
    nugget = 0.2,
    sill   = 2.0
  )
  expect_equal(rlang::eval_tidy(spec$args$covariance_function), "spherical")
  expect_equal(rlang::eval_tidy(spec$args$range),  50)
  expect_equal(rlang::eval_tidy(spec$args$nugget), 0.2)
  expect_equal(rlang::eval_tidy(spec$args$sill),   2.0)
})

test_that("gaussian_process_spatial() NULL args are NULL quosures", {
  spec <- gaussian_process_spatial()
  expect_true(rlang::quo_is_null(spec$args$covariance_function))
  expect_true(rlang::quo_is_null(spec$args$range))
})

test_that("gaussian_process_spatial() print method works", {
  spec <- gaussian_process_spatial() |> parsnip::set_engine("gstat")
  expect_output(print(spec), "Gaussian Process Model")
})

test_that("update.gaussian_process_spatial() changes args", {
  spec  <- gaussian_process_spatial(covariance_function = "exponential")
  spec2 <- update(spec, covariance_function = "matern", range = 20)
  expect_equal(
    rlang::eval_tidy(spec2$args$covariance_function), "matern"
  )
  expect_equal(rlang::eval_tidy(spec2$args$range), 20)
})

test_that("update fresh = TRUE resets unspecified args to NULL", {
  spec  <- gaussian_process_spatial(
    covariance_function = "exponential",
    range = 10
  )
  spec2 <- update(spec, covariance_function = "gaussian", fresh = TRUE)
  expect_equal(
    rlang::eval_tidy(spec2$args$covariance_function), "gaussian"
  )
  expect_true(rlang::quo_is_null(spec2$args$range))
})
