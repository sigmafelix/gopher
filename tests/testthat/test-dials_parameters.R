# Tests for dials parameter functions ------------------------------------

skip_if_not_installed("dials")

test_that("covariance_function() returns a dials parameter", {
  param <- covariance_function()
  expect_s3_class(param, "param")
  expect_equal(param$type, "character")
})

test_that("covariance_function() includes expected values", {
  param  <- covariance_function()
  vals   <- param$values
  expect_true("exponential" %in% vals)
  expect_true("spherical"   %in% vals)
  expect_true("gaussian"    %in% vals)
  expect_true("matern"      %in% vals)
})

test_that("gp_range() returns a dials parameter", {
  param <- gp_range()
  expect_s3_class(param, "param")
  expect_equal(param$type, "double")
})

test_that("gp_nugget() returns a dials parameter", {
  param <- gp_nugget()
  expect_s3_class(param, "param")
  expect_equal(param$type, "double")
})

test_that("gp_sill() returns a dials parameter", {
  param <- gp_sill()
  expect_s3_class(param, "param")
  expect_equal(param$type, "double")
})

test_that("gp_range() lower bound is positive", {
  param <- gp_range()
  bounds <- param$range
  expect_gt(bounds[[1]], 0)
})
