# Tests for engine registration -------------------------------------------

test_that("all six engines are registered with parsnip", {
  engines <- parsnip::show_engines("gaussian_process_spatial")
  registered <- engines$engine
  expect_true("gstat"      %in% registered)
  expect_true("fields"     %in% registered)
  expect_true("GPvecchia"  %in% registered)
  expect_true("spNNGP"     %in% registered)
  expect_true("PrestoGP"   %in% registered)
  expect_true("sdmTMB"     %in% registered)
})

test_that("model mode is regression for all engines", {
  engines <- parsnip::show_engines("gaussian_process_spatial")
  expect_true(all(engines$mode == "regression"))
})

test_that("set_engine('gstat') works without error", {
  spec <- gaussian_process_spatial() |> parsnip::set_engine("gstat")
  expect_equal(spec$engine, "gstat")
})

test_that("set_engine('fields') works without error", {
  spec <- gaussian_process_spatial() |> parsnip::set_engine("fields")
  expect_equal(spec$engine, "fields")
})

test_that("set_engine('GPvecchia') works without error", {
  spec <- gaussian_process_spatial() |> parsnip::set_engine("GPvecchia")
  expect_equal(spec$engine, "GPvecchia")
})

test_that("set_engine('spNNGP') works without error", {
  spec <- gaussian_process_spatial() |> parsnip::set_engine("spNNGP")
  expect_equal(spec$engine, "spNNGP")
})

test_that("set_engine('PrestoGP') works without error", {
  spec <- gaussian_process_spatial() |> parsnip::set_engine("PrestoGP")
  expect_equal(spec$engine, "PrestoGP")
})

test_that("set_engine('sdmTMB') works without error", {
  spec <- gaussian_process_spatial() |> parsnip::set_engine("sdmTMB")
  expect_equal(spec$engine, "sdmTMB")
})

test_that("tunable parameters are recognised", {
  skip_if_not_installed("parsnip")
  spec <- gaussian_process_spatial(
    covariance_function = parsnip::tune(),
    range               = parsnip::tune(),
    nugget              = parsnip::tune(),
    sill                = parsnip::tune()
  ) |>
    parsnip::set_engine("gstat")

  tun <- if ("tunable" %in% getNamespaceExports("parsnip")) {
    getExportedValue("parsnip", "tunable")(spec)
  } else {
    getFromNamespace("tunable.model_spec", "parsnip")(spec)
  }
  expect_true("covariance_function" %in% tun$name)
  expect_true("range"               %in% tun$name)
  expect_true("nugget"              %in% tun$name)
  expect_true("sill"                %in% tun$name)
})
