# Tests for utility functions --------------------------------------------

# ---- extract_coords -------------------------------------------------

test_that("extract_coords works with sf object", {
  sf_data <- make_spatial_sf(10L)
  coords  <- gopher:::extract_coords(sf_data)
  expect_true(is.matrix(coords))
  expect_equal(ncol(coords), 2L)
  expect_equal(nrow(coords), 10L)
  expect_true(all(c("X", "Y") %in% colnames(coords)))
})

test_that("extract_coords works with x/y columns", {
  df     <- make_spatial_df(10L)
  coords <- gopher:::extract_coords(df)
  expect_true(is.matrix(coords))
  expect_equal(ncol(coords), 2L)
  expect_equal(colnames(coords), c("X", "Y"))
})

test_that("extract_coords works with explicit coord_cols", {
  df        <- make_spatial_df(10L)
  names(df)[names(df) == "x"] <- "lon"
  names(df)[names(df) == "y"] <- "lat"
  coords    <- gopher:::extract_coords(df, coord_cols = c("lon", "lat"))
  expect_equal(colnames(coords), c("X", "Y"))
})

test_that("extract_coords errors when columns are missing", {
  df <- data.frame(a = 1:3, b = 4:6)
  expect_error(gopher:::extract_coords(df), "Cannot determine coordinate")
})

test_that("extract_coords errors when coord_cols length != 2", {
  df <- make_spatial_df(5L)
  expect_error(gopher:::extract_coords(df, coord_cols = "x"), "length 2")
})

# ---- extract_st_coords ----------------------------------------------

test_that("extract_st_coords appends numeric time from Date column", {
  df <- make_spatial_df(8L)
  df$day <- as.Date("2020-01-01") + seq_len(nrow(df))
  coords <- gopher:::extract_st_coords(df, time_col = "day")
  expect_true(is.matrix(coords))
  expect_equal(ncol(coords), 3L)
  expect_equal(colnames(coords), c("X", "Y", "T"))
})

test_that("extract_st_coords applies time scaling", {
  df <- make_spatial_df(5L)
  df$t <- as.numeric(seq_len(nrow(df)))
  coords1 <- gopher:::extract_st_coords(df, time_col = "t", time_scale = 1)
  coords2 <- gopher:::extract_st_coords(df, time_col = "t", time_scale = 2)
  expect_equal(coords2[, "T"], coords1[, "T"] / 2)
})

test_that("extract_st_coords errors for missing time column", {
  df <- make_spatial_df(4L)
  expect_error(
    gopher:::extract_st_coords(df, time_col = "not_a_col"),
    "Time column"
  )
})

# ---- translate_covariance -------------------------------------------

test_that("translate_covariance maps canonical names to gstat names", {
  expect_equal(gopher:::translate_covariance("exponential", "gstat"), "Exp")
  expect_equal(gopher:::translate_covariance("spherical",   "gstat"), "Sph")
  expect_equal(gopher:::translate_covariance("gaussian",    "gstat"), "Gau")
  expect_equal(gopher:::translate_covariance("matern",      "gstat"), "Mat")
  expect_equal(gopher:::translate_covariance("stein_matern","gstat"), "Ste")
})

test_that("translate_covariance maps canonical names to fields names", {
  expect_equal(
    gopher:::translate_covariance("exponential", "fields"), "Exponential"
  )
  expect_equal(
    gopher:::translate_covariance("matern", "fields"), "Matern"
  )
})

test_that("translate_covariance returns NULL default for NULL input", {
  result <- gopher:::translate_covariance(NULL, "gstat", default = NULL)
  expect_null(result)
})

test_that("translate_covariance returns specified default for NULL input", {
  result <- gopher:::translate_covariance(NULL, "gstat", default = "Exp")
  expect_equal(result, "Exp")
})

test_that("translate_covariance warns and falls back on unknown name", {
  expect_warning(
    res <- gopher:::translate_covariance("unknown_cov", "gstat", default = "Exp"),
    regexp = "not recognised"
  )
  expect_equal(res, "Exp")
})

# ---- parse_formula --------------------------------------------------

test_that("parse_formula extracts response correctly", {
  df  <- make_spatial_df(10L)
  out <- gopher:::parse_formula(z ~ 1, df)
  expect_equal(out$response,      df$z)
  expect_equal(out$response_name, "z")
  expect_false(out$has_covariates)
})

test_that("parse_formula detects covariates", {
  df  <- make_spatial_df(10L)
  out <- gopher:::parse_formula(z ~ cov1, df)
  expect_true(out$has_covariates)
  expect_true("cov1" %in% colnames(out$X))
})

test_that("parse_formula model matrix includes intercept for ~ 1", {
  df  <- make_spatial_df(10L)
  out <- gopher:::parse_formula(z ~ 1, df)
  expect_true("(Intercept)" %in% colnames(out$X))
})

# ---- drop_geometry --------------------------------------------------

test_that("drop_geometry removes geometry from sf", {
  sf_data <- make_spatial_sf(5L)
  plain   <- gopher:::drop_geometry(sf_data)
  expect_false(inherits(plain, "sf"))
  expect_false("geometry" %in% names(plain))
})

test_that("drop_geometry is a no-op for plain data.frame", {
  df    <- make_spatial_df(5L)
  plain <- gopher:::drop_geometry(df)
  expect_identical(plain, df)
})
