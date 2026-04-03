# Helper: synthetic spatial data for testing
#
# All test data is generated in this file and reused across test files.
# Functions are prefixed with `make_` to make them clearly test-helpers.

# ---- Purely spatial data ------------------------------------------------

#' Create a small synthetic spatial dataset (data.frame with x/y)
make_spatial_df <- function(n = 40L, seed = 42L) {
  set.seed(seed)
  x   <- runif(n, 0, 10)
  y   <- runif(n, 0, 10)
  # Spatial signal + noise
  z   <- sin(x) + cos(y) + rnorm(n, 0, 0.3)
  cov1 <- x * 0.5 + rnorm(n, 0, 0.1)   # a covariate
  data.frame(x = x, y = y, z = z, cov1 = cov1)
}

#' Convert the spatial df to sf
make_spatial_sf <- function(n = 40L, seed = 42L) {
  df <- make_spatial_df(n, seed)
  sf::st_as_sf(df, coords = c("x", "y"), crs = 4326L)
}

#' New-data sf for prediction (held-out grid)
make_pred_sf <- function(n_side = 5L) {
  grid <- expand.grid(
    x = seq(1, 9, length.out = n_side),
    y = seq(1, 9, length.out = n_side)
  )
  grid$cov1 <- runif(nrow(grid), 0, 5)
  sf::st_as_sf(grid, coords = c("x", "y"), crs = 4326L)
}

#' New-data data.frame for prediction
make_pred_df <- function(n_side = 5L) {
  grid <- expand.grid(
    x = seq(1, 9, length.out = n_side),
    y = seq(1, 9, length.out = n_side)
  )
  grid$cov1 <- runif(nrow(grid), 0, 5)
  grid
}

# ---- Spatiotemporal data ------------------------------------------------

#' Minimal spatiotemporal dataset
make_st_df <- function(n_loc = 10L, n_time = 3L, seed = 99L) {
  set.seed(seed)
  locs  <- data.frame(
    x = runif(n_loc, 0, 10),
    y = runif(n_loc, 0, 10)
  )
  times <- seq(as.Date("2020-01-01"), by = "month", length.out = n_time)
  df    <- do.call(rbind, lapply(times, function(t) {
    row <- locs
    row$z    <- sin(locs$x) + cos(locs$y) + rnorm(n_loc, 0, 0.3)
    row$time <- t
    row
  }))
  sf::st_as_sf(df, coords = c("x", "y"), crs = 4326L)
}
