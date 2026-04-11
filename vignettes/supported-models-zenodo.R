## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE,
  warning = FALSE,
  message = FALSE
)


## -----------------------------------------------------------------------------
library(gopher)
library(parsnip)
library(sf)

# zenodo_url <- paste0(
#   "https://zenodo.org/records/10120281/files/",
#   "gridded_OSM_GSV.RDS?download=1"
# )

# data_path <- file.path(tempdir(), "gridded_OSM_GSV.RDS")

# if (!file.exists(data_path)) {
#   download.file(zenodo_url, destfile = data_path, mode = "wb")
# }

data_path <- "~/Downloads/gridded_gsv.rds"
gsv_sf <- readRDS(data_path)

stopifnot(inherits(gsv_sf, "sf"))
stopifnot(all(c("NO", "NO2") %in% names(gsv_sf)))

covariate_vars <- grep("^dens_", names(gsv_sf), value = TRUE)
response_vars <- c("NO", "NO2")

length(covariate_vars)
covariate_vars


## -----------------------------------------------------------------------------
gsv_pt <- gsv_sf |>
  sf::st_transform(32610) |>
  sf::st_point_on_surface()

xy <- sf::st_coordinates(gsv_pt)

gsv_df <- sf::st_drop_geometry(gsv_pt)
gsv_df$x <- xy[, "X"] / 1000
gsv_df$y <- xy[, "Y"] / 1000


## -----------------------------------------------------------------------------
set.seed(2026)

prepare_split <- function(data, response, covariates, prop = 0.8, max_n = 300) {
  keep <- stats::complete.cases(
    data[, c(response, covariates, "x", "y"), drop = FALSE]
  )
  data_complete <- data[keep, c(response, covariates, "x", "y")]
  data_complete <- data_complete[!duplicated(data_complete[, c("x", "y")]), ]

  if (nrow(data_complete) > max_n) {
    data_complete <- data_complete[
      sample.int(nrow(data_complete), size = max_n),
      ,
      drop = FALSE
    ]
  }

  n <- nrow(data_complete)
  train_id <- sample.int(n, size = floor(prop * n))

  list(
    data = data_complete,
    train = data_complete[train_id, , drop = FALSE],
    test = data_complete[-train_id, , drop = FALSE],
    formula = stats::reformulate(covariates, response = response)
  )
}

no_split <- prepare_split(gsv_df, "NO", covariate_vars)
no2_split <- prepare_split(gsv_df, "NO2", covariate_vars)


## -----------------------------------------------------------------------------
make_engine_specs <- function(train_data) {
  n_train <- nrow(train_data)

  list(
    gstat = gaussian_process_spatial(covariance_function = "matern") |>
      set_engine(
        "gstat",
        range = 500,
        psill = 120
      ),

    fields = gaussian_process_spatial(covariance_function = "matern") |>
      set_engine("fields"),

    GPvecchia = gaussian_process_spatial(
      covariance_function = NULL
    ) |>
      set_engine(
        "GPvecchia",
        m = max(50L, min(30L, n_train - 1L))
      ),

    spNNGP = gaussian_process_spatial(
      covariance_function = NULL,
      range = 2000
    ) |>
      set_engine(
        "spNNGP",
        covariance_function = "matern",
        n_neighbors = max(50L, min(15L, n_train - 1L)),
        n_samples = 1000L,
        n_burnin = 500L
      ),

    PrestoGP = gaussian_process_spatial(covariance_function = "matern") |>
      set_engine(
        "PrestoGP",
        n_neighbors = max(50L, min(20L, n_train - 1L)),
        model_type = "vecchia",
        quiet = TRUE
      ),

    sdmTMB = gaussian_process_spatial(covariance_function = "matern") |>
      set_engine(
        "sdmTMB",
        spatial = "on",
        n_knots = 80
      )
  )
}


## -----------------------------------------------------------------------------
engine_specs_no <- make_engine_specs(no_split$train)

fits_no <- lapply(engine_specs_no, function(spec) {
  fit(spec, no_split$formula, data = no_split$train)
})

preds_no <- lapply(fits_no, function(mod) {
  predict(mod, new_data = no_split$test)
})

results_no <- data.frame(
  engine = names(preds_no),
  rmse = vapply(
    preds_no,
    function(x) {
      sqrt(mean((no_split$test$NO - x$.pred)^2))
    },
    numeric(1)
  )
)

results_no[order(results_no$rmse), ]


## -----------------------------------------------------------------------------
pred_int_no_gstat <- predict(
  fits_no$sdmTMB,
  new_data = no_split$test,
  type = "pred_int"
)

head(pred_int_no_gstat)


## -----------------------------------------------------------------------------
engine_specs_no2 <- make_engine_specs(no2_split$train)

fits_no2 <- lapply(engine_specs_no2, function(spec) {
  fit(spec, no2_split$formula, data = no2_split$train)
})

preds_no2 <- lapply(fits_no2, function(mod) {
  predict(mod, new_data = no2_split$test)
})

results_no2 <- data.frame(
  engine = names(preds_no2),
  rmse = vapply(
    preds_no2,
    function(x) {
      sqrt(mean((no2_split$test$NO2 - x$.pred)^2))
    },
    numeric(1)
  )
)

results_no2[order(results_no2$rmse), ]


## -----------------------------------------------------------------------------
fit_all_supported_models <- function(data, response, covariates) {
  split <- prepare_split(data, response, covariates)
  specs <- make_engine_specs(split$train)

  fits <- lapply(specs, function(spec) {
    fit(spec, split$formula, data = split$train)
  })

  preds <- lapply(fits, function(mod) {
    predict(mod, new_data = split$test)
  })

  data.frame(
    response = response,
    engine = names(preds),
    rmse = vapply(
      preds,
      function(x) {
        sqrt(mean((split$test[[response]] - x$.pred)^2))
      },
      numeric(1)
    ),
    mae = vapply(
      preds,
      function(x) {
        mean(abs(split$test[[response]] - x$.pred))
      },
      numeric(1)
    ),
    row.names = NULL
  )
}

benchmark_table <- do.call(
  rbind,
  lapply(response_vars, function(resp) {
    fit_all_supported_models(gsv_df, resp, covariate_vars)
  })
)


benchmark_table[order(benchmark_table$response, benchmark_table$rmse), ]


## -----------------------------------------------------------------------------
# plot benchmark results
library(ggplot2)
ggplot(benchmark_table, aes(x = engine, y = rmse, fill = engine)) +
  geom_point() +
  geom_pointrange(aes(ymin = 0, ymax = rmse), size = 0.5) +
  facet_wrap(~ response, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = "RMSE of Supported GP Engines on the Zenodo Street-View Density Data",
    x = "GP Engine",
    y = "Root Mean Squared Error (RMSE)"
  )

