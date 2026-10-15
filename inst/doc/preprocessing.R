## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4,
  out.width = "100%",
  dpi = 96
)

## ----requirements, include = FALSE--------------------------------------------
needed <- c("recipes", "parsnip", "plsmod", "mixOmics", "workflows", "tune", "rsample", "yardstick")
if (!all(vapply(needed, requireNamespace, logical(1), quietly = TRUE))) {
  knitr::knit_exit("This vignette needs the tidymodels packages recipes, parsnip, plsmod (with mixOmics), workflows, tune, rsample and yardstick.")
}

## ----setup, message = FALSE---------------------------------------------------
suppressPackageStartupMessages(library(mixOmics))
library(specProc)
library(recipes)
library(parsnip)
library(plsmod)
library(workflows)
library(tune)
library(rsample)
library(yardstick)
library(dplyr)
library(purrr)
library(ggplot2)
# mixOmics (the PLS engine) and MASS, which it attaches, have functions
# named like those of the tidyverse and tidymodels: use the latter
select <- dplyr::select
map <- purrr::map
tune <- tune::tune

data("forageLIBS")
spectra_id <- forageLIBS |> select(1:2) |> names()
minerals <- forageLIBS |> select(3:14) |> names()
spectra <- forageLIBS |> select(-all_of(c(spectra_id, minerals)))
wl <- as.numeric(names(spectra))

## ----recipe-------------------------------------------------------------------
base <- recipe(K ~ ., data = forageLIBS) |>
  update_role(all_of(spectra_id), all_of(setdiff(minerals, "K")), new_role = "id")

## ----baseline, fig.height = 3.5-----------------------------------------------
region <- wl > 276 & wl < 290
fitted_baselines <- map(c(`1e3` = 1e3, `1e5` = 1e5, `1e7` = 1e7), \(lambda) {
  unlist(baseline_arpls(spectra[1, ], lambda = lambda, max.iter = 20)$background)[region]
})
first_spectrum <- tibble(wavelength = wl[region], counts = unlist(spectra[1, region]))

tibble(wavelength = wl[region], !!!fitted_baselines) |>
  tidyr::pivot_longer(-wavelength, names_to = "lambda", values_to = "counts") |>
  ggplot(aes(wavelength, counts)) +
  geom_line(data = first_spectrum, colour = "grey55") +
  geom_line(aes(colour = lambda), linewidth = 0.8) +
  coord_cartesian(ylim = c(800, 2500)) +
  labs(x = "Wavelength (nm)", y = "Counts (lines truncated)") +
  theme_bw() +
  theme(legend.position = "top")

## ----step-baseline------------------------------------------------------------
corrected <- base |>
  step_baseline(all_predictors(), lambda = 1e5, options = list(max.iter = 20))

## ----corrected-data-----------------------------------------------------------
corrected_data <- corrected |> prep() |> bake(new_data = NULL)
base_corrected <- recipe(K ~ ., data = corrected_data) |>
  update_role(all_of(spectra_id), all_of(setdiff(minerals, "K")), new_role = "id")

## ----normalizations-----------------------------------------------------------
lines <- c(K = 769.896, Mg = 285.213, Ca = 317.933, Na = 589.592, Mn = 257.610, Fe = 259.940)
normalizations <- list(
  none = identity,
  area = \(r) step_spectral_norm(r, all_predictors(), method = "area"),
  SNV = \(r) step_snv(r, all_predictors()),
  MSC = \(r) step_msc(r, all_predictors()),
  `C I 247.86` = \(r) step_line_ratio(r, all_predictors(), reference = 247.856, window = 0.15)
)

imap(normalizations, \(normalize, name) {
  intensities <- base_corrected |>
    normalize() |>
    step_line_intensities(all_predictors(), lines = lines) |>
    prep() |>
    bake(new_data = NULL)
  map_dbl(names(lines), \(e) cor(intensities[[paste0("line_", e)]], forageLIBS[[e]],
                                 use = "complete.obs")) |>
    set_names(names(lines))
}) |>
  bind_rows(.id = "normalization") |>
  mutate(across(-normalization, \(r) round(r, 2)))

## ----wavelet------------------------------------------------------------------
base_corrected |>
  step_wavelet(all_predictors(), level = 3) |>
  prep() |>
  bake(new_data = NULL) |>
  select(starts_with("wav_")) |>
  ncol()

## ----split--------------------------------------------------------------------
set.seed(1)
split <- group_initial_split(corrected_data, group = Sample, prop = 0.8)
training_set <- training(split)
folds <- group_vfold_cv(training_set, group = Sample, v = 5)

## ----recipes------------------------------------------------------------------
carbon_step <- \(r) step_line_ratio(r, all_predictors(), reference = 247.856, window = 0.15)

recipes <- list(
  SNV = base_corrected |> step_snv(all_predictors()),
  carbon = base_corrected |> carbon_step(),
  `carbon + wavelet` = base_corrected |> carbon_step() |>
    step_wavelet(all_predictors(), level = 3)
)

# mixOmics scales every predictor by default; scale = FALSE only centers
pls_model <- pls(num_comp = tune()) |>
  set_mode("regression") |>
  set_engine("mixOmics", scale = FALSE)

## ----tune---------------------------------------------------------------------
tuned <- map(recipes, \(rec) {
  tune_grid(workflow(rec, pls_model), resamples = folds,
            grid = tibble(num_comp = c(2, 4, 6, 8, 10, 12, 15)),
            metrics = metric_set(rmse, rsq))
})

map(tuned, \(res) show_best(res, metric = "rmse", n = 1)) |>
  bind_rows(.id = "recipe") |>
  select(recipe, num_comp, rmse = mean, std_err)

## ----tune-plot, fig.height = 3.5----------------------------------------------
map(tuned, collect_metrics) |>
  bind_rows(.id = "recipe") |>
  filter(.metric == "rmse") |>
  ggplot(aes(num_comp, mean, colour = recipe)) +
  geom_line() +
  geom_point() +
  labs(x = "PLS components", y = "Cross-validated RMSE (% K)", colour = NULL) +
  theme_bw()

## ----best, include = FALSE----------------------------------------------------
scores <- map_dbl(tuned, \(res) show_best(res, metric = "rmse", n = 1)$mean)
best_recipe <- names(which.min(scores))

## ----final--------------------------------------------------------------------
final <- workflow(recipes[[best_recipe]], pls_model) |>
  finalize_workflow(select_best(tuned[[best_recipe]], metric = "rmse")) |>
  last_fit(split, metrics = metric_set(rmse, rsq))
collect_metrics(final)

## ----test-numbers, include = FALSE--------------------------------------------
test_rmse <- collect_metrics(final)$.estimate[1]
cv_rmse <- min(scores)

## ----final-plot, fig.height = 4, fig.width = 4.5, out.width = "60%"-----------
collect_predictions(final) |>
  ggplot(aes(K, .pred)) +
  geom_abline(colour = "grey50", linetype = "dashed") +
  geom_point(alpha = 0.6) +
  coord_equal() +
  labs(x = "Laboratory K (%)", y = "Predicted K (%)") +
  theme_bw()

