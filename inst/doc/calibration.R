## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4,
  out.width = "100%",
  dpi = 96
)

## ----setup, message = FALSE---------------------------------------------------
library(specProc)
library(dplyr)
library(tidyr)
library(ggplot2)

data("forageLIBS")
spectra_id <- forageLIBS |> select(1:2) |> names()
minerals <- forageLIBS |> select(3:14) |> names()
spectra <- forageLIBS |> select(-all_of(c(spectra_id, minerals)))

## ----intensities--------------------------------------------------------------
corrected <- baseline_arpls(spectra, lambda = 1e5, max.iter = 20)$correction

lines <- c(C = 247.856, K = 769.896, Fe = 259.940)
intensities <- line_intensities(corrected, lines) |>
  select(spectrum, line, intensity) |>
  pivot_wider(names_from = line, values_from = intensity)

data <- forageLIBS |>
  select(all_of(spectra_id), K, Fe) |>
  mutate(K_signal = intensities$K / intensities$C,
         Fe_signal = intensities$Fe / intensities$C)

## ----split--------------------------------------------------------------------
set.seed(1)
test_samples <- sample(unique(data$Sample), round(0.25 * n_distinct(data$Sample)))
calibration <- filter(data, !Sample %in% test_samples)
test <- filter(data, Sample %in% test_samples)
c(calibration = nrow(calibration), test = nrow(test))

## ----k-curve------------------------------------------------------------------
k_curve <- calibration_curve(calibration, K_signal, K)
k_curve

## ----k-numbers, include = FALSE-----------------------------------------------
k_fom <- k_curve$figures_of_merit
mandel <- k_curve$linearity |> filter(test == "Mandel")

## ----k-plot, fig.height = 4.5-------------------------------------------------
plot_calibration(k_curve)

## ----k-predict----------------------------------------------------------------
k_predicted <- predict(k_curve, test$K_signal)
head(k_predicted)
sqrt(mean((k_predicted$concentration - test$K)^2))   # RMSEP, % K

## ----fe-curves----------------------------------------------------------------
fe_calibration <- filter(calibration, !is.na(Fe))
fe_unweighted <- calibration_curve(fe_calibration, Fe_signal, Fe)
fe_weighted <- calibration_curve(fe_calibration, Fe_signal, Fe, weights = "1/x",
                                 lod_method = "intercept")

bind_rows(unweighted = fe_unweighted$figures_of_merit,
          weighted = fe_weighted$figures_of_merit, .id = "fit") |>
  select(fit, n, r_squared, lod, loq)

## ----fe-plot, fig.height = 4.5------------------------------------------------
plot_calibration(fe_weighted, interval = "prediction")

## ----fe-predict---------------------------------------------------------------
fe_test <- filter(test, !is.na(Fe))
fe_predicted <- predict(fe_weighted, fe_test$Fe_signal)
count(fe_predicted, below_loq)
median(abs(fe_predicted$concentration - fe_test$Fe) / fe_test$Fe)   # relative error

## ----pls----------------------------------------------------------------------
x <- as.matrix(corrected) / intensities$C
is_calibration <- !data$Sample %in% test_samples
x_calibration <- x[is_calibration, ]

cv <- pls::plsr(K ~ x, data = data.frame(K = data$K[is_calibration], x = I(x_calibration)),
                ncomp = 15, validation = "CV", segments = 10)
ncomp <- pls::selectNcomp(cv, method = "onesigma")
ncomp

## ----nas----------------------------------------------------------------------
repeated <- data$Sample[duplicated(data$Sample)]
differences <- sapply(repeated, \(s) {
  rows <- which(data$Sample == s)
  x[rows[1], ] - x[rows[2], ]
})
noise <- sd(differences) / sqrt(2)

k_nas <- nas(x_calibration, data$K[is_calibration], ncomp = ncomp, noise = noise)
k_nas$figures_of_merit

## ----nas-predict--------------------------------------------------------------
nas_predicted <- predict(k_nas, x[!is_calibration, ])
sqrt(mean((nas_predicted$predicted - test$K)^2))   # RMSEP, % K
summary(nas_predicted$selectivity)

## ----nas-numbers, include = FALSE---------------------------------------------
fom <- k_nas$figures_of_merit
rmsep_pls <- sqrt(mean((nas_predicted$predicted - test$K)^2))
rmsep_line <- sqrt(mean((k_predicted$concentration - test$K)^2))

## ----comparison-plot, fig.height = 4, fig.width = 7---------------------------
bind_rows(
  `K I 769.90 nm curve` = tibble(laboratory = test$K, predicted = k_predicted$concentration),
  `PLS (NAS)` = tibble(laboratory = test$K, predicted = nas_predicted$predicted),
  .id = "model"
) |>
  ggplot(aes(laboratory, predicted)) +
  geom_abline(colour = "grey50", linetype = "dashed") +
  geom_point(alpha = 0.6) +
  facet_wrap(~ model) +
  coord_equal() +
  labs(x = "Laboratory K (%)", y = "Predicted K (%)") +
  theme_bw()

