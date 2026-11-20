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
corrected <- baseline_arpls(spectra, lambda = 1e5, max_iter = 20)$correction

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
plot_calibration(k_curve, point_alpha = 0.4, xlab = "K (%)",
                 ylab = "K I 769.90 nm / C I 247.86 nm")

## ----k-predict----------------------------------------------------------------
k_predicted <- predict(k_curve, test$K_signal)
head(k_predicted)
sqrt(mean((k_predicted$concentration - test$K)^2))   # RMSEP, % K

## ----k-iso--------------------------------------------------------------------
k_iso <- calibration_curve(calibration, K_signal, K, lod_convention = "iso11843")
k_iso$figures_of_merit |> select(critical_value, lod, loq)

## ----k-york-------------------------------------------------------------------
noise <- line_intensities(corrected, lines[c("C", "K")]) |>
  mutate(relative_se = intensity_se / intensity) |>
  select(spectrum, line, relative_se) |>
  pivot_wider(names_from = line, values_from = relative_se)
data <- data |>
  mutate(K_se = 0.05 * K, K_signal_se = K_signal * sqrt(noise$K^2 + noise$C^2))
calibration <- filter(data, !Sample %in% test_samples)
test <- filter(data, Sample %in% test_samples)

k_york <- calibration_curve(calibration, K_signal, K, regression = "york",
                            concentration_se = K_se, signal_se = K_signal_se)
k_york

## ----york-numbers, include = FALSE--------------------------------------------
york_predicted <- predict(k_york, test$K_signal)
rmsep <- \(p) sqrt(mean((p$concentration - test$K)^2))
slope <- \(curve) signif(curve$coefficients$estimate[2], 3)

## ----fe-curves----------------------------------------------------------------
fe_calibration <- filter(calibration, !is.na(Fe))
fe_unweighted <- calibration_curve(fe_calibration, Fe_signal, Fe)
fe_weighted <- calibration_curve(fe_calibration, Fe_signal, Fe, weights = "1/x",
                                 lod_method = "intercept")

bind_rows(unweighted = fe_unweighted$figures_of_merit,
          weighted = fe_weighted$figures_of_merit, .id = "fit") |>
  select(fit, n, r_squared, lod, loq)

## ----fe-plot, fig.height = 4.5------------------------------------------------
plot_calibration(fe_weighted, interval = "prediction", point_alpha = 0.4, xlab = "Fe (ppm)",
                 ylab = "Fe I 259.94 nm / C I 247.86 nm")

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

## ----domain-------------------------------------------------------------------
domain <- applicability_domain(x_calibration, ncomp = ncomp)
domain
outside <- predict(domain, x[!is_calibration, ])
count(outside, outlier)

## ----domain-plot, fig.height = 4.5--------------------------------------------
plot_influence(outside)

## ----domain-numbers, include = FALSE------------------------------------------
errors <- abs(nas_predicted$predicted - test$K)

## ----secondary----------------------------------------------------------------
wl_x <- as.numeric(colnames(x))
response <- 1 + 0.3 * (wl_x - mean(wl_x)) / diff(range(wl_x))
secondary <- sweep(cbind(x[, -1], x[, ncol(x)]), 2, response, "*") + 0.002

## ----transfer-----------------------------------------------------------------
predict_k <- \(spectra) drop(predict(cv, newdata = data.frame(x = I(spectra)), ncomp = ncomp))
set.seed(2)
transfer_samples <- sample(which(is_calibration), 20)
secondary_test <- secondary[!is_calibration, ]

by_ds <- ds(x[transfer_samples, ], secondary[transfer_samples, ])
by_sst <- sst(x[transfer_samples, ], secondary[transfer_samples, ])
by_slope_bias <- slope_bias(predict_k(secondary[transfer_samples, ]), data$K[transfer_samples])

prediction_error <- \(predicted) sqrt(mean((predicted - test$K)^2))
transfer_errors <- c(
  primary = prediction_error(predict_k(x[!is_calibration, ])),
  secondary = prediction_error(predict_k(secondary_test)),
  ds = prediction_error(predict_k(predict(by_ds, secondary_test))),
  sst = prediction_error(predict_k(predict(by_sst, secondary_test))),
  slope_bias = prediction_error(predict(by_slope_bias, predict_k(secondary_test)))
)
round(transfer_errors, 2)   # RMSEP, % K

## ----drift, fig.height = 4----------------------------------------------------
carbon <- control_chart(intensities$C, order = data$Measurement, method = "ewma",
                        baseline = 100)
plot_control_chart(carbon, xlab = "Measurement (acquisition order)",
                   ylab = "EWMA of the C I 247.86 nm intensity")

## ----drift-numbers, include = FALSE-------------------------------------------
runs <- rle(carbon$signal)
longest <- which.max(ifelse(runs$values, runs$lengths, 0))
run_end <- cumsum(runs$lengths)[longest]
run_start <- run_end - runs$lengths[longest] + 1

