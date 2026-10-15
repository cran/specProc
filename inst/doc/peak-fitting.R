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
if (!requireNamespace("purrr", quietly = TRUE)) {
  knitr::knit_exit("This vignette needs the purrr package, which is not installed.")
}

## ----setup, message = FALSE---------------------------------------------------
library(specProc)
library(dplyr)
library(purrr)
library(ggplot2)

data("forageLIBS")
spectra_id <- forageLIBS |> select(1:2) |> names()
minerals <- forageLIBS |> select(3:14) |> names()
spectra <- forageLIBS |> select(-all_of(c(spectra_id, minerals)))
wl <- as.numeric(names(spectra))

## ----spectra, fig.height = 3.5------------------------------------------------
forageLIBS |>
  slice(1:3) |>
  select(Measurement, all_of(names(spectra)[wl > 400 & wl < 410])) |>
  plot_spectra(id = Measurement)

## ----identify, eval = FALSE---------------------------------------------------
# candidates <- libs_lines(c("K I", "Ca I", "Fe I", "Mn I"), wavelength = c(400, 410))
# plot_lines(colMeans(spectra), candidates)
# line_finder(spectra)

## ----calibration--------------------------------------------------------------
first_detector <- spectra[seq_len(which(diff(wl) < 0)[1])]

reference <- c(`C I` = 247.856, `Mg I` = 285.213, `Ca II` = 317.933, `Ca II` = 370.603,
               `Ca II` = 373.690, `K I` = 404.414, `Ca I` = 430.253, `Ca I` = 445.478,
               `Mg I` = 516.732, `Mg I` = 517.268, `Mg I` = 518.360, `Na I` = 588.995,
               `Na I` = 589.592, `Ca I` = 616.217, `Ca I` = 643.907)
wavecal <- wavelength_calibration(first_detector, reference)
wavecal$segments

## ----calibration-plot, fig.height = 3.5---------------------------------------
plot_wavelength_calibration(wavecal)

## ----corrected----------------------------------------------------------------
calibrated <- apply_calibration(first_detector, wavecal)
corrected <- baseline_arpls(calibrated, lambda = 1e5, max.iter = 20)$correction
wl_corrected <- as.numeric(names(corrected))
mean_spectrum <- colMeans(corrected)

# The channels between two wavelengths, as a one-row table
window <- function(from, to, x = mean_spectrum) {
  keep <- wl_corrected > from & wl_corrected < to
  if (is.null(dim(x))) as_tibble(as.list(x[keep])) else x[keep]
}

## ----halpha-------------------------------------------------------------------
halpha <- window(653.5, 659.5)
fits <- c("gaussian", "lorentzian", "pseudo_voigt", "voigt") |>
  set_names() |>
  map(\(p) peak_fit(halpha, profile = p))

tibble(
  profile = names(fits),
  AIC = map_dbl(fits, \(f) AIC(f$fit[[1]])),
  area = map_dbl(fits, \(f) coef(f$fit[[1]])[["A"]])
) |>
  mutate(across(c(AIC, area), round))

## ----halpha-voigt-------------------------------------------------------------
fits$voigt$tidied[[1]]

## ----halpha-plot, fig.height = 5----------------------------------------------
plot_fit(fits$lorentzian, title = "H-alpha 656.28 nm, Lorentzian profile")

## ----doublet------------------------------------------------------------------
doublet <- window(404.1, 405.0)
k_voigt <- multipeak_fit(doublet, peaks = c(404.414, 404.721), profiles = "voigt")
k_voigt$tidied[[1]]

## ----doublet-gaussian---------------------------------------------------------
k_gauss <- multipeak_fit(doublet, peaks = c(404.414, 404.721), profiles = "gaussian")
k_gauss$tidied[[1]]
c(voigt = AIC(k_voigt$fit[[1]]), gaussian = AIC(k_gauss$fit[[1]]))

## ----doublet-plot, fig.height = 5---------------------------------------------
plot_fit(k_gauss, title = "K I 404.41 and 404.72 nm, Gaussian profiles")

## ----batch, warning = FALSE---------------------------------------------------
doublets <- bind_cols(
  select(forageLIBS, Measurement),
  window(404.1, 405.0, x = corrected)
)
k_fits <- multipeak_fit(doublets, peaks = c(404.414, 404.721), profiles = "gaussian",
                        id = "Measurement")

k_areas <- k_fits |>
  mutate(area = map_dbl(tidied, \(t) if (is.null(t)) NA else t$estimate[t$term == "A_1"])) |>
  select(Measurement, area) |>
  left_join(select(forageLIBS, Measurement, K), by = "Measurement")
sum(is.na(k_areas$area))

## ----batch-plot, fig.height = 3.5, fig.width = 5, out.width = "70%", warning = FALSE----
ggplot(k_areas, aes(K, area)) +
  geom_point(alpha = 0.5) +
  labs(x = "K (%)", y = "Area of K I 404.41 nm") +
  theme_bw()

## ----intensities--------------------------------------------------------------
lines <- c(`K I 404.41` = 404.414, `Mg I 518.36` = 518.360, `Fe II 259.94` = 259.940)
line_intensities(corrected[1:2, ], lines, half_width = 0.1)

