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
library(tibble)
library(ggplot2)

data("forageLIBS")
spectra_id <- forageLIBS |> select(1:2) |> names()
minerals <- forageLIBS |> select(3:14) |> names()
spectra <- forageLIBS |> select(-all_of(c(spectra_id, minerals)))
wl <- as.numeric(names(spectra))

corrected <- baseline_arpls(spectra, lambda = 1e5, max.iter = 20)$correction
mean_spectrum <- colMeans(corrected)

## ----saturation---------------------------------------------------------------
saturated <- saturation_summary(spectra, limit = 65535, tolerance = 655)

# adjacent saturated channels belong to the same line
saturated$channels |>
  mutate(line = cumsum(c(1, diff(index) > 2))) |>
  group_by(line) |>
  summarise(wavelength = round(mean(wavelength), 2), spectra = max(n_spectra)) |>
  select(-line) |>
  arrange(desc(spectra))

## ----density------------------------------------------------------------------
halpha <- as_tibble(as.list(mean_spectrum[wl > 653.5 & wl < 659.5])) |>
  peak_fit(profile = "lorentzian")
width <- halpha$tidied[[1]] |> filter(term == "wL")
width
ne <- electron_density(width$estimate, method = "halpha")
ne

## ----nist, eval = FALSE-------------------------------------------------------
# ca_lines <- nist_lines(c("Ca I", "Ca II"), wavelength = c(310, 660))
# nist_ionization_energy("Ca I")

## ----atomic-------------------------------------------------------------------
atomic <- tribble(
  ~species, ~stage, ~wavelength, ~Aki,    ~accuracy, ~Ek,       ~gk,
  "Ca I",   1,      428.301,     4.34e7,  "C+",      4.779784,  5,
  "Ca I",   1,      430.253,     1.36e8,  "C+",      4.779784,  5,
  "Ca I",   1,      431.865,     7.40e7,  "C+",      4.7690284, 3,
  "Ca I",   1,      443.496,     6.70e7,  "C",       4.680635,  5,
  "Ca I",   1,      445.478,     8.70e7,  "C",       4.681327,  7,
  "Ca I",   1,      558.876,     4.90e7,  "D",       4.7435268, 7,
  "Ca I",   1,      610.272,     9.60e6,  "C",       3.910399,  3,
  "Ca I",   1,      612.222,     2.87e7,  "C",       3.910399,  3,
  "Ca I",   1,      616.217,     4.77e7,  "C",       3.910399,  3,
  "Ca I",   1,      643.907,     5.30e7,  "D",       4.450647,  9,
  "Ca I",   1,      646.257,     4.70e7,  "D",       4.4409544, 7,
  "Ca I",   1,      649.378,     4.40e7,  "D",       4.4300117, 5,
  "Ca II",  2,      315.887,     3.10e8,  "C",       7.047169,  4,
  "Ca II",  2,      317.933,     3.60e8,  "C",       7.049551,  6,
  "Ca II",  2,      370.603,     8.80e7,  "C",       6.467875,  2,
  "Ca II",  2,      373.690,     1.70e8,  "C",       6.467875,  2
)
ionization_energy <- 6.1131549   # Ca I, eV

## ----intensities--------------------------------------------------------------
lines <- line_intensities(mean_spectrum, atomic)
lines |> select(species, wavelength, peak_wavelength, intensity, snr)

## ----boltzmann----------------------------------------------------------------
neutral <- boltzmann(filter(lines, stage == 1))
neutral

## ----saha---------------------------------------------------------------------
saha <- saha_boltzmann(lines, ionization_energy = ionization_energy, electron_density = ne)
saha

## ----saha-plot, fig.height = 4------------------------------------------------
plot_boltzmann(saha)

## ----mcwhirter----------------------------------------------------------------
lte <- mcwhirter_criterion(temperature = saha$temperature, delta_e = 3.15, electron_density = ne)
lte

## ----self-absorption----------------------------------------------------------
k_lines <- tribble(
  ~species, ~wavelength, ~Aki,      ~Ei,      ~Ek,      ~gk,
  "K I",    404.414,     1.150e+06, 0,        3.064907, 4,
  "K I",    404.721,     1.070e+06, 0,        3.062581, 2,
  "K I",    691.108,     2.500e+06, 1.609958, 3.403454, 2,
  "K I",    693.877,     4.956e+06, 1.617113, 3.403454, 2,
  "K I",    766.490,     3.779e+07, 0,        1.617113, 4,
  "K I",    769.896,     3.734e+07, 0,        1.609958, 2
)
line_intensities(mean_spectrum, k_lines) |>
  correct_self_absorption(temperature = saha$temperature) |>
  select(wavelength, measured_intensity, SA, reference)

## ----cf-fetch, eval = FALSE---------------------------------------------------
# species <- c("Ca I", "Ca II", "Mg I", "Mg II", "K I", "K II", "Na I", "Na II")
# levels <- nist_levels(species)
# ionization <- nist_ionization_energy(c("Ca I", "Mg I", "K I", "Na I"))

## ----cf-atomic----------------------------------------------------------------
cf_atomic <- tribble(
  ~species, ~wavelength, ~Aki,      ~accuracy, ~Ek,      ~gk,
  "Ca I",   428.301,     4.340e+07, "C+",      4.779784, 5,
  "Ca I",   430.253,     1.360e+08, "C+",      4.779784, 5,
  "Ca I",   431.865,     7.400e+07, "C+",      4.769028, 3,
  "Ca I",   443.496,     6.700e+07, "C",       4.680635, 5,
  "Ca I",   445.478,     8.700e+07, "C",       4.681327, 7,
  "Ca I",   616.217,     4.770e+07, "C",       3.910399, 3,
  "Ca II",  315.887,     3.100e+08, "C",       7.047169, 4,
  "Ca II",  317.933,     3.600e+08, "C",       7.049551, 6,
  "Ca II",  370.603,     8.800e+07, "C",       6.467875, 2,
  "Ca II",  373.690,     1.700e+08, "C",       6.467875, 2,
  "Mg I",   383.829,     1.610e+08, "B+",      5.945916, 7,
  "Mg I",   516.732,     1.130e+07, "B+",      5.107827, 3,
  "Mg I",   517.268,     3.370e+07, "B+",      5.107827, 3,
  "Mg I",   518.360,     5.610e+07, "A",       5.107827, 3,
  "Mg II",  279.078,     4.010e+08, "A",       8.863762, 4,
  "Mg II",  279.800,     4.790e+08, "A",       8.863654, 6,
  "Mg II",  292.863,     1.150e+08, "A",       8.654711, 2,
  "Mg II",  293.651,     2.300e+08, "A",       8.654711, 2,
  "K I",    404.414,     1.150e+06, "A",       3.064907, 4,
  "K I",    404.721,     1.070e+06, "B+",      3.062581, 2,
  "K I",    691.108,     2.500e+06, "AA",      3.403454, 2,
  "K I",    693.877,     4.956e+06, "AA",      3.403454, 2,
  "Na I",   818.326,     4.290e+07, "A+",      3.616977, 4,
  "Na I",   819.482,     5.140e+07, "A+",      3.616971, 6
)
cf_ionization <- c(Ca = 6.113155, Mg = 7.646236, K = 4.340664, Na = 5.139077)

partition_table <- tribble(
  ~species, ~`6000`, ~`7000`, ~`8000`, ~`9000`, ~`10000`, ~`11000`, ~`12000`,
  "Ca I",   1.407,   1.798,   2.401,   3.281,   4.499,    6.184,    8.356,
  "Ca II",  2.389,   2.633,   2.917,   3.229,   3.560,    3.905,    4.262,
  "Mg I",   1.049,   1.109,   1.207,   1.362,   1.597,    1.959,    2.463,
  "Mg II",  2.001,   2.004,   2.010,   2.020,   2.036,    2.058,    2.087,
  "K I",    2.488,   3.000,   3.858,   5.078,   6.743,    8.938,    11.503,
  "K II",   1,       1,       1,       1,       1,        1,        1,
  "Na I",   2.166,   2.401,   2.868,   3.635,   4.779,    6.548,    8.724,
  "Na II",  1,       1,       1,       1,       1,        1,        1
)
# log-linear interpolation of the table, in the form cf_libs() accepts
partition <- function(species, temperature) {
  grid <- as.numeric(names(partition_table)[-1])
  vapply(species, function(s) {
    u <- unlist(partition_table[partition_table$species == s, -1])
    exp(stats::approx(grid, log(u), xout = temperature)$y)
  }, numeric(1))
}

## ----cf-mean------------------------------------------------------------------
cf_mean <- line_intensities(mean_spectrum, cf_atomic) |>
  cf_libs(electron_density = ne, partition = partition, ionization_energy = cf_ionization)
cf_mean

## ----cf-plot, fig.height = 4.5------------------------------------------------
plot_boltzmann(cf_mean)

## ----cf-compare---------------------------------------------------------------
laboratory <- forageLIBS |>
  summarise(across(c(Ca, Mg, K, Na), \(v) mean(v, na.rm = TRUE))) |>
  unlist()
cf_mean$composition |>
  select(element, cf_libs = mass_fraction) |>
  mutate(laboratory = laboratory[element] / sum(laboratory),
         ratio = cf_libs / laboratory)

## ----cf-numbers, include = FALSE----------------------------------------------
ratios <- cf_mean$composition |>
  mutate(laboratory = laboratory[element] / sum(laboratory), ratio = mass_fraction / laboratory)
ratio_of <- \(e) round(ratios$ratio[ratios$element == e], 1)

