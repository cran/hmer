## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "demo-"
)

## ----setup--------------------------------------------------------------------
library(hmer)
library(ggplot2)
set.seed(1)

## ----train-ems-1--------------------------------------------------------------
input_ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
output_names <- c('nS', 'nI', 'nR')
wave1_emulators <- emulator_from_data(SIRSample$training, output_names, input_ranges)

## ----print-ems----------------------------------------------------------------
wave1_emulators$nS
wave1_emulators$nS$corr

## ----custom-ems---------------------------------------------------------------
custom_emulators <- emulator_from_data(SIRSample$training, output_names, input_ranges,
specified_priors = list(hyper_p = c(0.6, 0.7, 0.85)))
custom_emulators$nR$corr

## ----matern-ems---------------------------------------------------------------
matern_emulators <- emulator_from_data(SIRSample$training, output_names, input_ranges, corr_name = "matern")

## ----set-targets--------------------------------------------------------------
targets <- list(
  nS = c(580, 651),
  nI = list(val = 169, sigma = 8.45),
  nR = c(199, 221)
)

## ----validate, fig.width = 7, fig.height = 7----------------------------------
invalid_points <- validation_diagnostics(wave1_emulators, targets, SIRSample$validation, plt = TRUE)

## ----em-plot, fig.width = 7, fig.height = 7-----------------------------------
emulator_plot(wave1_emulators)
emulator_plot(wave1_emulators$nS, params = c('aSI', 'aIR'), fixed_vals = list(aSR = 0.03))
emulator_plot(wave1_emulators$nI, plot_type = 'var')
emulator_plot(wave1_emulators$nI, plot_type = 'sd')

## ----em-plot-augment, fig.width = 7, fig.height = 7---------------------------
emulator_plot(wave1_emulators$nS, params = c('aIR', 'aSI')) + ggplot2::geom_point(data = SIRSample$training, ggplot2::aes(x = aSI, y = aIR))

## ----em-plot-imp, fig.width = 7, fig.height = 7-------------------------------
emulator_plot(wave1_emulators, plot_type = 'imp', targets = targets)
emulator_plot(wave1_emulators, plot_type = 'nimp', targets = targets, cb = TRUE)

## ----lattice-plot, fig.width = 7, fig.height = 7------------------------------
plot_lattice(wave1_emulators, targets)

## ----space-removed, fig.width = 7, fig.height = 7-----------------------------
space_removed(wave1_emulators, targets, ppd = 15) + geom_vline(xintercept = 3, lty = 2)

## ----gen-new-runs, fig.height = 7, fig.width = 7------------------------------
new_points <- generate_new_design(wave1_emulators, 90, targets)
plot(rbind(rbind(SIRSample$training, SIRSample$validation)[,names(input_ranges)], new_points), pch = 16, cex = 0.8, col = rep(c('black', 'blue'), each = 90))

## ----create-waves-------------------------------------------------------------
wave.points <- list(SIRSample$training[,c('aSI', 'aIR', 'aSR')], new_points)

## ----plot-waves, fig.height = 7, fig.width = 7--------------------------------
wave_points(wave.points, c('aSI', 'aIR', 'aSR'))

## ----multi-wave-plots, fig.height = 7, fig.width = 7--------------------------
wave_values(SIRMultiWaveData, targets)
wave_dependencies(SIRMultiWaveData, targets)
simulator_plot(SIRMultiWaveData, targets, barcol = 'white')

## ----full-wave----------------------------------------------------------------
f_w <- suppressWarnings(full_wave(do.call('rbind.data.frame', SIRSample), list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)), targets))

