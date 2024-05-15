## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "handbook-"
)

## ----setup, include=FALSE-----------------------------------------------------
library(hmer)
library(lhs)
library(deSolve)
library(ggplot2)

ode_results <- function(parms, end_time = 25) {
  des = function(time, state, parms) {
    with(as.list(c(state, parms)), {
      dS <- aSR*R-aSI*I*S/(S+I+R)
      dI <- aSI*I*S/(S+I+R)-aIR*I
      dR <- aIR*I-aSR*R
      return(list(c(dS, dI, dR)))
    })
  }
  yini = c(S = 950, I = 50, R = 0)
  times = seq(0, end_time, by = 1)
  out = ode(yini, times, des, parms)
  return(out)
}

get_res <- function(inputs) {
  ode_out <- data.frame(t(apply(inputs, 1, function(x) ode_results(x)[11, c('S', 'I', 'R')])))
  return(setNames(cbind(inputs, ode_out), c(names(inputs), c('nS', 'nI', 'nR'))))
}
set.seed(12)

## ----bad-design---------------------------------------------------------------
bad_design <- data.frame(aSI = runif(30, 0.6, 0.8),
                         aIR = runif(30, 0, 0.1),
                         aSR = runif(30, 0, 0.05))

bad_results <- get_res(bad_design)
bad_emulators <- emulator_from_data(bad_results, c('nS', 'nI', 'nR'),
                                    ranges = list(aSI = c(0.6, 0.8), aIR = c(0, 0.1), aSR = c(0, 0.05)))

## ----bad-design-diagnostic, fig.height = 7, fig.width = 7---------------------
bad_design_v <- data.frame(aSI = runif(30, 0.6, 0.8),
                         aIR = runif(30, 0, 0.1),
                         aSR = runif(30, 0, 0.05))
bad_validation <- get_res(bad_design_v)
bad_targets <- list(
 nS = c(150, 700),
 nI = list(val = 169, sigma = 50),
 nR = c(160, 250)
)
invalid <- validation_diagnostics(bad_emulators, bad_targets, bad_validation)

## ----bad_design_generation----------------------------------------------------
try_generate <- generate_new_design(bad_emulators, 100, bad_targets, verbose = TRUE,
 opts = list(points.factor = 100))

## ----train-emulators-problem--------------------------------------------------
problem_ems <- emulator_from_data(problem_data$data, names(problem_data$targets)[7:9], problem_data$ranges, targets = problem_data$targets, more_verbose = FALSE)

## ----check-output, fig.height = 7, fig.width = 7------------------------------
plot(x = problem_data$data[,names(problem_data$ranges)[9]], y = problem_data$data[,names(problem_data$targets)[7]],
     pch = 16, xlim = problem_data$ranges[[9]], ylim = c(0, 1), xlab = "Parameter", ylab = "Output")
abline(h = problem_data$targets[[7]], lty = 2)

## ----behaviour-plot, fig.width = 7, fig.height = 7----------------------------
behaviour_plot(ems = list(problem_ems[[1]]), targets = problem_data$targets)

## ----check-output-extra, fig.height = 7, fig.width = 7------------------------
plot_data <- rbind(problem_data$data, problem_data$extra)
plot(x = plot_data[,names(problem_data$ranges)[9]], y = plot_data[,names(problem_data$targets)[7]],
     pch = 16, xlim = problem_data$ranges[[9]], ylim = c(0, 1), xlab = "Parameter", ylab = "Output",
     col = c(rep('black', nrow(problem_data$data)), rep('blue', nrow(problem_data$extra))))
abline(h = problem_data$targets[[7]], lty = 2)

## ----logit-example, fig.height = 7, fig.width = 7-----------------------------
logit_data_first <- problem_data$data
logit_data_second <- problem_data$extra
logit_data_first$output7 <- log(logit_data_first$output7/(1-logit_data_first$output7))
logit_data_second$output7 <- log(logit_data_second$output7/(1-logit_data_second$output7))
logit_data_first <- logit_data_first[is.finite(logit_data_first$output7),]
logit_data <- rbind(logit_data_first, logit_data_second)
plot(x = logit_data[,names(problem_data$ranges)[9]], y = logit_data[,names(problem_data$targets)[7]],
     pch = 16, xlim = problem_data$ranges[[9]], ylim = range(logit_data[,names(problem_data$targets)[7]]), xlab = "Parameter", ylab = "Output",
     col = c(rep('black', nrow(logit_data_first)), rep('blue', nrow(logit_data_second))))
abline(h = log(problem_data$targets[[7]]/(1-problem_data$targets[[7]])), lty = 2)

## ----problem-imp, fig.height = 7, fig.width = 7-------------------------------
plot(problem_ems[[1]], plot_type = 'imp', params = c('input9', 'input15'), fixed_vals = problem_data$extra[1,(1:21)[-c(9, 15)]], targets = problem_data$targets)

## ----transformed-imp, fig.height = 7, fig.width = 7---------------------------
new_emulator <- emulator_from_data(logit_data_first, c('output7'), problem_data$ranges)
plot(new_emulator$output7, plot_type = 'imp', params = c('input9', 'input15'), fixed_vals = logit_data_second[1, (1:21)[-c(9, 15)]], targets = list(output7 = log(problem_data$targets$output7/(1-problem_data$targets$output7))))

## ----plot-active, fig.height = 7, fig.width = 7-------------------------------
plot_actives(new_emulator)

## ----output-active-1, fig.height = 7, fig.width = 7---------------------------
plot(new_emulator$output7, params = c('input1', 'input8'),
              fixed_vals = problem_data$data[1,(1:21)[-c(1,8)]]) +
  geom_point(x = problem_data$data[1,1], y = problem_data$data[1,8])
plot(new_emulator$output7, plot_type = 'sd', params = c('input1', 'input8'),
              fixed_vals = problem_data$data[1,(1:21)[-c(1,8)]]) +
  geom_point(x = problem_data$data[1,1], y = problem_data$data[1,8])

## ----output-active-2, fig.height = 7, fig.width = 7---------------------------
plot(new_emulator$output7, params = c('input18', 'input20'),
              fixed_vals = problem_data$data[1,(1:21)[-c(18, 20)]]) +
  geom_point(x = problem_data$data[1,2], y = problem_data$data[1,3])
plot(new_emulator$output7, plot_type = 'sd', params = c('input18', 'input20'),
              fixed_vals = problem_data$data[1,(1:21)[-c(18,20)]]) +
  geom_point(x = problem_data$data[1,2], y = problem_data$data[1,3])

## ----dirac-delta, fig.height = 5, fig.width = 5-------------------------------
dm <- matrix(0, ncol = length(problem_data$ranges), nrow = 1000)
for (i in 1:length(problem_data$ranges)) dm[,i] <- rep(problem_data$data[1,i], 1000)
dm[,18] <- sort(c(problem_data$data[1,18], seq(problem_data$ranges[[18]][1], problem_data$ranges[[18]][2], length.out = 999)))
data1d <- setNames(data.frame(dm), names(problem_data$ranges))
vars <- new_emulator$output7$get_cov(data1d)
plot(x = data1d$input18, y = vars, type = 'l', xlab = "input18", ylab = 'Variance')

