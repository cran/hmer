## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "lowdim-"
)

## ----setup--------------------------------------------------------------------
library(hmer)
library(lhs)
library(ggplot2)
set.seed(1)

## ----1d-function--------------------------------------------------------------
func <- function(x) {
  2*x + 3*x*sin(5*pi*(x-0.1)/0.4)
}

## ----1d-data------------------------------------------------------------------
data1d <- data.frame(x = seq(0.05, 0.5, by = 0.05), f = func(seq(0.05, 0.5, by = 0.05)))

## ----train-1d-emulator--------------------------------------------------------
ranges1d <- list(x = c(0, 0.6))
em1d <- emulator_from_data(data1d, c('f'), ranges1d)
em1d

## ----print-1d-emulator--------------------------------------------------------
emulator_from_data(data1d, c('f'), ranges1d, adjusted = FALSE)$f
em1d$f$o_em

## ----get-1d-exp-and-cov-------------------------------------------------------
em1d$f$get_exp(data1d) - data1d$f
em1d$f$get_cov(data1d)

## ----predict-1d-vals----------------------------------------------------------
test1d <- data.frame(x = seq(0, 0.6, by = 0.001))
em1d_exp <- em1d$f$get_exp(test1d)
em1d_var <- em1d$f$get_cov(test1d)

## ----plot-1d-first, fig.width = 7, fig.height = 7-----------------------------
#> Define a data.frame for the plotting
plot1d <- data.frame(
  x = test1d$x,
  f = func(test1d$x),
  E = em1d_exp,
  max = em1d_exp + 3*sqrt(abs(em1d_var)),
  min = em1d_exp - 3*sqrt(abs(em1d_var))
)
plot(data = plot1d, f ~ x, ylim = c(min(plot1d[,-1]), max(plot1d[-1])),
     type = 'l', main = "Emulation of a Simple 1d Function", xlab = "x", ylab = "f(x)")
lines(data = plot1d, E ~ x, col = 'blue')
lines(data = plot1d, max ~ x, col = 'red', lty = 2)
lines(data = plot1d, min ~ x, col = 'red', lty = 2)
points(data = data1d, f ~ x, pch = 16, cex = 1)
legend('bottomleft', inset = c(0.05, 0.05), legend = c("Function value", "Emulated value", "Uncertainty Bounds"),
       col = c('black', 'blue', 'red'), lty = c(1,1,2))

## ----define-1d-target---------------------------------------------------------
target1d <- list(f = list(val = 0, sigma = 0.05))

## ----define-1d-alt-target-----------------------------------------------------
target1d_alt <- list(f = c(-0.1, 0.1))

## ----gen-1d-points-first------------------------------------------------------
new_points1d <- generate_new_design(em1d, 10, target1d, method = 'lhs')

## ----plot-1d-results-first, fig.width = 7, fig.height = 7---------------------
col_func <- function(x, max_val) {
  if (x <= 3) return(rgb(221*x/3, 255, 0, maxColorValue = 256))
  return(rgb(221+34*(x-3)/(max_val-3), 255*(1-(x-3)/(max_val-3)), 0, maxColorValue = 256))
}
imps <- em1d$f$implausibility(plot1d$x, target1d$f)
col.pal <- purrr::map_chr(imps, col_func, max(imps))

plot(data = plot1d, f ~ x, ylim = c(min(plot1d[,-1]-1), max(plot1d[,-1])),
     type = 'l', main = "Emulation of a Simple 1d Function", xlab = "x", ylab = "f(x)")
lines(data = plot1d, E ~ x, col = 'blue')
lines(data = plot1d, max ~ x, col = 'red', lty = 2)
lines(data = plot1d, min ~ x, col = 'red', lty = 2)
points(data = data1d, f ~ x, pch = 16, cex = 1)
abline(h = target1d$f$val, lty = 2)
abline(h = target1d$f$val + 3*target1d$f$sigma, lty = 2)
abline(h = target1d$f$val - 3*target1d$f$sigma, lty = 2)
points(unlist(new_points1d, use.names = F), y = func(unlist(new_points1d, use.names = F)), pch = 16, col = 'blue')
points(plot1d$x, rep(min(plot1d[,-1])-0.5, length(plot1d$x)), col = col.pal, pch = 16)
legend('bottomleft', inset = c(0.05, 0.1), legend = c("Function value", "Emulated value", "Uncertainty Bounds"),
       col = c('black', 'blue', 'red'), lty = c(1,1,2))

## ----wave-2-1d, fig.width = 7, fig.height = 7---------------------------------
new_data1d <- data.frame(x = unlist(new_points1d, use.names = F), f = func(unlist(new_points1d, use.names = F)))

em1d_2 <- emulator_from_data(new_data1d, c('f'), ranges1d)

plot1d_2 <- data.frame(
  x = plot1d$x, f = plot1d$f,
  E = em1d_2$f$get_exp(test1d),
  min = em1d_2$f$get_exp(test1d) - 3*sqrt(abs(em1d_2$f$get_cov(test1d))),
  max = em1d_2$f$get_exp(test1d) + 3*sqrt(abs(em1d_2$f$get_cov(test1d)))
)
plot(data = plot1d_2, f ~ x, ylim = c(min(plot1d_2[,-1]), max(plot1d_2[,-1])),
     type = 'l', main = "Emulator of a Simple 1-dimensional Function: Wave 2", xlab = "Parameter value", ylab = "Function value")
lines(data = plot1d_2, E ~ x, col = 'blue')
lines(data = plot1d_2, max ~ x, col = 'red', lty = 2)
lines(data = plot1d_2, min ~ x, col = 'red', lty = 2)
points(data = new_data1d, f ~ x, pch = 16, cex = 1)
legend('topleft', inset = c(0.05, 0.05), legend = c("Function value", "Emulated value", "Uncertainty Bounds"), col = c('black', 'blue', 'red'), lty = c(1,1,2))

## ----gen-1d-points-2, fig.width = 7, fig.height = 7---------------------------
new_new_points1d <- generate_new_design(c(em1d_2, em1d), 10, z = target1d, method = 'lhs')

plot(data = plot1d_2, f ~ x, ylim = c(min(plot1d_2[,-1]), max(plot1d_2[,-1])),
     type = 'l', main = "Emulator of a Simple 1-dimensional Function: Wave 2", xlab = "Parameter value", ylab = "Function value")
lines(data = plot1d_2, E ~ x, col = 'blue')
lines(data = plot1d_2, max ~ x, col = 'red', lty = 2)
lines(data = plot1d_2, min ~ x, col = 'red', lty = 2)
points(data = new_data1d, f ~ x, pch = 16, cex = 1)
legend('topleft', inset = c(0.05, 0.05), legend = c("Function value", "Emulated value (wave 2)", "Uncertainty Bounds"), col = c('black', 'blue', 'red'), lty = c(1,1,2))
abline(h = target1d$f$val, lty = 2)
abline(h = target1d$f$val + 3*target1d$f$sigma, lty = 2)
abline(h = target1d$f$val - 3*target1d$f$sigma, lty = 2)
points(x = unlist(new_new_points1d, use.names = F), y = func(unlist(new_new_points1d, use.names = F)), pch = 16, col = 'blue')

## ----2d-introduction----------------------------------------------------------
func1 <- function(x) {
  2*cos(1.2*x[[1]]-2) + 3*sin(-0.8*x[[2]]+1)
}
func2 <- function(x) {
  x[[2]]*sin(x[[1]]) - 3*cos(x[[1]]*x[[2]])
}
initial_data <- setNames(data.frame(sweep(maximinLHS(20, 2) - 1/2, 2, c(pi, 2), "*")), c('x', 'y'))
validation_data <- setNames(data.frame(sweep(maximinLHS(20, 2) - 1/2, 2, c(pi, 2), "*")), c('x', 'y'))
initial_data$f1 <- apply(initial_data, 1, func1)
initial_data$f2 <- apply(initial_data, 1, func2)
validation_data$f1 <- apply(validation_data, 1, func1)
validation_data$f2 <- apply(validation_data, 1, func2)

## ----set-2d-ranges-and-targets------------------------------------------------
ranges2d <- list(x = c(-pi/2, pi/2), y = c(-1, 1))
targets2d <- list(
  f1 = list(val = func1(c(0.1, 0.4)), sigma = sqrt(0.1)),
  f2 = list(val = func2(c(0.1, 0.4)), sigma = sqrt(0.005))
)

## ----train-2d-ems-------------------------------------------------------------
ems2d <- emulator_from_data(initial_data, c('f1','f2'), ranges2d)
ems2d

## ----validate-2d-ems, fig.width = 7, fig.height = 5---------------------------
invalid_points <- validation_diagnostics(ems2d, targets2d, validation_data, plt = TRUE, row = 2)

## ----changed-2d-validation, fig.width = 7, fig.height = 5---------------------
modified_em <- ems2d$f1$mult_sigma(2)$set_hyperparams(list(theta = 0.5))
new_invalid <- validation_diagnostics(list(f1 = modified_em, f2 = ems2d$f2), targets2d, validation_data, plt = TRUE, row = 2)

## ----plot-2d-ems, fig.width = 7, fig.height = 4-------------------------------
emulator_plot(ems2d, include_legend = FALSE)
emulator_plot(ems2d, plot_type = 'sd', include_legend = FALSE)

## ----plot-2d-ems-with-augment, fig.width = 7, fig.height = 6.5----------------
emulator_plot(ems2d$f1, plot_type = 'sd') + geom_point(data = initial_data, aes(x = x, y = y))

## ----plot-2d-ems-corr, fig.width = 7, fig.height = 7--------------------------
inflated_corr <- ems2d$f1$set_hyperparams(list(theta = 1))
emulator_plot(inflated_corr, 'sd') + geom_point(data = initial_data, aes(x = x, y = y))

## ----plot-2d-ems-imp, fig.width = 7, fig.height = 5---------------------------
emulator_plot(ems2d, plot_type = 'imp', targets = targets2d)
emulator_plot(ems2d, plot_type = 'nimp', targets = targets2d)

## ----gen-2d-points-first, fig.width = 7, fig.height = 7-----------------------
new_points2d <- generate_new_design(ems2d, 40, targets2d, resample = 0)
new_data2d <- data.frame(x = new_points2d$x, y = new_points2d$y,
                         f1 = apply(new_points2d, 1, func1), f2 = apply(new_points2d, 1, func2))
plot_wrap(new_data2d[,1:2], ranges2d)

## ----next-2d-points-----------------------------------------------------------
sample2d <- sample(40, 20)
train2d <- new_data2d[sample2d,]
valid2d <- new_data2d[-sample2d,]

## ----next-2d-ems--------------------------------------------------------------
new_ranges2d <- list(x = c(-pi/2, 0.6), y = c(-1, 1))
ems2d_2 <- emulator_from_data(train2d, c('f1', 'f2'), new_ranges2d)
ems2d_2

## ----validate-2d-ems-2, fig.width = 7, fig.height = 5-------------------------
invalid_points2 <- validation_diagnostics(ems2d_2, targets2d, valid2d, plt = TRUE, row = 2)

## ----modify-2d----------------------------------------------------------------
ems2d_2$f1 <- ems2d_2$f1$mult_sigma(1.4)$set_hyperparams(list(theta = 1/3))

## ----plot-2d-imp-2, fig.width = 7, fig.height = 5-----------------------------
emulator_plot(ems2d_2, 'imp', targets = targets2d)

## ----gen-2d-points-2, fig.width = 7, fig.height = 7---------------------------
new_new_points2d <- generate_new_design(c(ems2d_2, ems2d), 40, targets2d, resample = 0)
plot_wrap(new_new_points2d, ranges2d)

## ----final-2d-points----------------------------------------------------------
new_new_data2d <- data.frame(x = new_new_points2d$x, y = new_new_points2d$y,
                             f1 = apply(new_new_points2d, 1, func1), f2 = apply(new_new_points2d, 1, func2))

## ----combine-2d-waves---------------------------------------------------------
all_waves <- list(rbind(initial_data, validation_data), new_data2d, new_new_data2d)

## ----plot-2d-waves, fig.width = 7, fig.height = 7-----------------------------
simulator_plot(all_waves, z = targets2d)
wave_points(all_waves, names(ranges2d))
wave_values(all_waves, targets2d, l_wid = 0.8)

## ----uncertainty-compare------------------------------------------------------
c(ems2d_2$f1$u_sigma, ems2d_2$f2$u_sigma)
c(targets2d$f1$sigma, targets2d$f2$sigma)

