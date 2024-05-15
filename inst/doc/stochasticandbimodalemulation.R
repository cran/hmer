## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "stoch-"
)

## ----setup, include = FALSE---------------------------------------------------
library(hmer)
library(ggplot2)
set.seed(123)

gillespied=function (N, T=400, dt=1, ...)
{
  tt=0
  n=T%/%dt
  x=N$M
  S=t(N$Post-N$Pre)
  u=nrow(S)
  v=ncol(S)
  xmat=matrix(0,ncol=u,nrow=n)
  i=1
  target=0
  repeat {
    h=N$h(x, tt, ...)
    h0=sum(h)
    if (h0<1e-10)
      tt=1e99
    else
      tt=tt+rexp(1,h0)
    while (tt>=target) {
      xmat[i,]=x
      i=i+1
      target=target+dt
      if (i>n)
        #			        return(ts(xmat,start=0,deltat=dt))
        return(xmat)
    }
    j=sample(v,1,prob=h)
    x=x+S[,j]
  }
}

Num <- 1000
N=list()
N$M=c(995,5,0)
N$Pre = matrix(c(1,0,0,
                 0,1,0,
                 0,0,1), ncol = 3, byrow = TRUE)
N$Post = matrix(c(0,1,0,
                  0,0,1,
                  1,0,0), ncol = 3, byrow = TRUE)

N$h=function(x,t,th=rep(1,3))
{
  Num = x[1]+x[2]+x[3]
  return(
    c(th[1]*x[1]*x[2]/Num,
      th[2]*x[2],
      th[3]*x[3])
  )
}

get_results <- function(params, nreps = 100, outs = c("I", "R"), times = seq(0, 60, by = 5), raw = FALSE) {
  tseq <- 0:max(times)
  arra <- array(0, dim = c(max(tseq)+1, 3, nreps))
  for(i in 1:nreps) arra[,,i] <- gillespied(N, T=max(times) + 1 + 0.001, dt=1, th=params)
  if (raw) return(arra)
  collected <- list()
  for (i in 1:nreps) {
    relev <- c(arra[times+1, which(c("S", "I", "R") %in% outs), i])
    names <- unlist(purrr::map(outs, ~paste0(., times, sep = "")))
    relev <- setNames(relev, names)
    collected[[i]] <- relev
  }
  input_dat <- setNames(data.frame(matrix(rep(params, nreps), ncol = length(params), byrow = TRUE)), names(params))
  return(cbind(input_dat, do.call('rbind', collected)))
}

## ----birth-death-plot, fig.height = 7, fig.width = 7--------------------------
plot(unique(BirthDeath$training[,1:2]), pch = 16, col = c(rep('red', 5), rep('black', 15)))
legend('top', legend = c("High rep", "Low rep"), col = c('red', 'black'), pch = 16, inset = c(0.1, 0.1))

## ----define-ranges------------------------------------------------------------
ranges <- list(
  lambda = c(0, 0.08),
  mu = c(0.04, 0.13)
)
targets = list(
  Y = c(90, 110)
)

## ----make-variance-em---------------------------------------------------------
stochastic_em <- emulator_from_data(BirthDeath$training, names(targets), ranges, emulator_type = "variance")
stochastic_em

## ----print-class--------------------------------------------------------------
class(stochastic_em$expectation$Y)

## ----plot-stoch-ems, fig.width = 7, fig.height = 7----------------------------
emulator_plot(stochastic_em$expectation$Y)
emulator_plot(stochastic_em$variance$Y)

## ----plot-stoch-var, fig.width = 7, fig.height = 7----------------------------
emulator_plot(stochastic_em$expectation$Y, 'sd') + geom_point(data = BirthDeath$training, aes(x = lambda, y = mu))

## ----plot-stoch-imp, fig.width = 7, fig.height = 7----------------------------
emulator_plot(stochastic_em$expectation$Y, 'imp', targets = targets)

## ----validate-stoch-em, fig.height = 7/3, fig.width = 7-----------------------
# Mean emulator: familiar code
mean_em_validation <- validation_diagnostics(stochastic_em$expectation, targets, BirthDeath$validation, plt = TRUE, row = 1)
# Variance emulator (no targets; see below)
var_em_validation <- suppressWarnings(validation_diagnostics(stochastic_em$variance, validation = BirthDeath$validation, plt = TRUE, row = 1))
# No targets: LOO cross-validation performed
no_validation <- validation_diagnostics(stochastic_em$expectation, targets, plt = TRUE, row = 1)

## ----validate-var-em, fig.width = 7, fig.height = 7/3-------------------------
fake_var_targ <- list(Y = c(30, 50))
var_em_fake_target <- validation_diagnostics(stochastic_em$variance, fake_var_targ, BirthDeath$validation, plt = TRUE, row = 1)

## ----propose-var-em, fig.height = 7, fig.width = 7----------------------------
new_points <- generate_new_design(stochastic_em, 30, targets, resample = 0)
plot_wrap(new_points, ranges)

## ----multistate-demo, fig.height = 7, fig.width = 7---------------------------
example_params <- c(0.45, 0.25, 0.025)
model_result <- get_results(example_params, raw = TRUE)
plot(0:60, ylim=c(0,600), ty="n")
for(j in 2:3) for(i in 1:100) lines(0:60, model_result[,j,i], col=(2:4)[j], lwd=0.3, xlab = "Time", ylab = "Number")
legend('topleft', legend = c('Recovered', "Infected"), lty = 1, col = c(4, 3), inset = c(0.05, 0.05))

## ----define-multistate--------------------------------------------------------
SIR_ranges <- list(
  aSI = c(0.1, 0.8),
  aIR = c(0, 0.5),
  aSR = c(0, 0.05)
)
training <- SIR_stochastic$training
validation <- SIR_stochastic$validation
bim_targets <- list(
  I10 = list(val = 35, sigma = 3.5),
  I25 = list(val = 147, sigma = 14.7),
  I50 = list(val = 55, sigma = 5.5),
  R10 = list(val = 29, sigma = 2.9),
  R25 = list(val = 276, sigma = 27.6),
  R50 = list(val = 579, sigma = 57.9)
)

## ----train-multistate---------------------------------------------------------
bim_ems <- emulator_from_data(training, names(bim_targets), SIR_ranges, emulator_type = "multistate")
names(bim_ems)

## ----plot-multistate, fig.width = 7, fig.height = 7---------------------------
emulator_plot(bim_ems$mode1$expectation$R50)
emulator_plot(bim_ems$mode2$expectation$R50)

## ----validate-multistate, fig.width = 7, fig.height = 7-----------------------
bim_invalid <- validation_diagnostics(bim_ems, bim_targets, validation)

## ----plot-modes, fig.width = 7, fig.height = 7--------------------------------
emulator_plot(bim_ems$mode1)
emulator_plot(bim_ems$mode2)

## ----multistate-imp, fig.width = 7, fig.height = 7----------------------------
emulator_plot(bim_ems, 'nimp', targets = bim_targets)

## ----multistate-indiv-imp, fig.width = 7, fig.height = 7----------------------
emulator_plot(bim_ems$mode1, 'nimp', targets = bim_targets)
emulator_plot(bim_ems$mode2, 'nimp', targets = bim_targets)

## ----propose-multistate, fig.width = 7, fig.height = 7------------------------
new_points <- generate_new_design(bim_ems, 100, bim_targets, resample = 0)
plot_wrap(new_points, SIR_ranges)

