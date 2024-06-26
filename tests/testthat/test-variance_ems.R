v_em <- emulator_from_data(BirthDeath$training,
                                    c('Y'),
                                    list(lambda = c(0, 0.08), mu = c(0.04, 0.13)),
                                    verbose = FALSE, beta.var = TRUE,
                           emulator_type = "variance")

test_that("Variance Emulators", {
  expect_equal(
    class(v_em$expectation$Y),
    c("Hierarchical", "Emulator", "R6")
  )
  expect_equal(
    v_em$expectation$Y$em_type,
    "mean"
  )
  expect_equal(
    v_em$variance$Y$em_type,
    "variance"
  )
})

test_that("Batch runs", {
  many_points <- data.frame(lambda = runif(1400, 0, 0.08),
                            mu = runif(1400, 0.04, 0.13))
  expect_equal(
    length(c(v_em$variance$Y$get_exp(many_points))),
    1400
  )
  expect_equal(
    length(c(v_em$expectation$Y$get_cov(many_points))),
    1400
  )
  expect_equal(
    length(c(v_em$expectation$Y$implausibility(many_points,
                                               list(Y = c(90, 110))$Y))),
    1400
  )
})

em <- v_em$variance$Y
test_train <- unique(BirthDeath$training[,c('lambda', 'mu')])[1:5,]
test_points <- unique(BirthDeath$validation[,c('lambda', 'mu')])[1:5,]
test_that("Modifying priors and functional sigma", {
  em_2 <- em$set_sigma(2)
  expect_equal(
    em_2$u_sigma,
    2
  )
  em_3 <- em_2$mult_sigma(2)
  expect_equal(
    em_3$u_sigma,
    4
  )
  em_4 <- em_2$set_hyperparams(
    hp = list(theta = 0.75),
    nugget = 0.1
  )
  expect_equal(
    em_4$corr$hyper_p$theta,
    0.75
  )
  expect_equal(
    em_4$corr$nugget,
    0.1
  )
  em_sigma <- em$set_sigma(function(x) x[[1]]*5)
  expect_false(
    all(em_sigma$get_cov(test_train) == 0)
  )
  expect_equal(
    dim(em_sigma$get_cov(test_train[1:3,],
                         test_train[2:5,],
                         full = TRUE, check_neg = FALSE)),
    c(3, 4)
  )
  em_sigma_2 <- em_sigma$mult_sigma(2)
  expect_equal(
    em_sigma_2$u_sigma(c(0.01, 0)),
    0.1
  )
})

test_that("Modifying priors and functional sigma - untrained", {
  em_o <- em$o_em
  em_o2 <- em_o$set_sigma(2)
  expect_equal(
    em_o2$u_sigma,
    2
  )
  em_o3 <- em_o2$mult_sigma(2)
  expect_equal(
    unname(em_o3$get_cov(test_points[1,,drop=FALSE])),
    359.2923,
    tolerance = 1e-4
  )
  expect_equal(
    em_o3$u_sigma,
    2
  )
  em_o4 <- em_o2$set_hyperparams(
    hp = list(theta = 0.7),
    nugget = 0.3
  )
  expect_equal(
    em_o4$corr$hyper_p$theta,
    0.7
  )
  expect_equal(
    em_o4$corr$nugget,
    0.3
  )
})

test_that("Printing works", {
  expect_output(
    print(em),
    "Parameters and ranges"
  )
  expect_output(
    print(em),
    "Regression surface Variance"
  )
  expect_output(
    print(em),
    "Bayes-adjusted emulator - prior specifications listed"
  )
})

### Covariance Emulation
test_that("Covariance emulation building - basic", {
  cov_ems <- emulator_from_data(
    SIR_stochastic$training, c("I10", "I25", "R10", "R25"),
    list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)),
    emulator_type = "covariance", verbose = FALSE
  )
  expect_equal(
    dim(cov_ems$variance$get_matrix()),
    c(4,4)
  )
  expect_equal(
    class(cov_ems$variance),
    c("EmulatorMatrix", "R6")
  )
})

test_that("Covariance emulation building - specified covariance elements", {
  cov_ems_spec <- emulator_from_data(
    SIR_stochastic$training, c("I10", "I25", "R10", "R25"),
    list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05)),
    verbose = FALSE,
    emulator_type = "covariance", covariance_opts = list(
      matrix = matrix(c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE,
                        TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE), nrow = 4))
  )
  expect_equal(
    class(cov_ems_spec$variance$get_matrix()[[2,3]]),
    c("EmProto", "Emulator", "R6")
  )
  cov_preds <- cov_ems_spec$variance$get_exp(unique(SIR_stochastic$training[,1:3])[1:3,])
  expect_equal(
    dim(cov_preds),
    c(4,4,3)
  )
  expect_true(
    all(apply(cov_preds, 3, function(x) all(round(eigen(x)$values, 6) >= 0)))
  )
  cov_covs <- cov_ems_spec$variance$get_cov(unique(SIR_stochastic$training[,1:3])[1:3,])
  expect_equal(
    dim(cov_covs),
    c(4,4,3)
  )
  cov_uncert <- cov_ems_spec$variance$get_uncertainty(unique(SIR_stochastic$training[,1:3])[1:3,],
                                              cov_ems_spec$expectation)
  expect_equal(
    dim(cov_uncert),
    c(4,4,3)
  )
})

test_that("Variance emulation - point proposal", {
  pts <- generate_new_design(v_em, 100, list(Y = c(90, 105)), verbose = FALSE)
  expect_equal(
    nrow(pts),
    100
  )
})

test_that("Variance emulation - point proposal with seek_good", {
  skip_on_cran()
  pts <- generate_new_design(v_em, 100, list(Y = c(90, 105)), verbose = FALSE,
                             opts = list(seek = 10))
  expect_equal(
    nrow(pts),
    100
  )
})
