# ----------------------
# Emulator Documentation
# ----------------------
#' @title Bayes Linear Emulator
#'
#' @description Creates a univariate emulator object.
#'
#'     The structure of the emulator is \code{f(x) = g(x) * beta + u(x)}, for
#'     regression functions \code{g(x)}, regression coefficients \code{beta},
#'     and correlation structure \code{u(x)}. An emulator can be created with
#'     or without data; the preferred method is to create an emulator based on
#'     prior specifications in the absence of data, then use that emulator
#'     with data to generate a new one (see examples).
#'
#' @name Emulator
#'
#' @importFrom R6 R6Class
#'
#' @section Constructor: \code{Emulator$new(basis_f, beta, u, ranges, ...)}
#'
#' @section Arguments:
#'
#'   Required:
#'
#'   \code{basis_f} A list of basis functions to be used. The constant function
#'   \code{function(x) 1} should be provided as the first element.
#'
#'   \code{beta} The specification for the regression parameters. This should
#'   be provided in the form \code{list(mu, sigma)}, where \code{mu} are the
#'   expectations of the coefficients (aligning with the ordering of \code{basis_f})
#'   and \code{sigma} the corresponding covariance matrix.
#'
#'   \code{u} The specifications for the correlation structure. This should
#'   be specified in the form \code{list(sigma, corr)}, where \code{sigma} is
#'   a single-valued object, and \code{corr} is a Correlator object.
#'
#'   \code{ranges} A named list of ranges for the input parameters, provided as
#'   a named list of length-two numeric vectors.
#'
#'   Optional:
#'
#'   \code{data} A \code{data.frame} consisting of the data with which to adjust
#'   the emulator, consisting of input values for each parameter and the output.
#'
#'   \code{out_name} The name of the output variable.
#'
#'   \code{a_vars} A logical vector indicating which variables are active for
#'   this emulator.
#'
#'   \code{discs} Model discrepancies: does not include observational error. Ideally
#'   split into \code{list(internal = ..., external = ...)}.
#'
#'   Internal:
#'
#'   \code{model} If a linear model, or otherwise, has been fitted to the data,
#'   it lives here.
#'
#'   \code{original_em} If the emulator has been adjusted, the unadjusted
#'   \code{Emulator} object is stored, for use of \code{set_sigma} or similar.
#'
#'   \code{multiplier} A multiplicative factor to be applied to u_sigma. Typically
#'   equal to 1, unless changes have been made by, for example, \code{mult_sigma}.
#'
#' @section Constructor Details:
#'
#'    The constructor must take, as a minimum: a list of vectorised basis
#'    functions, whose length is equal to the number of regression
#'    coefficients; a correlation structure, which can be non-stationary;
#'    and the parameter ranges, used to scale all inputs to the range [-1,1].
#'
#'    The construction of a correlation structure is detailed in the documentation
#'    for Correlator.
#'
#' @section Accessor Methods:
#'
#'    \code{get_exp(x, include_c)} Returns the emulator expectation at a point,
#'    or at a collection of points. If \code{include_c = FALSE}, the contribution
#'    made by the correlation structure is not included.
#'
#'    \code{get_cov(x, xp = NULL, full = FALSE, include_c)} Returns the covariance between
#'    collections of points \code{x} and \code{xp}. If \code{xp} is not supplied,
#'    then this is equivalent to \code{get_cov(x, x, ...)}; if \code{full = TRUE},
#'    then the full covariance matrix is calculated - this is FALSE by default
#'    due to most built-in uses requiring only the diagonal terms, and allows us
#'    to take advantage of computational tricks for efficiency.
#'
#'    \code{implausibility(x, z, cutoff = NULL)} Returns the implausibility for a
#'    collection of points \code{x}. The implausibility is the distance between the
#'    emulator expectation and a desired output value, weighted by the emulator
#'    variance and any external uncertainty. The target, z, should be specified
#'    as a named pair \code{list(val, sigma)}, or a single numeric value.
#'    If \code{cutoff = NULL}, the output is a numeric \code{I}; if \code{cutoff}
#'    is a numeric value, then the output is boolean corresponding to
#'    \code{I <= cutoff}.
#'
#'    \code{get_exp_d(x, p)} Returns the expectation of the derivative of the emulated
#'    function, E[f'(x)]. Similar in structure to \code{get_exp} but for the additional
#'    parameter \code{p}, which indicates which of the input dimensions the derivative
#'    is performed with respect to.
#'
#'    \code{get_cov_d(x, p1, xp = NULL, p2 = NULL, full = FALSE)} Returns the variance of
#'    the derivative of the emulated function, Var[f'(x)]. The arguments are similar to
#'    that of \code{get_cov}, but for the addition of parameters \code{p1} and \code{p2},
#'    which indicate the derivative directions. Formally, the output of this function is
#'    equivalent to Cov[df/dp1, df/dp2].
#'
#'    \code{print(...)} Returns a summary of the emulator specifications.
#'
#'    \code{plot(...)} A wrapper for \code{\link{emulator_plot}} for a single Emulator object.
#'
#' @section Object Methods:
#'
#'    \code{adjust(data, out_name)} Performs Bayes Linear Adjustment, given \code{data}.
#'    The data should contain all input parameters, even inactive ones, and the
#'    single output that we wish to emulate. \code{adjust} creates a new \code{Emulator}
#'    object with the adjusted expectation and variance resulting from Bayes
#'    Linear adjustment, allowing for the requisite predictions to be made using
#'    \code{get_exp} and \code{get_cov}.
#'
#'    \code{set_sigma(sigma)} Modifies the (usually constant) global variance of
#'    the correlation structure, \code{Var[u(X)]}. If the emulator has been trained,
#'    the original emulator is modified and Bayes Linear adjustment is again performed.
#'
#'    \code{mult_sigma(m)} Modifies the global variance of the correlation structure via
#'    a multiplicative factor. As with \code{set_sigma}, this change will chain through
#'    any prior emulators if the emulator in question is Bayes Linear adjusted.
#'
#'    \code{set_hyperparams(hp, nugget)} Modifies the underlying correlator for \code{u(x)}.
#'    Behaves in a similar way to \code{set_sigma} as regards trained emulators. See
#'    the Correlator documentation for details of \code{hp} and \code{nugget}.
#'
#'
#' @references
#'  Goldstein & Wooff (2007) <ISBN: 9780470065662>
#'
#'  Craig, Goldstein, Seheult & Smith (1998) <doi:10.1111/1467-9884.00115>
#'
#' @export
#'
#' @examples
#' basis_functions <- list(function(x) 1, function(x) x[[1]], function(x) x[[2]])
#' beta <- list(mu = c(1,2,3),
#'              sigma = matrix(c(0.5, -0.1, 0.2, -0.1, 1, 0, 0.2, 0, 1.5), nrow = 3))
#' u <- list(mu = function(x) 0, sigma = 3, corr = Correlator$new('exp_sq', list(theta = 0.1)))
#' ranges <- list(a = c(-0.5, 0.5), b = c(-1, 2))
#' em <- Emulator$new(basis_functions, beta, u, ranges)
#' em
#' # Individual evaluations of points
#' # Points should still be declared in a data.frame
#' em$get_exp(data.frame(a = 0.1, b = 0.1)) #> 0.6
#' em$get_cov(data.frame(a = 0.1, b = 0.1)) #> 9.5
#' # 4x4 grid of points
#' sample_points <- expand.grid(a = seq(-0.5, 0.5, length.out = 4), b = seq(-1, 2, length.out = 4))
#' em$get_exp(sample_points) # Returns 16 expectations
#' em$get_cov(sample_points) # Returns 16 variances
#' sample_points_2 <- expand.grid(a = seq(-0.5, 0.5, length.out = 3),
#'                                b = seq(-1, 2, length.out = 4))
#' em$get_cov(sample_points, xp = sample_points_2, full = TRUE) # Returns a 16x12 matrix of covariances
#'
#'
#' fake_data <- data.frame(a = runif(10, -0.5, 0.5), b = runif(10, -1, 2))
#' fake_data$c <- fake_data$a + 2*fake_data$b
#' newem <- em$adjust(fake_data, 'c')
#' all(round(newem$get_exp(fake_data[,names(ranges)]),5) == round(fake_data$c,5)) #>TRUE
#'
#' matern_em <- Emulator$new(basis_f = c(function(x) 1, function(x) x[[1]], function(x) x[[2]]),
#'  beta = list(mu = c(1, 0.5, 2), sigma = diag(0, nrow = 3)),
#'  u = list(corr = Correlator$new('matern', list(nu = 1.5, theta = 0.4))),
#'  ranges = list(x = c(-1, 1), y = c(0, 3)))
#' matern_em$get_exp(data.frame(x = 0.4, y = 2.3))
#'
#' newem_data <- Emulator$new(basis_functions, beta, u, ranges, data = fake_data)
#' all(round(newem$get_exp(fake_data[,names(ranges)]),5)
#'    == round(newem_data$get_exp(fake_data[,names(ranges)]), 5)) #>TRUE
#' newem$get_exp_d(sample_points, 'a')
#' newem$get_cov_d(sample_points, 'b', p2 = 'a')
NULL
