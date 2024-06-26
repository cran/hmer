% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlations.R
\name{matern}
\alias{matern}
\title{Matern correlation function}
\usage{
matern(x, xp, hp)
}
\arguments{
\item{x}{A data.frame of rows corresponding to position vectors}

\item{xp}{A data.frame of rows corresponding to position vectors}

\item{hp}{The hyperparameters nu (smoothness) and theta (correlation length), as a named list}
}
\value{
The Matern correlation between x and xp.
}
\description{
For points \code{x}, \code{xp}, and a pair of hyperparameters \code{nu} and \code{theta}, gives
the Matern correlation between the two points.
}
\details{
At present, only half-integer arguments for nu are supported.
}
\examples{
matern(data.frame(a=1), data.frame(a=2), list(nu = 1.5, theta = 0.1))
#> 5.504735e-07
matern(data.frame(a=1,b=2,c=-1),data.frame(a=1.5,b=2.9,c=-0.7), list(nu = 1.5, theta = 0.2))
#> 0.0009527116
}
\references{
Rasmussen & Williams (2005) <ISBN: 9780262182539>
}
