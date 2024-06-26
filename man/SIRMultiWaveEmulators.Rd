% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SIRMultiWaveEmulators}
\alias{SIRMultiWaveEmulators}
\title{Sample Multi-wave Emulators}
\format{
A list containing \code{\link{Emulator}} objects:
\describe{
  \item{Wave 1}{Emulators trained on Wave 0, generating wave 1 points}
  \item{Wave 2}{Emulators trained on the results of the above wave 1 points}
  \item{Wave 3}{Emulators trained on the results of the above wave 2 points}
}
}
\usage{
SIRMultiWaveEmulators
}
\description{
An rda object containing three waves of emulators applied to
SIR model (described in \code{\link{SIRSample}}). The corresponding points
(both training and validation) are stored in \code{\link{SIRMultiWaveData}}.
}
\keyword{datasets}
