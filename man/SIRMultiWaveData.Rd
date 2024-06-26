% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{SIRMultiWaveData}
\alias{SIRMultiWaveData}
\title{Sample Multi-wave Results}
\format{
A list of data.frame objects:
\describe{
  \item{Wave 0}{The initial points used in other examples}
  \item{Wave 1}{Points generated from the wave 1 emulators}
  \item{Wave 2}{Points generated from the wave 2 emulators}
  \item{Wave 3}{Points generated from the wave 3 emulators}
}
}
\usage{
SIRMultiWaveData
}
\description{
An rda object containing four data.frames: an initial set of points
also provided in \code{\link{SIRSample}}, and
the 90 points generated at each of three subsequent waves. The trained
emulators are provided in \code{\link{SIRMultiWaveEmulators}}.
}
\keyword{datasets}
