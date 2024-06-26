% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiwaveplot.R
\name{diagnostic_wrap}
\alias{diagnostic_wrap}
\title{Diagnostic plots for wave outputs}
\usage{
diagnostic_wrap(
  waves,
  targets,
  output_names = names(targets),
  input_names = names(waves[[1]])[!names(waves[[1]]) \%in\% names(targets)],
  directory = NULL,
  s.heights = rep(1000, 4),
  s.widths = s.heights,
  include.norm = TRUE,
  include.log = TRUE,
  ...
)
}
\arguments{
\item{waves}{The wave points, as a list of data.frames.}

\item{targets}{The output targets.}

\item{output_names}{The outputs to plot.}

\item{input_names}{The inputs to plot.}

\item{directory}{The location of files to be saved (if required).}

\item{s.heights}{The heights of the saved pngs (if directory is not NULL).}

\item{s.widths}{The widths of the saved pngs (if directory is not NULL).}

\item{include.norm}{Should normalized versions of simulator_plot and wave_dependencies be made?}

\item{include.log}{Should the log-scale version of simulator_plot be made?}

\item{...}{Optional parameters (eg \code{p_size}, \code{l_wid}, ...)}
}
\value{
The set of plots (either into console or saved).
}
\description{
A wrapper function for the set of diagnostic plots for multiple waves.
}
\details{
The functions \code{\link{simulator_plot}}, \code{\link{wave_points}}, \code{\link{wave_points}},
and \code{\link{wave_dependencies}} are called, one after the other, to allow diagnosis of waves
of emulation.

The \code{directory} option should be used as follows. If the desired location is in fact
a folder, it should end in "/"; if instead the structure requires each plot to be saved with a
prefix, then it should be provided. For example, \code{directory = "Plots/"} in the first event
or \code{directory = "Plots/unique-identifier"} in the second event.
}
\examples{
\donttest{
 diagnostic_wrap(SIRMultiWaveData, SIREmulators$targets)
 diagnostic_wrap(SIRMultiWaveData, SIREmulators$targets,
  input_names = c('aSI', 'aIR'), output_names = c('nI', 'nR'),
  p_size = 0.8, l_wid = 0.8, wave_numbers = 1:3, zero_in = FALSE, surround = TRUE)
  }
}
\seealso{
Other visualisation tools: 
\code{\link{behaviour_plot}()},
\code{\link{effect_strength}()},
\code{\link{emulator_plot}()},
\code{\link{hit_by_wave}()},
\code{\link{output_plot}()},
\code{\link{plot_actives}()},
\code{\link{plot_lattice}()},
\code{\link{plot_wrap}()},
\code{\link{simulator_plot}()},
\code{\link{space_removed}()},
\code{\link{validation_pairs}()},
\code{\link{wave_dependencies}()},
\code{\link{wave_points}()},
\code{\link{wave_values}()}
}
\concept{visualisation tools}
