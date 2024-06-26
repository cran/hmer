% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiwaveplot.R
\name{wave_values}
\alias{wave_values}
\title{Multiple Wave Output Plotting}
\usage{
wave_values(
  waves,
  targets,
  output_names = names(targets),
  ems = NULL,
  surround = FALSE,
  restrict = FALSE,
  p_size = 1.5,
  l_wid = 1.5,
  zero_in = TRUE,
  wave_numbers = ifelse(zero_in, 0, 1):(length(waves) - ifelse(zero_in, 0, 1)),
  which_wave = ifelse(zero_in, 0, 1),
  upper_scale = 1,
  ...
)
}
\arguments{
\item{waves}{The list of data.frames, one for each set of outputs at that wave.}

\item{targets}{The output targets.}

\item{output_names}{The outputs to plot.}

\item{ems}{If provided, plots the emulator expectations and 3-standard deviations.}

\item{surround}{As in \code{\link{wave_points}}.}

\item{restrict}{Should the plotting automatically restrict to failing target windows?}

\item{p_size}{As in \code{\link{wave_points}}.}

\item{l_wid}{The width of the lines that create the target boxes.}

\item{zero_in}{Is a wave 0 included in the waves list?}

\item{wave_numbers}{Which waves to plot.}

\item{which_wave}{Scaling for lower plots (see description)}

\item{upper_scale}{Scaling for upper plots (ibid)}

\item{...}{Optional parameters (not to be used directly)}
}
\value{
A ggplot object.
}
\description{
Given multiple waves of points, produces pairs plots of the outputs.
}
\details{
This function operates in a similar fashion to \code{\link{wave_points}} - the main
difference is that the output values are plotted. Consequently, the set of targets is required
to overlay the region of interest onto the plot.

To ensure that the wave numbers provided in the legend match, one should provide waves
as a list of data.frames with the earliest wave at the start of the list.

The parameters \code{which_wave} and \code{upper_scale} control the level of `zoom' on
each of the lower-triangular and upper-triangular plots, respectively. For the lower
plots, \code{which_wave} determines which of the provided waves is to be used to determine
the output ranges to plot with respect to: generally, higher \code{which_wave} values
result in a more zoomed-in plot. For the upper plots, \code{upper_scale} determines the
plot window via a multiple of the target bounds: higher values result in a more zoomed-out
plot. If not provided, these default to \code{which_wave=0} (or 1 if no wave 0 is given)
and \code{upper_scale = 1}. If the value provided to \code{which_wave} does not correspond
to a provided wave (or one explicitly not included in \code{wave_numbers}), it defaults to
the closest available wave to the value of \code{which_wave}.

If \code{ems} is provided, it should follow the same structure as \code{waves}: at the very
least, it should contain all emulators trained over the course of the waves. The emulator
predictions for a target are made by the emulator for that target whose ranges are the
smallest such that contain the point.
}
\examples{
 wave_values(SIRMultiWaveData, SIREmulators$targets, surround = TRUE, p_size = 1)
 \donttest{
   wave_values(SIRMultiWaveData, SIREmulators$targets, c('nS', 'nI'), l_wid = 0.8)
     wave_values(SIRMultiWaveData, SIREmulators$targets, l_wid = 0.8,
      wave_numbers = c(0, 1, 3), which_wave = 2, upper_scale =  1.5)
     # For many plots, it may be helpful to manually modify the font size
     wave_values(SIRMultiWaveData, SIREmulators$targets) +
      ggplot2::theme(text = ggplot2::element_text(size = 5))
 }
}
\seealso{
Other visualisation tools: 
\code{\link{behaviour_plot}()},
\code{\link{diagnostic_wrap}()},
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
\code{\link{wave_points}()}
}
\concept{visualisation tools}
