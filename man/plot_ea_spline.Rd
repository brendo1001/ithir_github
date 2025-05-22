\name{plot_ea_spline}
\docType{methods}
\alias{plot_ea_spline}
\title{Plot Soil Profile Outputs from \code{ea_spline}}
\description{
A plotting function for visualizing the fitted spline outputs from the \code{ea_spline} function. It displays observed horizon data alongside continuous spline fits and/or spline-averaged values at harmonised depth intervals.
}
\usage{
plot_ea_spline(splineOuts, d = t(c(0, 5, 15, 30, 60, 100, 200)), maxd, 
               type = 1, label = "", plot.which = 1)
}
\arguments{
  \item{splineOuts}{A list; output from the \code{ea_spline} function containing harmonised spline results and observed/predicted values.} 
  \item{d}{Numeric vector; harmonised standard depths used in the spline fitting.}
  \item{maxd}{Numeric; maximum soil depth (cm) to be displayed on the plot.}
  \item{type}{Integer (1, 2, or 3); controls which elements to display:
    \describe{
      \item{1}{Observed values and the continuous spline (default).}
      \item{2}{Observed values and spline averages at standard depths.}
      \item{3}{Observed values, spline averages, and the continuous spline.}
    }
  }
  \item{label}{Character string; optional label for the x-axis of the plot.}
  \item{plot.which}{Integer; index of the profile to plot (when multiple profiles were processed by \code{ea_spline}).}
}
\value{
Generates a base R plot showing a single soil profile with user-specified spline elements overlaid. The function is designed to visualize one profile at a time.
}

\note{
This is a companion function for \code{ea_spline}. For plotting multiple profiles, this function should be called iteratively or wrapped in a custom loop or apply function.
}

\author{Brendan Malone}

\examples{
library(ithir)
library(aqp)

# Load example data and convert to SoilProfileCollection
data(oneProfile)
depths(oneProfile) <- Soil.ID ~ Upper.Boundary + Lower.Boundary

# Fit spline to the profile
eaFit <- ea_spline(oneProfile, var.name = "C.kg.m3.", 
                   d = t(c(0, 5, 15, 30, 60, 100, 200)), 
                   lam = 0.1, vlow = 0, show.progress = FALSE)

# Plot the fitted spline
plot_ea_spline(splineOuts = eaFit, 
               d = t(c(0, 5, 15, 30, 60, 100, 200)), 
               maxd = 200, 
               type = 1, 
               label = "Carbon Density")
}

\keyword{methods}
