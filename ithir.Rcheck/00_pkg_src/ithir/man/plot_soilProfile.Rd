\name{plot_soilProfile}
\docType{methods}
\alias{plot_soilProfile}
\title{Plot Soil Profile Data}
\description{
A simple plotting function to visualize soil profile data from a \code{data.frame}. Each horizon is plotted as a bar, with its vertical extent representing the depth interval and its length representing the value of a soil property.
}
\usage{
plot_soilProfile(data, vals, depths, label = "")
}
\arguments{
  \item{data}{A \code{data.frame}; typically a soil profile table containing at least upper and lower depth columns, along with one or more soil property columns.}
  \item{vals}{A numeric vector; the soil property values to be plotted, usually a column extracted from \code{data}.}
  \item{depths}{A \code{data.frame} or \code{matrix} with two columns specifying the upper and lower boundaries of each horizon in the profile.}
  \item{label}{Character string; optional label to display on the x-axis of the plot.}
}
\value{
Returns a plot representing the soil profile, with bars for each horizon whose length corresponds to the supplied variable values, and height corresponds to the depth interval.
}

\note{
This function is intended to plot one soil profile at a time. To visualize multiple profiles, it should be embedded in a loop or wrapper function.
}

\author{Brendan Malone}

\examples{
library(ithir)

# Load example profile data
data(oneProfile)

# Visualize the profile
plot_soilProfile(
  data = oneProfile,
  vals = oneProfile$C.kg.m3.,
  depths = oneProfile[, 2:3],
  label = names(oneProfile)[4]
)
}
\keyword{methods}
