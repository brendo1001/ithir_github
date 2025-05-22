\name{ea_spline}
\docType{methods}
\alias{ea_spline}
\title{Fit a Mass-Preserving Spline to Soil Profile Data}
\description{
Fits a continuous mass-preserving spline to numeric soil profile data, typically horizon-structured, such as soil organic carbon or pH. The spline ensures that the integral (or mass) of the fitted curve over each horizon matches the observed values, enabling harmonisation to standard depths.
}
\usage{
ea_spline(obj, var.name, lam = 0.1, d = c(0, 5, 15, 30, 60, 100, 200),
          vlow = 0, vhigh = 1000, show.progress = TRUE)
}
\arguments{
  \item{obj}{An object of class \code{data.frame} or \code{SoilProfileCollection}. Must contain top and bottom depth values and a target numeric variable.}
  \item{var.name}{Character string; name of the target numeric variable to fit the spline on.}
  \item{lam}{Numeric; smoothing parameter (lambda). Lower values give smoother splines.}
  \item{d}{Numeric vector of standard depths at which to estimate the harmonised values.}
  \item{vlow}{Numeric; lower bound. Values below this will be set to \code{vlow}.}
  \item{vhigh}{Numeric; upper bound. Values above this will be set to \code{vhigh}.}
  \item{show.progress}{Logical; if \code{TRUE}, displays a progress bar.}
}
\value{
A list with the following components:
\describe{
  \item{\code{harmonised}}{A data frame of spline-estimated values at the specified standard depths.}
  \item{\code{obs.preds}}{A data frame of observed values and corresponding spline predictions by profile and depth.}
  \item{\code{var.1cm}}{A matrix of fitted values at 1 cm increments, for each profile.}
  \item{\code{splineFitError}}{A data frame with root mean square error (RMSE) and RMSE divided by IQR for each profile.}
}
}
\note{
The function requires at least two horizons with valid numeric values for spline fitting. Profiles with only one horizon are handled using value replication (no spline fitting). Only positive values for top and bottom depths are accepted. Horizon values are assumed to be block-supported averages.
}
\author{Brendan Malone}
\references{
\itemize{
  \item Bishop, T.F.A., McBratney, A.B., Laslett, G.M. (1999). \href{http://dx.doi.org/10.1016/S0016-7061(99)00003-8}{Modelling soil attribute depth functions with equal-area quadratic smoothing splines}. Geoderma, 91(1–2), 27–45.
  \item Malone, B.P., McBratney, A.B., Minasny, B., Laslett, G.M. (2009). \href{http://dx.doi.org/10.1016/j.geoderma.2009.10.007}{Mapping continuous depth functions of soil carbon storage and available water capacity}. Geoderma, 154(1–2), 138–152.
}
}
\examples{
# Example using a simple data.frame
data(oneProfile)
str(oneProfile)
sp_fit <- ea_spline(obj = oneProfile, var.name = "C.kg.m3.")

# Example using a SoilProfileCollection from the aqp package
# library(aqp)
# library(plyr)
# lon <- 3.90; lat <- 7.50; id <- "ISRIC:NG0017"
# top <- c(0, 18, 36, 65, 87, 127)
# bottom <- c(18, 36, 65, 87, 127, 181)
# ORCDRC <- c(18.4, 4.4, 3.6, 3.6, 3.2, 1.2)
# munsell <- c("7.5YR3/2", "7.5YR4/4", "2.5YR5/6", "5YR5/8", "5YR5/4", "10YR7/3")
# prof1 <- join(data.frame(id, top, bottom, ORCDRC, munsell),
#               data.frame(id, lon, lat), type = 'inner')
# depths(prof1) <- id ~ top + bottom
# site(prof1) <- ~ lon + lat
# ORCDRC.s <- ea_spline(prof1, var.name = "ORCDRC")
# str(ORCDRC.s)
}
\keyword{methods}
