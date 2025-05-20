\name{fit_mpspline_optimized}
\docType{methods}
\alias{fit_mpspline_optimized}
\title{Fit Mass-Preserving Spline to a Single Soil Profile}
\description{
Applies the mass-preserving spline algorithm to a numeric vector of soil profile values and returns averaged values over specified output depth intervals. 
It uses precomputed matrix structures for efficiency.
}
\usage{
fit_mpspline_optimized(vals, spline_info, dOut, vlow, vhigh, depth_res = 1)
}
\arguments{
  \item{vals}{numeric vector of input values across input soil depth intervals (e.g., from one raster pixel).}
  \item{spline_info}{list output from \code{\link{precompute_spline_structures}} that contains spline matrix components.}
  \item{dOut}{numeric vector defining output depth intervals (e.g., \code{c(0,30,60)}).}
  \item{vlow}{minimum bound to truncate spline predictions (e.g., 0).}
  \item{vhigh}{maximum bound to truncate spline predictions (e.g., 100).}
  \item{depth_res}{numeric; resolution for interpolating spline (e.g., 1 = 1cm steps).}
}
\value{
Returns a numeric vector of spline-averaged values for each specified output depth interval.
}
\author{Brendan Malone}
\examples{
# Not run on CRAN due to external raster data size
\dontrun{
library(terra)

# Define SLGA clay raster URLs
clay_urls <- c(
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif',
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_005_015_EV_N_P_AU_TRN_N_20210902.tif',
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_015_030_EV_N_P_AU_TRN_N_20210902.tif',
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_030_060_EV_N_P_AU_TRN_N_20210902.tif',
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_060_100_EV_N_P_AU_TRN_N_20210902.tif',
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_100_200_EV_N_P_AU_TRN_N_20210902.tif'
)

# Load and crop raster stack
clay_stack <- rast(clay_urls)
aoi <- ext(149.00, 149.10, -36.00, -35.90)
clay_crop <- crop(clay_stack, aoi)

# Extract a single profile (pixel)
vals <- terra::extract(clay_crop, cbind(149.05, -35.95))[1, -1]

# Precompute spline structures
dIn <- c(0, 5, 15, 30, 60, 100, 200)
spline_info <- precompute_spline_structures(dIn, lam = 0.1)

# Fit spline to single profile
fit <- fit_mpspline_optimized(
  vals = vals,
  spline_info = spline_info,
  dOut = c(0, 30, 60),
  vlow = 0,
  vhigh = 100,
  depth_res = 1
)

fit
}
}
\seealso{\code{\link{precompute_spline_structures}}, \code{\link{ea_rasSp_fast}}}
\keyword{methods}
