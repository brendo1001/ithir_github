\name{ea_rasSp_fast}
\docType{methods}
\alias{ea_rasSp_fast}
\title{Fast Mass-Preserving Spline on Raster Soil Data}
\description{
Fits a mass-preserving spline model to a multi-layer soil property \code{SpatRaster} and returns raster layers containing spline-averaged values over user-defined depth intervals.
This function uses precomputed spline matrices for efficiency and applies the spline cell-wise using \code{terra::app()}.
}
\usage{
ea_rasSp_fast(obj, lam = 0.1, dIn, dOut, vlow = 0, vhigh = 100, depth_res = 1)
}
\arguments{
  \item{obj}{object of class \code{"SpatRaster"} with one layer per input depth interval.}
  \item{lam}{spline stiffness parameter (\eqn{\lambda}); smaller values = smoother splines.}
  \item{dIn}{numeric vector of input depth boundaries (e.g. \code{c(0,5,15,...)}).}
  \item{dOut}{numeric vector defining output depth intervals (e.g. \code{c(0,30,60)}).}
  \item{vlow}{minimum bound to truncate spline predictions (default = 0).}
  \item{vhigh}{maximum bound to truncate spline predictions (default = 100).}
  \item{depth_res}{numeric; depth interpolation resolution (e.g. 1 = 1cm, 5 = 5cm).}
}
\value{
Returns a \code{SpatRaster} with one layer per output depth interval.
}
\author{Brendan Malone}
\examples{
# Not run on CRAN due to file size and download time
\dontrun{
library(terra)

# Define SLGA V2 clay layer URLs (0–200 cm depth range)
clay_urls <- c(
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif',
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_005_015_EV_N_P_AU_TRN_N_20210902.tif',
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_015_030_EV_N_P_AU_TRN_N_20210902.tif',
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_030_060_EV_N_P_AU_TRN_N_20210902.tif',
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_060_100_EV_N_P_AU_TRN_N_20210902.tif',
  '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_100_200_EV_N_P_AU_TRN_N_20210902.tif'
)

# Load and crop a small extent near Canberra
clay_stack <- rast(clay_urls)
aoi <- ext(149.00, 149.10, -36.00, -35.90)
clay_crop <- crop(clay_stack, aoi)

# Fit spline and generate interpolated output
out <- ea_rasSp_fast(
  obj = clay_crop,
  lam = 0.1,
  dIn = c(0, 5, 15, 30, 60, 100, 200),
  dOut = c(0, 30, 60),
  depth_res = 2
)

# Plot the result
plot(out)
}
}
\keyword{methods}
