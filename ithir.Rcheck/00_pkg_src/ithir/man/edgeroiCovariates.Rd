\name{edgeroiCovariates}
\docType{data}
\alias{edgeroiCovariates}
\alias{edgeroiCovariates_elevation.tif}
\alias{edgeroiCovariates_twi.tif}
\alias{edgeroiCovariates_radK.tif}
\alias{edgeroiCovariates_landsat_b3.tif}
\alias{edgeroiCovariates_landsat_b4.tif}
\title{Environmental Covariate Rasters for the Full Edgeroi District, NSW}
\description{
GeoTIFF raster files of selected environmental covariates covering the full Edgeroi District in northern New South Wales, Australia. These spatial layers are commonly used in digital soil mapping as predictor variables.
}
\format{
This dataset consists of five raster layers (GeoTIFF format) with a pixel resolution of 90 m. The following files are included:
\describe{
  \item{\code{edgeroiCovariates_elevation.tif}}{Ground elevation (bare earth DEM)}
  \item{\code{edgeroiCovariates_twi.tif}}{Topographic wetness index derived from the DEM}
  \item{\code{edgeroiCovariates_radK.tif}}{Gamma radiometric potassium (K) dose rate}
  \item{\code{edgeroiCovariates_landsat_b3.tif}}{Band 3 reflectance from Landsat 7 ETM+}
  \item{\code{edgeroiCovariates_landsat_b4.tif}}{Band 4 reflectance from Landsat 7 ETM+}
}
}
\details{
The Edgeroi District is a well-studied alluvial cropping area situated on the Namoi River plain. It has been the focus of numerous soil investigations, including systematic sampling on a triangular grid (McGarry et al., 1989) and subsequent digital soil mapping applications (Malone et al., 2009). These rasters provide spatial covariates that are commonly used for soil attribute modelling.
}
\note{
Raster files are stored in the \code{inst/extdata} directory of the package and can be accessed with \code{system.file()}. For example:

\code{system.file("extdata/edgeroiCovariates_elevation.tif", package = "ithir")}

All rasters use WGS84 / UTM Zone 55.
}
\references{
\itemize{
  \item Malone, B.P., McBratney, A.B., Minasny, B. (2009). \href{http://dx.doi.org/10.1016/j.geoderma.2009.10.007}{Mapping continuous depth functions of soil carbon storage and available water capacity}. Geoderma, 154, 138–152.
  \item McGarry, D., Ward, W.T., McBratney, A.B. (1989). \emph{Soil Studies in the Lower Namoi Valley: Methods and Data. The Edgeroi Data Set}. (2 vols). CSIRO Division of Soils: Adelaide.
}
}
\examples{
library(terra)

# Load elevation layer
elev_path <- system.file("extdata/edgeroiCovariates_elevation.tif", package = "ithir")
elev_rast <- rast(elev_path)
plot(elev_rast, main = "Edgeroi Elevation Map")
}
\keyword{datasets}
