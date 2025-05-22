\name{hunterCovariates}
\docType{data}
\alias{hunterCovariates}
\title{Environmental Covariate Rasters for the Lower Hunter Valley, NSW}
\description{
A set of GeoTIFF rasters representing environmental covariates across the full extent of the Lower Hunter Valley in New South Wales, Australia. These covariates are commonly used in digital soil mapping and terrain analysis.
}
\format{
These raster layers (GeoTIFF format) are stored in the \code{inst/extdata/} directory of the package. Each raster has a spatial resolution of 25 m × 25 m and uses the WGS 84 UTM Zone 56 coordinate reference system.

Available files, all prefixed with \code{hunterCovariates_A_}, include:

\describe{
  \item{\code{hunterCovariates_A_AACN.tif}}{Elevation above channel network base level. Requires known stream network for interpretation.}
  \item{\code{hunterCovariates_A_elevation.tif}}{Bare-earth digital elevation model (DEM), in meters above sea level.}
  \item{\code{hunterCovariates_A_Hillshading.tif}}{DEM-derived hillshade raster with fixed illumination angle.}
  \item{\code{hunterCovariates_A_light_insolation.tif}}{Annual solar radiation potential derived from DEM (5-day timestep resolution).}
  \item{\code{hunterCovariates_A_MRVBF.tif}}{Multi-resolution valley bottom flatness index for identifying depositional zones.}
  \item{\code{hunterCovariates_A_Slope.tif}}{Slope angle in degrees (first derivative of elevation).}
  \item{\code{hunterCovariates_A_TRI.tif}}{Topographic Ruggedness Index — measures terrain heterogeneity.}
  \item{\code{hunterCovariates_A_TWI.tif}}{Topographic Wetness Index — indicates potential for water accumulation.}
}
}
\details{
These rasters cover the Hunter Wine Country Private Irrigation District (HWCPID), located in the Lower Hunter Valley (approx. 32.83°S, 151.35°E), ~140 km north of Sydney. The area spans about 220 km², features a temperate humid climate, and is dominated by viticulture and dryland grazing systems. These layers are used in digital soil mapping studies and training material.
}
\note{
This raster stack is used in the "Using R for Digital Soil Mapping" training course and is provided for demonstration and educational purposes.
}
\references{
\itemize{
  \item Gallant, J.C., Dowling, T.I. (2003). \href{http://dx.doi.org/10.1029/2002WR001426}{A multiresolution index of valley bottom flatness for mapping depositional areas}. \emph{Water Resources Research}, 39(12), 1347.
  \item Malone, B.P., Hughes, P., McBratney, A.B., Minasny, B. (2014). \href{http://dx.doi.org/10.1016/j.geodrs.2014.08.001}{A model for the identification of terrons in the Lower Hunter Valley, Australia}. \emph{Geoderma Regional}, 1, 31–47.
}
}
\examples{
library(ithir)
library(terra)

# Load and plot the slope raster layer
slope <- rast(system.file("extdata/hunterCovariates_A_Slope.tif", package = "ithir"))
plot(slope, main = "Hunter Valley - Slope")
}
\keyword{datasets}
