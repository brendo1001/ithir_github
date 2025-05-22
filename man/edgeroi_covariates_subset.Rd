\name{edgeroi_covariates_subset}
\docType{data}
\alias{edgeroi_covariates_subset}
\title{Selected Subset of Environmental Covariates for the Edgeroi District, NSW}
\description{
A set of GeoTIFF rasters representing selected environmental covariates for a small area in the Edgeroi District, New South Wales, Australia.
}
\format{
The rasters are stored as GeoTIFF files with a spatial resolution of 90 m × 90 m. The projection is WGS 84 UTM Zone 55. The following rasters are available:
\itemize{
  \item \code{edgeGrids_Doserate.tif} — Gamma radiometric data.
  \item \code{edgeGrids_Elevation.tif} — Ground elevation derived from a digital elevation model (DEM).
  \item \code{edgeGrids_Panchromat.tif} — Panchromatic band from Landsat 7 ETM+ imagery.
  \item \code{edgeGrids_Slope.tif} — Slope gradient derived from the DEM.
  \item \code{edgeGrids_TWI.tif} — Topographic wetness index (TWI), derived from the DEM.
}
}
\details{
The Edgeroi District is located on the alluvial Namoi River plain and is a highly productive agricultural region. It has been the focus of numerous soil investigations, including those described in McGarry et al. (1989). This dataset supports digital soil mapping (DSM) applications and is a spatial subset of a more extensive covariate dataset for the district.
}
\note{
The original spatial data used to generate these rasters were sourced from public repositories maintained by CSIRO, Geoscience Australia, and NASA. All rasters use the WGS 84 UTM Zone 55 coordinate reference system.
}
\references{
\itemize{
  \item Malone, B.P., McBratney, A.B., Minasny, B. (2009). \href{http://dx.doi.org/10.1016/j.geoderma.2009.10.007}{Mapping continuous depth functions of soil carbon storage and available water capacity}. Geoderma, 154, 138–152.
  \item McGarry, D., Ward, W.T., McBratney, A.B. (1989). Soil Studies in the Lower Namoi Valley: Methods and Data. The Edgeroi Data Set. CSIRO Division of Soils.
}
}
\examples{
library(ithir)
library(terra)

# Load and plot the elevation raster
elevation <- rast(system.file("extdata/edgeGrids_Elevation.tif", package = "ithir"))
plot(elevation, main = "Edgeroi Elevation Map")
}
\keyword{datasets}
