\name{edgeroi covariates (whole district)}
\docType{data}
\alias{edgeroiCovariates}
\title{Suite of selected environmental covariates for the Edgeroi District, NSW}
\description{Geotiffs of selected environmental covariates with near complete coverage of the Edgeroi district in NSW, Australia.}
\usage{system.file("extdata/edgeCovariates_NAME.tif", package="ithir")}
\format{
\code{edgeroiCovariates} consists of 5 Geotiffs of selected environmental covariates for the near entire Edgeroi district NSW, Australia. The grids have a pixel resolution of 90m x 90m. The available rasters are:
  \describe{
  \item{\code{elevation.tif}}{numeric; topographic variable of bare earth ground elevation. Derived from digital elevation model}
  \item{\code{twi.tif}}{numeric; topographic wetness index. Secondary derivative of the digital elevation model}
  \item{\code{radK.tif}}{numeric; gamma radiometric data}
  \item{\code{landsat_b3.tif}}{numeric; band 3 reflectance of the Landsat 7 satelite }  
  \item{\code{landsat_b4.tif}}{numeric; band 4 reflectance of the Landsat 7 satelite }  
}
}
\details{The Edgeroi District, NSW is an intensive cropping area upon the fertile alluvial Namoi River plain. The District has been the subject of many soil invetigations, namely McGarry et al. (1989) whom describe an extensive soil data set collected from the area. More recently, digital soil mapping studies of the area have been conducted, for example, Malone et al. (2009).}
\note{The raw spatial data that contributed to the creation of the Geotiffs were sourced from publically accessable repositories hosted by various Australian Government and international agencies including CSIRO (for the DEM), Geosciences Australia (for the radiometric data) and NASA (for the Landsat 7 ETM+ data). The projection for each raster is WGS 84 Zone 55.}

\references{
\itemize{
\item Malone, B.P., McBratney, A.B., Minasny, B. (2009) \href{http://dx.doi.org/10.1016/j.geoderma.2009.10.007}{Mapping continuous depth functions of soil carbon storage and available water capacity}. Geoderma 154, 138-152.
\item McGarry, D., Ward, W.T., McBratney, A.B. (1989) Soil Studies in the Lower Namoi Valley: Methods and Data. The Edgeroi Data Set. (2 vols) (CSIRO Division of Soils: Adelaide).
}
}
\examples{

# library(ithir)
# library(terra)

# example load the elevation grid
# elevation <- rast(system.file("extdata/edgeroiCovariates_elevation.tif", package="ithir"))

# simple plot
#plot(elevation, main= "Edgeroi Elevation Map")

}
\keyword{datasets}
