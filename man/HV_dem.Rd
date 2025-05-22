\name{HV_dem}
\docType{data}
\alias{HV_dem}
\title{Digital Elevation Model of the Hunter Valley, NSW}
\description{
A \code{data.frame} containing spatial coordinates and elevation values representing a regular 100-metre resolution grid over the Lower Hunter Valley region of New South Wales, Australia. The coordinates are in WGS 84 / UTM Zone 56. When converted to a raster, the data form a digital elevation model (DEM).
}
\usage{data(HV_dem)}
\format{
A \code{data.frame} with three columns:
\describe{
  \item{\code{x}}{Easting (UTM Zone 56)}
  \item{\code{y}}{Northing (UTM Zone 56)}
  \item{\code{elevation}}{Elevation in metres above sea level}
}
}
\details{
This dataset can be readily converted to a raster using the \code{terra::rast()} function with \code{type = "xyz"}. It can support terrain analysis and serve as a covariate for digital soil mapping.
}
\references{
\itemize{
  \item Malone, B.P., Minasny, B., McBratney, A.B. (2017). \href{https://link.springer.com/book/10.1007/978-3-319-44327-0}{Using R for Digital Soil Mapping}. Springer, Cham.
}
}
\examples{
library(terra)
data(HV_dem)

# Convert to raster and plot
dem_rast <- terra::rast(x = HV_dem, type = "xyz")
plot(dem_rast, main = "Hunter Valley DEM")
}
\keyword{datasets}
