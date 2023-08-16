\name{Hunter Valley, NSW 100m digital elevation model}
\docType{data}
\alias{Hunter Valley DEM}
\title{Hunter Valley DEM}
\description{A \code{dataframe} that is easily converted to a grid raster. It has 3 columns, with the first two being spatial coordinates, and the third being ground elevation information. The coordinates are a regular grid point pattern with 100m spacing, and when converted to raster, resolve to be a digital elevation model for the Lower Hunter Valley region, NSW, Australia. The CRS of the coordinates is WGS 84 UTM Zone 56. }
\usage{data(HV_dem)}
\format{
\code{HV_dem} is a large \code{dataframe} with spatial coordinates, and elevation information labelled in respective columns. 
}
\details{This \code{dataframe} stores information that depicts a regular raster grid}

\references{
\itemize{
\item Malone, B.P., Minasny, B., McBratney, A.B. (2017) \href{https://link.springer.com/book/10.1007/978-3-319-44327-0}{Using R for Digital Soil Mapping}. Springer Cham. 262 Pages.
}
}

\examples{

## HV DEM
# library(ithir)
# library(terra)
 
# data(HV_dem)
# map<- terra::rast(x = HV_dem, type = "xyz")
# plot(map, main = "Hunter Valley DEM") 

}
\keyword{datasets}
