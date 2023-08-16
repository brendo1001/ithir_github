\name{bbRaster}
\docType{methods}
\alias{bbRaster}
\title{Get the bounding box information of a \code{SpatRaster} from its extents.}
\description{This is a simple function that returns a 4 by 2 matrix of a \code{SpatRaster} bounding box}

\usage{bbRaster(obj)
}
\arguments{
  \item{obj}{object of class \code{"SpatRaster"}} 
}
\value{Returns a 4  x 2 matrix with each row indicating a coordinate pair of the bounding box
}

\author{Brendan Malone}


\examples{

# library(terra)
# target <- rast(system.file("extdata/edgeTarget_C.tif", package="ithir"))
# target
# bbRaster(target)

}
\keyword{methods}
