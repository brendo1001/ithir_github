\name{bbRaster}
\docType{methods}
\alias{bbRaster}
\title{Get Bounding Box Coordinates from a \code{SpatRaster}}
\description{
Returns the bounding box of a \code{SpatRaster} object as a 4 x 2 matrix, with each row representing one corner of the extent in coordinate pairs.
}
\usage{
bbRaster(obj)
}
\arguments{
  \item{obj}{An object of class \code{SpatRaster}, from the \pkg{terra} package.}
}
\value{
A 4 x 2 numeric matrix. Each row corresponds to a corner of the bounding box, in clockwise order starting from the lower left:
\describe{
  \item{Row 1}{Lower-left (xmin, ymin)}
  \item{Row 2}{Lower-right (xmax, ymin)}
  \item{Row 3}{Upper-right (xmax, ymax)}
  \item{Row 4}{Upper-left (xmin, ymax)}
}
}
\author{Brendan Malone}
\examples{
library(terra)

# Example raster
target <- rast(system.file("extdata/edgeTarget_C.tif", package = "ithir"))

# Get bounding box coordinates
bbRaster(target)
}
\keyword{methods}
