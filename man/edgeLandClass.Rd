\name{edgeLandClass}
\docType{data}
\alias{edgeLandClass}
\title{Land Classification Points from the Edgeroi District, NSW}
\description{
A \code{data.frame} of 500 point locations in the Edgeroi District, NSW, Australia (approx. 30.11°S, 149.66°E), each assigned to one of six estimated land classes. The land classification was derived via unsupervised classification of Landsat 7 ETM+ spectral data (acquisition date unknown). These points represent a random sample from a classified land use map.
}
\usage{data(edgeLandClass)}
\format{
A \code{data.frame} with 500 rows and 3 columns:
\describe{
  \item{\code{x}}{Easting (UTM Zone 55)}
  \item{\code{y}}{Northing (UTM Zone 55)}
  \item{\code{LandClass}}{Integer values from 1 to 6 representing land classes: 
    \enumerate{
      \item Dense forest
      \item Open forest
      \item Water bodies
      \item Woody vegetation and native grassland
      \item Irrigated cropping
      \item Dryland cropping
    }
  }
}
}
\details{
This dataset represents typical point-based environmental classification data and is suitable for supervised modelling, classification accuracy assessment, or land use interpretation.
}
\references{
\itemize{
  \item Soil Security Laboratory, 2015. \emph{Use R for Digital Soil Mapping Manual}. The University of Sydney, Sydney, Australia.
}
}
\examples{
data(edgeLandClass)

# View land class summary
summary(edgeLandClass$LandClass)

# Plot the locations by class
plot(edgeLandClass$x, edgeLandClass$y, col = edgeLandClass$LandClass,
     pch = 20, main = "Edgeroi Land Classification Points")
}
\keyword{datasets}
