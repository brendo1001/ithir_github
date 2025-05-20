\name{hvTerronDat}
\docType{data}
\alias{hvTerronDat}
\title{Soil Point Data with Terron Class Labels from the Hunter Valley, NSW, Australia}
\description{
A \code{data.frame} of 1000 sites containing terron class information from various locations in the Hunter Valley, NSW, Australia (32.83°S, 151.35°E). Terrons are soil–landscape entities, conceptually similar to soil classes, but derived from a bottom-up clustering of relevant soil and environmental features. This approach is particularly useful for assessing viticultural suitability and delineating wine-growing subregions. The terron concept is discussed in Carré and McBratney (2005), and implemented in Malone et al. (2014).
}
\usage{
data(hvTerronDat)
}
\format{
\code{hvTerronDat} is a 1000-row \code{data.frame} with the following columns:
\itemize{
  \item \code{east}: Easting coordinate (GDA94 / MGA Zone 56)
  \item \code{north}: Northing coordinate (GDA94 / MGA Zone 56)
  \item \code{terron_class}: A factor with 12 different Terron class labels
}
}
\details{
The data are a random sample of Terron class assignments from the regional classification map published by Malone et al. (2014). Each entry represents a unique site with a corresponding terron label. These data are suitable for use in classification mapping, supervised learning, and regional planning.
}
\references{
\itemize{
  \item Carré, F., McBratney, A.B. (2005). \href{http://dx.doi.org/10.1016/j.geoderma.2005.04.012}{Digital terron mapping}. \emph{Geoderma}, 128(3–4), 340–353.
  \item Malone, B.P., Hughes, P., McBratney, A.B., Minasny, B. (2014). \href{http://dx.doi.org/10.1016/j.geodrs.2014.08.001}{A model for the identification of terrons in the Lower Hunter Valley, Australia}. \emph{Geoderma Regional}, 1, 31–47.
}
}

\examples{
library(ithir)
data(hvTerronDat)
head(hvTerronDat)
table(hvTerronDat$terron_class)
}
\keyword{datasets}
