\name{topo_dem}
\docType{data}
\alias{topo_dem}
\title{Example Digital Elevation Model as a Matrix}
\description{
A small synthetic digital elevation model (DEM) represented as a numeric \code{matrix}. Row and column indices act as proxy spatial coordinates. Values represent elevation in metres.
}
\usage{data(topo_dem)}
\format{
A numeric \code{matrix} with 109 rows and 110 columns. Each value corresponds to a ground elevation.
}
\details{
This dataset is used to exemplify procedures for generating random catenas or toposequences, as described in the book "Using R for Digital Soil Mapping".
}
\references{
\itemize{
  \item Malone, B.P., Minasny, B., McBratney, A.B. (2017). \href{https://link.springer.com/book/10.1007/978-3-319-44327-0}{Using R for Digital Soil Mapping}. Springer, Cham.
}
}
\examples{
data(topo_dem)

# Basic inspection
str(topo_dem)
image(topo_dem, main = "Synthetic DEM", col = terrain.colors(20))
}
\keyword{datasets}
