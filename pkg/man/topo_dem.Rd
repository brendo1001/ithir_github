\name{topo_dem}
\docType{data}
\alias{digital elevation}
\title{matrix of digital elevation}
\description{A small matrix that is representative of a digital elevation model. Row and column positions take the place of spatial coordinates}
\usage{data(topo_dem)}
\format{
\code{topo_dem} is an 109 row \code{matrix} with 110 columns. Values correspond to elevation. 
}
\details{This dataset is used to exemplfy the procedure for generating random catena or toposequences which is described in the Using R for Digital Soil Mapping book.}

\references{
\itemize{
\item This data is copyright of the Soil Security Lab, The University of Sydney. If any part of this data is to be used in any publication or report, please provide a citation:
Soil Security Laboratory, 2016. Use R for Digital Soil Mapping Manual. The University of Sydney, Sydney, Australia.
}
}

\examples{
library(ithir)


data(topo_dem)
str(topo_dem)

}
\keyword{datasets}
