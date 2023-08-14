\name{Hunter Valley terron data}
\docType{data}
\alias{hvTerronDat}

\title{Soil point data from the Hunter Valley, NSW, Australia}

\description{A \code{dataframe} of 1000 sites containing information of terron classes from various locations of the Hunter Valley, NSW, Australia (32.83S 151.35E). Terrons are soil and landscape entities, similar to a soil class, yet have been created from a bottom up approach using soil and landscape information deemed necessary for evaluating areas for suitability and/or  differentiation of wine growing regions. The terron concept is described in XXX and brought into implementation in Malone et al. (2014). Specifically, the data are Terron classes as sampled from the map presented in Malone et al. (2014).The sample data contains 1000 entries of which there are 12 different Terron classes. Locations have a recorded coordinate location. The CRS of the points is WGS84 UTM Zone 56.}

\usage{data(hvTerronDat)}
\format{
\code{terron.dat} is a 1000 row \code{dataframe} with locational and terron information labelled in respective columns. 
}
\details{This data frame is a typical soil information table}

\references{
\itemize{
\item Carre, F., McBratney, A.B. (2005) \href{http://dx.doi.org/10.1016/j.geoderma.2005.04.012}{Digital terron mapping}, Geoderma, Volume 128, Issues 3-4, October 2005, Pages 340-353.

\item Malone, B.P., Hughes, P.,  McBratney, A.B., Minasny, B. (2014) \href{http://dx.doi.org/10.1016/j.geodrs.2014.08.001}{A model for the identification of terrons in the Lower Hunter
Valley, Australia}. Geoderma Regional 1, 31-47.


}
}

\examples{
library(ithir)


data(hvTerronDat)
head(hvTerronDat)

}
\keyword{datasets}
