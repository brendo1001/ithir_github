\name{hvPoints250}
\docType{data}
\alias{hvPoints250}
\title{Hunter Valley Soil Observation Points (n = 250)}
\description{
A \code{data.frame} containing coordinates for 250 randomly selected locations in the Lower Hunter Valley, NSW, Australia.
}
\usage{data(hvPoints250)}
\format{
A \code{data.frame} with 250 rows and 2 columns:
\describe{
  \item{\code{x}}{Easting (UTM Zone 56)}
  \item{\code{y}}{Northing (UTM Zone 56)}
}
}
\details{
The area sampled is the Hunter Wine Country Private Irrigation District (HWCPID), located in the Lower Hunter Valley (approx. 32.83°S, 151.35°E), covering about 220 km^2. The HWCPID lies roughly 140 km north of Sydney in a temperate zone with warm, humid summers and relatively mild, humid winters. Average annual rainfall is approximately 750 mm. Land use includes extensive viticulture, followed by dryland grazing.
}
\references{
\itemize{
  \item Malone, B.P., Hughes, P., McBratney, A.B., Minasny, B. (2009). \href{http://dx.doi.org/10.1016/j.geodrs.2014.08.001}{A model for the identification of terrons in the Lower Hunter Valley, Australia}. Geoderma Regional, 1, 31–47.
}
}
\examples{
data(hvPoints250)

# Summary
summary(hvPoints250)

# Simple plot
plot(hvPoints250, pch = 20, col = "darkgreen", main = "Hunter Valley Sample Points")
}
\keyword{datasets}
