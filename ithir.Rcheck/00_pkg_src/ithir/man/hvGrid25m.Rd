\name{hvGrid25m}
\docType{data}
\alias{hvGrid25m}
\title{Raster Grid of the Lower Hunter Valley, NSW, Australia}
\description{
A raster layer (GeoTIFF) representing pixel indices over the Lower Hunter Valley, NSW, Australia. The grid is aligned at 25 m × 25 m resolution and is commonly used as a spatial index for simulations and soil mapping exercises.
}
\format{
\code{hvGrid25m} is a single-layer raster (GeoTIFF format) with:
\itemize{
  \item 860 rows
  \item 676 columns
  \item Pixel resolution of 25 m × 25 m
}
Each pixel contains a numeric index value. The coordinate reference system is \code{GDA94 / MGA Zone 56 (EPSG:28356)}.
}
\details{
The grid covers the Hunter Wine Country Private Irrigation District (HWCPID) in the Lower Hunter Valley (32.83°S, 151.35°E), approximately 140 km north of Sydney. The region spans about 220 km² and features a temperate humid climate with ~750 mm annual rainfall. Dominant land uses include viticulture and dryland grazing.

This grid is used in digital soil mapping exercises, including those taught through the "Using R for Digital Soil Mapping" course. It is useful for simulations, uncertainty mapping, and raster-based modeling workflows.
}
\references{
\itemize{
  \item Malone, B.P., Hughes, P., McBratney, A.B., Minasny, B. (2014). \href{http://dx.doi.org/10.1016/j.geodrs.2014.08.001}{A model for the identification of terrons in the Lower Hunter Valley, Australia}. \emph{Geoderma Regional}, 1, 31–47.
}
}

\examples{
library(ithir)
library(terra)

# Load and plot the raster grid
hv.grid <- rast(system.file("extdata/hvGrid25m_grid.tif", package = "ithir"))
plot(hv.grid, main = "Hunter Valley 25m Grid Index")
}
\keyword{datasets}
