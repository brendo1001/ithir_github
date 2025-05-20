\name{Hunter Valley covariates}
\docType{data}
\alias{hunterCovariates}
\title{Suite of selected environmental covariates for the the Lower Hunter Valley, NSW}
\description{A set of Geotiffs of selected environmental covariates covering the Lower Hunter Valley, NSW.}
\usage{system.file("extdata/hunterCovariates_NAME.tif", package="ithir")}

\format{
These data are Geotiffs of selected environmental covariates. The rasters have a pixel resolution of 25m x 25m. The CRS of each raster is WGS 84 UTM Zone 56. The following raster are available:
 \itemize{
\item \code{AACN}: Difference between elevation and an interpolation of a channel network base level elevation. Knowledge of the spatial distribution of channel networks (lines) is therefore necessary for this parameter.
\item \code{Drainage.Index}: A measure of soil water drainage capability. This index is derived from an empirical model of soil color and soil drainage. 
\item \code{Light_insolation}: Measure of potential incoming solar radiation, and used as a parameter for evaluating the positional aspect effect. Derived from digital elevation model, this parameter was evaluated over the duration of a single calendar year with a 5 day time step.
\item \code{TWI}: A secondary land form parameter which estimates for each pixel, its tendency to accumulate water.
\item \code{Gamma.Total.Count}: Gamma radiometric data. Here the variable is total gamma count im ppm.
}
}

\details{The area in question is the Hunter Wine Country Private Irrigation District (HWCPID), situated in the Lower Hunter Valley, NSW (32.83S 151.35E), and covers an area of approximately 220 km2. The HWCPID is approximately 140 km north of Sydney, NSW, Australia. Climatically, the HWCPID is situated in a temperate climatic zone, and experiences warm humid summers, and relatively cooler yet also humid winters. Rainfall is mostly uniformly distributed throughout the year. On average the HWCPID receives just over 750 mm of rainfall annually. In terms of land use, an expansive viticultural industry is situated in the area and is possibly most widespread
of rural industries, followed by dry land agricultural grazing systems.}

\note{These Geotiffs are used in the Using R for Digital Soil Mapping course of exercises.}

\references{
\itemize{
\item Malone, B.P., Hughes, P.,  McBratney, A.B., Minasny, B. (2014) \href{http://dx.doi.org/10.1016/j.geodrs.2014.08.001}{A model for the identification of terrons in the Lower Hunter
Valley, Australia}. Geoderma Regional 1, 31-47.
}
}

\examples{

# library(ithir)
# library(terra)

# example load the altitude above channel network grid
# aacn <- rast(system.file("extdata/hunterCovariates_hunterCovariates_AACN.tif", package="ithir"))

# simple plot
# plot(aacn, main= "Hunter Valley AACN Map")

}
\keyword{datasets}
