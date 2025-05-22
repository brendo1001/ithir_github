\name{HV_subsoilpH}
\docType{data}
\alias{HV_subsoilpH}
\title{Hunter Valley Subsoil pH Data with Environmental Covariates}
\description{
A \code{data.frame} containing 506 observations of soil pH from the Lower Hunter Valley, NSW, Australia. Each record represents the 60–100 cm depth interval and is associated with a set of intersected environmental covariates derived from digital elevation and Landsat data.
}
\usage{data(HV_subsoilpH)}

\format{
\code{HV_subsoilpH} is a 506-row \code{data.frame}. The first two columns represent spatial coordinates in WGS 84 / UTM Zone 56. Soil pH for the 60–100 cm depth interval is stored in the next column, followed by several environmental covariates:

\describe{
  \item{\code{Terrain_Ruggedness_Index}}{Quantifies topographic heterogeneity. High values indicate rugged terrain.}
  \item{\code{AACN}}{Elevation above channel network base level. Requires known stream network.}
  \item{\code{Landsat_Band1}}{Reflectance from Landsat 7 ETM+ Band 1 (0.45–0.52 µm).}
  \item{\code{Elevation}}{Ground elevation (m), derived from a DEM.}
  \item{\code{Hillshading}}{Hillshade based on DEM and fixed sun angle.}
  \item{\code{Light_insolation}}{Potential solar radiation, modeled over a year with 5-day intervals.}
  \item{\code{Mid_Slope_Position}}{Relative classification of slope position (valley to crest).}
  \item{\code{MRVBF}}{Multi-resolution Valley Bottom Flatness Index (Gallant & Dowling, 2003).}
  \item{\code{NDVI}}{Normalized Difference Vegetation Index from Landsat (B4-B3)/(B4+B3).}
  \item{\code{TWI}}{Topographic Wetness Index.}
  \item{\code{Slope}}{Slope angle in degrees (first derivative of elevation).}
}
}

\details{
The dataset covers the Hunter Wine Country Private Irrigation District (HWCPID), located in the Lower Hunter Valley (32.83°S, 151.35°E), about 140 km north of Sydney. The HWCPID spans ~220 km², with a temperate humid climate (~750 mm annual rainfall). The region supports extensive viticulture and dryland grazing. This dataset is used in exercises for quantifying prediction uncertainty in digital soil mapping.
}

\note{This dataset is used in the "Using R for Digital Soil Mapping" course for exercises on uncertainty quantification.}

\references{
\itemize{
  \item Gallant, J.C., Dowling, T.I. (2003). \href{http://dx.doi.org/10.1029/2002WR001426}{A multiresolution index of valley bottom flatness for mapping depositional areas}. \emph{Water Resources Research}, 39(12), 1347.
  \item Malone, B.P., Hughes, P., McBratney, A.B., Minasny, B. (2009). \href{http://dx.doi.org/10.1016/j.geodrs.2014.08.001}{A model for the identification of terrons in the Lower Hunter Valley, Australia}. \emph{Geoderma Regional}, 1, 31–47.
  \item Malone, B.P., Minasny, B., McBratney, A.B. (2017). \href{https://link.springer.com/book/10.1007/978-3-319-44327-0}{Using R for Digital Soil Mapping}. Springer Cham, 262 pp.
}
}

\examples{
library(ithir)
data(HV_subsoilpH)
summary(HV_subsoilpH)
}
\keyword{datasets}
