\name{edgeTarget_C}
\docType{data}
\alias{edgeTarget_C}
\title{1 km Resolution Soil Carbon Stock Map (Subset) – Edgeroi District, NSW}
\description{
A GeoTIFF raster representing predicted soil organic carbon (SOC) stock for the 0–30 cm depth interval over a small subset of the Edgeroi District, NSW, Australia. The map is at 1 km spatial resolution and was derived using digital soil mapping methods.
}
\format{
A single-layer raster (GeoTIFF) with 7 rows and 10 columns, each cell representing soil carbon stock in the 0–30 cm layer (e.g., Mg C/ha). The spatial resolution is 1 km x 1 km. The raster is projected in WGS84 / UTM Zone 55.
}
\details{
This raster was produced using digital soil mapping techniques described in McBratney et al. (2003), and uses soil profile data from McGarry et al. (1989) combined with environmental covariates. The extent of the raster matches that of the covariate subset rasters provided in this package (see \code{edgeGrids_*}), though the resolution is coarser (1 km).
}
\note{
The raster file is located in the \code{inst/extdata} directory and can be accessed using \code{system.file()}. It shares spatial coverage with the covariates in \code{edgeroi_covariates_subset} but has a coarser resolution.
}
\references{
\itemize{
  \item Malone, B.P., McBratney, A.B., Minasny, B. (2009). \href{http://dx.doi.org/10.1016/j.geoderma.2009.10.007}{Mapping continuous depth functions of soil carbon storage and available water capacity}. Geoderma, 154, 138–152.
  \item McBratney, A.B., Mendonça Santos, M.L., Minasny, B. (2003). \href{http://dx.doi.org/10.1016/S0016-7061(03)00223-4}{On digital soil mapping}. Geoderma, 117, 3–52.
  \item McGarry, D., Ward, W.T., McBratney, A.B. (1989). \emph{Soil Studies in the Lower Namoi Valley: Methods and Data. The Edgeroi Data Set}. (2 vols). CSIRO Division of Soils: Adelaide.
}
}
\examples{
library(terra)

# Load the raster from the package
soc_path <- system.file("extdata/edgeTarget_C.tif", package = "ithir")
soc_raster <- rast(soc_path)

# Plot the raster
plot(soc_raster, main = "Edgeroi SOC Stock (0–30 cm)", col = rev(terrain.colors(20)))
}
\keyword{datasets}
