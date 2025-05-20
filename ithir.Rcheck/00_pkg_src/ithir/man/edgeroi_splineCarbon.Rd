\name{edgeroi_splineCarbon}
\docType{data}
\alias{edgeroi_splineCarbon}
\title{Harmonised Soil Carbon Density Data from the Edgeroi District, NSW}
\description{
A \code{data.frame} containing soil carbon density estimates at 341 locations across the Edgeroi District, NSW, Australia (approx. 30.11°S, 149.66°E). Soil carbon values were derived using a pedotransfer function based on measured soil carbon concentration and texture. The resulting values were harmonised using mass-preserving splines according to the GlobalSoilMap specifications.
}
\usage{data(edgeroi_splineCarbon)}
\format{
A \code{data.frame} with 341 rows. Columns include:
\describe{
  \item{\code{id}}{Profile identifier}
  \item{\code{east}}{Easting (UTM Zone 55)}
  \item{\code{north}}{Northing (UTM Zone 55)}
  \item{\code{carbon_0_5cm}}{Estimated carbon density (0–5 cm)}
  \item{\code{carbon_5_15cm}}{Estimated carbon density (5–15 cm)}
  \item{\code{carbon_15_30cm}}{Estimated carbon density (15–30 cm)}
  \item{\code{carbon_30_60cm}}{Estimated carbon density (30–60 cm)}
  \item{\code{carbon_60_100cm}}{Estimated carbon density (60–100 cm)}
  \item{\code{carbon_100_200cm}}{Estimated carbon density (100–200 cm)}
  \item{\code{max_depth}}{Maximum spline depth fitted for each profile (cm)}
}
}
\details{
Sampling was based on a design described by McGarry et al. (1989), with 210 samples taken from a systematic equilateral triangular grid (2.8 km spacing), and 131 additional samples located along transects or more irregularly. Spline harmonisation was applied to align depths with GlobalSoilMap conventions.
}
\references{
\itemize{
  \item Arrouays, D., McKenzie, N., Hempel, J., Richer de Forges, A., and McBratney, A. (eds) (2014). \emph{GlobalSoilMap: Basis of the Global Spatial Soil Information System}. CRC Press.
  \item McGarry, D., Ward, W.T., McBratney, A.B. (1989). \emph{Soil Studies in the Lower Namoi Valley: Methods and Data. The Edgeroi Data Set}. (2 vols). CSIRO Division of Soils, Adelaide.
}
}
\examples{
library(ithir)
library(sf)

# Load the dataset
data(edgeroi_splineCarbon)

# Convert to sf object and plot
spat_edgeroi <- st_as_sf(edgeroi_splineCarbon, coords = c("east", "north"), crs = 32755)
plot(spat_edgeroi, pch = 16, cex = 0.5, main = "Edgeroi Soil Carbon Point Locations")
}
\keyword{datasets}
