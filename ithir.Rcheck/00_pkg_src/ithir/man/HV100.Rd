\name{HV100}
\docType{data}
\alias{HV100}
\title{Soil Point Data from the Hunter Valley, NSW, Australia}
\description{
A soil information \code{data.frame} with 100 observations from various locations in the Hunter Valley, NSW, Australia (approx. 32.83°S, 151.35°E). The data were collected using a stratified random sampling design in 2010, as described in Malone et al. (2011). Each row represents an observation at the 0–5 cm depth interval. Soil attributes include organic carbon, pH (1:5 soil:water), and electrical conductivity. Sites are labeled and georeferenced, with coordinates recorded in WGS84 / UTM Zone 56.
}
\usage{data(HV100)}
\format{
A \code{data.frame} with 100 rows and multiple columns:
\describe{
  \item{\code{Site}}{Site name or identifier}
  \item{\code{X}}{Easting (UTM Zone 56)}
  \item{\code{Y}}{Northing (UTM Zone 56)}
  \item{\code{SOC}}{Soil organic carbon (\%)}
  \item{\code{pH}}{Soil pH in water (1:5)}
  \item{\code{EC}}{Electrical conductivity (dS/m)}
  \item{...}{Other relevant soil or site variables}
}
}
\details{
This dataset represents a typical harmonized surface soil information table suitable for digital soil mapping and covariate analysis.
}
\references{
\itemize{
  \item Malone, B.P., de Gruijter, J.J., McBratney, A.B., Minasny, B., Brus, D.J. (2011). \href{http://dx.doi.org/10.2136/sssaj2010.0280}{Using Additional Criteria for Measuring the Quality of Predictions and Their Uncertainties in a Digital Soil Mapping Framework}. Soil Science Society of America Journal, 75(3), 1032–1043.
}
}
\examples{
# Load and inspect the dataset
data(HV100)

# Basic structure and summary
str(HV100)
summary(HV100)
}
\keyword{datasets}
