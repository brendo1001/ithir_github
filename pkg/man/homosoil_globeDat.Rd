\name{homosoil_globeDat}
\docType{data}
\alias{homosoil_globeDat}
\title{Global Environmental Data}
\description{
A data frame containing global-scale (0.5° resolution) climatic, lithological, and topographic variables.
}
\usage{data(homosoil_globeDat)}
\format{
\code{homosoil_globeDat} is a \code{data.frame} with 62,254 rows and 58 columns.
}
\details{
The dataset is structured on a global 0.5° resolution grid, integrating climate, topography, and lithology data.

\strong{Climate:} Sourced from ERA-40 reanalysis and CRU datasets, the climate variables include long-term monthly and seasonal averages of temperature, rainfall, solar radiation, and evapotranspiration. For each variable, 13 derived indicators were calculated: annual mean, wettest/driest month, annual range, wettest/driest quarter, hottest/coldest quarter, highest/lowest ET quarters, darkest/lightest quarters, and seasonality — yielding 52 total climatic predictors. More info: \url{http://www.ipcc-data.org/obs/get_30yr_means.html}

\strong{Topography:} Derived from the Hydro1k DEM (USGS), which includes elevation, slope, and the compound topographic index (CTI). Source: \url{https://lta.cr.usgs.gov/HYDRO1K}

\strong{Lithology:} Based on the global lithological map of Dürr et al. (2005), with 7 categorical values representing parent material groups:
\enumerate{
  \item Non-/semi-consolidated sediments
  \item Mixed consolidated sediments
  \item Siliciclastic sediments
  \item Acid volcanic rocks
  \item Basic volcanic rocks
  \item Metamorphic and igneous complexes
  \item Complex lithology
}
}
\references{
\itemize{
  \item Dürr, H.H., Meybeck, M., and Dürr, S.H. (2005). \href{http://dx.doi.org/10.1029/2005GB002515}{Lithologic composition of the Earth's continental surfaces derived from a new digital map emphasizing riverine material transfer}. \emph{Global Biogeochemical Cycles}, 19, GB4S10.
}
}
\examples{
library(ithir)
data(homosoil_globeDat)

# Structure of the dataset
str(homosoil_globeDat)

# Summary of the first few climatic variables
summary(homosoil_globeDat[, 1:6])
}
\keyword{datasets}
