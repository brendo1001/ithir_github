\name{homosoil_globeDat}
\docType{data}
\alias{homosoil_globeDat}
\title{Global Environmental Data }
\description{A data frame of global (0.5 degree resolution) climatic, litholgoy and topographic data.}
\usage{data(homosoil_globeDat)}
\format{
\code{homosoil_globeDat} is an 62254 row \code{data.frame} with 58 columns.  
}
\details{The basis is a global 0.5 degree resolution grid data of climate, topography, and lithology.

For climate, this consists of variables representing long-term mean monthly and seasonal temperature, rainfall, solar radiation and evapotranspiration data. We also use the DEM representing topography, and lithology, which gives broad information on the parent material. The climate data come from the ERA-40 reanalysis and Climate Research Unit (CRU) dataset. More details on the datasets are available on the website \url{http://www.ipcc-data.org/obs/get_30yr_means.html}. For each of the 4 climatic variables (rainfall, temperature, solar radiation and evapotranspiration), we calculated 13 indicators: annual mean, mean for the driest month, mean at the wettest month, annual range, driest quarter mean, wettest quarter mean, coldest quarter mean, hottest quarter mean, lowest ET quarter mean, highest ET quarter mean, darkest quarter mean, lightest quarter mean, and seasonality. From this analysis and including the acquired data, 52 global climatic variables were composed.

The DEM is from the Hydro1k dataset supplied from the USGS (\url{https://lta.cr.usgs.gov/HYDRO1K}), which includes the mean elevation, slope, and compound topgraphic index (CTI).

The lithology is from a global digital map (Durr et al., 2005)  with 7 values which represent the different broad groups of parent materials. The lithology classes are: non- or semi-consolidated sediments, mixed consolidated sediments, silic-clastic sediments, acid volcanic rocks, basic volcanic rocks, complex of metamorphic and igneous rocks, and complex lithology.}

\references{
\itemize{
  \item Durr, H.H., Meybeck, M., and Durr, S.H. (2005) \href{http://dx.doi.org/10.1029/2005GB002515}{Lithologic composition of the Earth's continental surfaces derived from a new digital map emphasizing riverine material transfer}. Global Biogeochemical Cycles, 19, GB4S10. 
}
}


\examples{
library(ithir)


data(homosoil_globeDat)
str(homosoil_globeDat)

}
\keyword{datasets}
