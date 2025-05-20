\name{USYD_soil1}
\docType{data}
\alias{USYD_soil1}
\title{Random Selection of Soil Profile Data from New South Wales}
\description{
A soil information \code{data.frame} containing 166 horizon observations across 29 soil profiles. Each row corresponds to a specific depth interval within a profile. Soil attributes such as organic carbon, pH, and electrical conductivity are included where available.
}
\usage{data(USYD_soil1)}
\format{
A \code{data.frame} with 166 rows and several columns:
\describe{
  \item{\code{id}}{Profile identifier}
  \item{\code{top}}{Upper depth (cm)}
  \item{\code{bottom}}{Lower depth (cm)}
  \item{\code{SOC}}{Soil organic carbon (\%)}
  \item{\code{pH}}{Soil pH in water (1:5)}
  \item{\code{EC}}{Electrical conductivity (dS/m)}
  \item{...}{Other soil or metadata fields}
}
}
\details{
This dataset is representative of a harmonized soil profile table commonly used for digital soil modelling, and is suitable for use with spline fitting tools such as \code{\link{ea_spline}}.
}
\references{
\itemize{
  \item Malone, B.P., Minasny, B., McBratney, A.B. (2017). \href{https://link.springer.com/book/10.1007/978-3-319-44327-0}{Using R for Digital Soil Mapping}. Springer, Cham. 262 pp.
}
}
\examples{
data(USYD_soil1)

# Overview
summary(USYD_soil1)

# Show profile IDs
table(USYD_soil1$id)
}
\keyword{datasets}
