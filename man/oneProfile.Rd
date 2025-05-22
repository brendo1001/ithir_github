\name{oneProfile}
\docType{data}
\alias{oneProfile}
\title{Example Soil Profile for Soil Carbon Density}
\description{
A soil information \code{data.frame} representing a single soil profile for the variable soil carbon density. Each row corresponds to a horizon (depth interval).
}
\usage{data(oneProfile)}
\format{
A \code{data.frame} with 8 rows and 4 columns:
\describe{
  \item{\code{id}}{Profile identifier}
  \item{\code{top}}{Upper depth (cm)}
  \item{\code{bottom}}{Lower depth (cm)}
  \item{\code{carbon_density}}{Soil carbon density (e.g., kg/m^2)}
}
}
\details{
This dataset is useful for illustrating profile spline fitting, e.g., with \code{\link{ea_spline}}, and is structured as a typical harmonized soil profile.
}
\references{
\itemize{
  \item Malone, B.P., Minasny, B., McBratney, A.B. (2017). \href{https://link.springer.com/book/10.1007/978-3-319-44327-0}{Using R for Digital Soil Mapping}. Springer, Cham. 262 pp.
}
}
\examples{
data(oneProfile)

# Show profile structure
str(oneProfile)

# Fit spline (if using ithir)
# result <- ea_spline(oneProfile, var.name = "carbon_density")
# result$harmonised
}
\keyword{datasets}
