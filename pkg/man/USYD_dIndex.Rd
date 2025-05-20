\name{USYD_dIndex}
\docType{data}
\alias{USYD_dIndex}
\title{Soil Drainage Index Observations and Predictions}
\description{
A \code{data.frame} of 446 records, each representing a location with observed and predicted values for a soil drainage index. Locations are not georeferenced in this dataset.
}
\usage{data(USYD_dIndex)}
\format{
A \code{data.frame} with 446 rows and 2 columns:
\describe{
  \item{\code{DI_observed}}{Observed soil drainage index (unitless or categorical)}
  \item{\code{DI_predicted}}{Predicted drainage index from a model}
}
}
\details{
This dataset is typical of model validation data in digital soil mapping, where observed values are compared to spatial model predictions. The dataset can be used to assess classification accuracy, correlation, or model bias.
}
\references{
\itemize{
  \item Malone, B.P., Minasny, B., McBratney, A.B. (2017). \href{https://link.springer.com/book/10.1007/978-3-319-44327-0}{Using R for Digital Soil Mapping}. Springer, Cham.
  \item Malone, B.P., McBratney, A.B., Minasny, B. (2018). \href{https://doi.org/10.7717/peerj.4659}{Description and spatial inference of soil drainage using matrix soil colours in the Lower Hunter Valley, New South Wales, Australia}. PeerJ, 6, e4659.
}
}
\examples{
data(USYD_dIndex)

# View summary statistics
summary(USYD_dIndex)

# Plot observed vs predicted
with(USYD_dIndex, plot(DI_observed, DI_predicted,
     main = "Observed vs Predicted Drainage Index",
     xlab = "Observed", ylab = "Predicted"))
}
\keyword{datasets}
