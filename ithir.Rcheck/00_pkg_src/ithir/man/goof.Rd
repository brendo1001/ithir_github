\name{goof}
\docType{methods}
\alias{goof}
\title{Goodness of Fit Measures}
\description{
Computes a suite of goodness-of-fit statistics for model evaluation, including root mean square error (RMSE), mean square error (MSE), prediction bias, coefficient of determination (R-squared), concordance correlation coefficient (CCC), ratio of performance to deviation (RPD), ratio of performance to interquartile distance (RPIQ), and residual variance estimates. These metrics are widely used in digital soil mapping and chemometric modelling to assess model accuracy and reliability.
}
\usage{
goof(observed, predicted, plot.it = FALSE, type = "DSM")
}
\arguments{
  \item{observed}{Numeric; a vector or matrix of observed values.}
  \item{predicted}{Numeric; a vector or matrix of predicted values.}
  \item{plot.it}{Logical; if \code{TRUE}, an xy-plot of observed vs predicted values is generated.}
  \item{type}{Character; choose either \code{"DSM"} (default) or \code{"spec"} to tailor the output to digital soil mapping or soil spectroscopy use cases.}
}
\value{
A \code{data.frame} containing the following goodness-of-fit statistics:
\describe{
  \item{\code{R2}}{Coefficient of determination}
  \item{\code{concordance}}{Lin’s concordance correlation coefficient}
  \item{\code{MSE}}{Mean square error}
  \item{\code{RMSE}}{Root mean square error}
  \item{\code{bias}}{Prediction bias (mean error)}
  \item{\code{MSEc}}{Mean residual variance}
  \item{\code{RMSEc}}{Root mean residual variance}
  \item{\code{RPD}}{Ratio of performance to standard deviation}
  \item{\code{RPIQ}}{Ratio of performance to interquartile range}
}
If \code{plot.it = TRUE}, a scatter plot of predicted vs observed values is also returned.
}
\note{
These goodness-of-fit statistics are commonly used across various modelling domains, including but not limited to digital soil mapping and chemometric spectroscopy.
}
\author{Brendan Malone}
\references{
\itemize{
  \item Bellon-Maurel, V., Fernandez-Ahumada, E., Palagos, B., Roger, J., McBratney, A.B. (2010). \href{http://dx.doi.org/10.1016/j.trac.2010.05.006}{Critical review of chemometric indicators commonly used for assessing the quality of the prediction of soil attributes by NIR spectroscopy}. \emph{Trends in Analytical Chemistry}, 29(9), 1073–1081.
  \item Hastie, T., Tibshirani, R., Friedman, J. (2009). \emph{The Elements of Statistical Learning}. Springer Series in Statistics.
}
}
\examples{
library(ithir)
library(MASS)

# Load sample soil data
data(USYD_soil1)

# Fit a linear model
mod.1 <- lm(CEC ~ clay, data = USYD_soil1, y = TRUE, x = TRUE)

# Calculate goodness-of-fit statistics and plot
goof(observed = mod.1$y, predicted = mod.1$fitted.values, plot.it = TRUE)
}
\keyword{methods}
