\name{goofcat}
\docType{methods}
\alias{goofcat}
\title{Goodness of Fit Measures for Categorical Models}
\description{
This function computes several diagnostic statistics for evaluating the performance of classification models where the target variable is categorical. The metrics include: Overall Accuracy, Producer’s Accuracy, User’s Accuracy, and the Kappa Statistic. These are widely used in remote sensing and classification modelling (see Congalton, 1991).
}
\usage{
goofcat(observed = NULL, predicted = NULL, conf.mat, imp = FALSE)
}
\arguments{
  \item{observed}{A vector of observed class labels (numeric or character).}
  \item{predicted}{A vector of predicted class labels (numeric or character). Must be the same length as \code{observed}.}
  \item{conf.mat}{Optional; a confusion matrix summarizing classification outcomes. Should be a square matrix.}
  \item{imp}{Logical; set to \code{TRUE} if providing a pre-computed confusion matrix via \code{conf.mat}. Default is \code{FALSE}.}
}
\value{
A named \code{list} containing:
\describe{
  \item{\code{OverallAccuracy}}{Proportion of correct classifications}
  \item{\code{ProducersAccuracy}}{Recall or sensitivity per class}
  \item{\code{UsersAccuracy}}{Precision or positive predictive value per class}
  \item{\code{Kappa}}{Kappa statistic (chance-corrected agreement)}
  \item{\code{ConfusionMatrix}}{Returned confusion matrix used in the calculation}
}
}
\note{
If \code{conf.mat} is not of class \code{matrix} or is not square in dimension, the function will halt with an error.
}
\author{Brendan Malone}
\references{
\itemize{
  \item Congalton, R. G. (1991). \href{http://dx.doi.org/10.1016/0034-4257(91)90048-B}{A review of assessing the accuracy of classifications of remotely sensed data}. \emph{Remote Sensing of Environment}, 37, 35–46.
}
}
\examples{
library(ithir)

# Using a pre-constructed confusion matrix
con.mat <- matrix(c(5, 0, 1, 2,
                    0, 15, 0, 5,
                    0, 1, 31, 0,
                    0, 10, 2, 11), nrow = 4, byrow = TRUE)
rownames(con.mat) <- colnames(con.mat) <- c("DE", "VE", "CH", "KU")
goofcat(conf.mat = con.mat, imp = TRUE)

# Using vectors of observed and predicted values
set.seed(123)
observed <- sample(1:5, 1000, replace = TRUE)
set.seed(321)
predicted <- sample(1:5, 1000, replace = TRUE)
goofcat(observed = observed, predicted = predicted)
}
\keyword{methods}
