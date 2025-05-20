\name{fuzzyEx}
\docType{methods}
\alias{fuzzyEx}
\title{Derivation of Fuzzy Memberships to Classes with Extragrades}
\description{
This function computes fuzzy membership values for a set of multivariate observations given a set of class centroids, based on the fuzzy k-means with extragrades algorithm (McBratney and de Gruijter, 1992). It is designed to complement outputs from the FuzMe software (Minasny and McBratney, 2002).
}
\usage{
fuzzyEx(data, centroid, cv, expon, alfa)
}
\arguments{
  \item{data}{A \code{data.frame} where the first column is an observation identifier, and the remaining columns contain numeric values for the variables used in the classification.}
  \item{centroid}{A numeric matrix of class centroids (excluding extragrade).}
  \item{cv}{A variance-covariance matrix used to compute Mahalanobis distances between observations and class centroids.}
  \item{expon}{A numeric value specifying the fuzzy exponent (e.g., 1.2 or 2).}
  \item{alfa}{A numeric scalar indicating the membership threshold for the extragrade class, as returned from FuzMe.}
}
\value{
A numeric matrix with one row per observation and \code{n + 2} columns, where \code{n} is the number of classes. The columns include:
\describe{
  \item{Class 1 to Class n}{Estimated fuzzy membership to each class}
  \item{Extragrade}{Estimated membership to the extragrade class}
  \item{MaxClass}{The class (including extragrade) with the highest membership for each observation}
}
}
\note{
Mahalanobis distance is used to evaluate dissimilarity between observations and centroids, incorporating variable covariance structure.
}
\author{Brendan Malone}
\references{
\itemize{
  \item McBratney, A.B., de Gruijter, J.J. (1992). \href{http://dx.doi.org/10.1111/j.1365-2389.1992.tb00127.x}{Continuum approach to soil classification by modified fuzzy k-means with extragrades}. Journal of Soil Science, 43, 159–175.
  \item Minasny, B., McBratney, A.B. (2002). \href{http://sydney.edu.au/agriculture/pal/software/fuzme.shtml}{FuzMe version 3.0}. Australian Centre for Precision Agriculture, The University of Sydney.
}
}
\examples{
## Not run: 
# Example (requires compatible centroid and covariance data)
# memberships <- fuzzyEx(data = input_data, 
#                        centroid = centroid_matrix, 
#                        cv = cov_matrix, 
#                        expon = 2, 
#                        alfa = 0.01)
## End(Not run)
}
\keyword{methods}
