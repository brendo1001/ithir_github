\name{precompute_spline_structures}
\docType{methods}
\alias{precompute_spline_structures}
\title{Precompute Spline Matrix Structures for Soil Profile Modeling}
\description{
Builds the matrix structures needed to fit a mass-preserving spline model, based on a given set of soil depth intervals and a smoothing parameter.
These matrices are reused across all soil profiles or raster cells for efficient evaluation.
}
\usage{
precompute_spline_structures(dIn, lam = 0.1)
}
\arguments{
  \item{dIn}{numeric vector of input soil depth boundaries (e.g., \code{c(0,5,15,30,60,100,200)}).}
  \item{lam}{numeric; spline smoothing parameter (\eqn{\lambda}). Smaller values produce smoother splines.}
}
\value{
Returns a \code{list} containing precomputed matrix components used in mass-preserving spline fitting:
\itemize{
  \item{\code{z}}{ - coefficient matrix used to solve for spline coefficients}
  \item{\code{rinv}}{ - inverse of spline smoothness matrix}
  \item{\code{q}}{ - difference matrix for first derivatives}
  \item{\code{u}}{ - upper depths of each interval}
  \item{\code{v}}{ - lower depths of each interval}
  \item{\code{delta}}{ - thickness of each interval}
}
}
\author{Brendan Malone}
\examples{
# Standard SLGA input depths
dIn <- c(0,5,15,30,60,100,200)

# Compute the spline matrices
spline_info <- precompute_spline_structures(dIn = dIn, lam = 0.1)

# View structure
str(spline_info)
}
\seealso{\code{\link{fit_mpspline_optimized}}, \code{\link{ea_rasSp_fast}}}
\keyword{methods}
