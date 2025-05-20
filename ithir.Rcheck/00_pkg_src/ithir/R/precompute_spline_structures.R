# Purpose        : Precompute spline matrix components for mass-preserving spline fitting
# Maintainer     : Brendan Malone (brendan.malone@csiro.au)
# Contributions  : 
# Status         : working
# Note           : Used to accelerate per-profile spline fitting (e.g., over rasters)
# Last Modified  : 2025-05-20

# Description:
# This function prepares the spline system matrices needed for mass-preserving spline interpolation
# based on specified soil layer boundaries (dIn) and a smoothing parameter (lam). These precomputed
# components are reused when fitting multiple profiles (e.g. from raster stacks) to avoid recomputation.

# Inputs:
#   dIn : numeric vector of input soil depth boundaries (must be increasing, e.g. c(0,5,15,...))
#   lam : smoothing parameter lambda controlling spline stiffness (default = 0.1)

# Output:
#   A list containing:
#     - z      : coefficient matrix for spline solution
#     - rinv   : inverse of the smoothness matrix r
#     - q      : first difference matrix
#     - u, v   : upper and lower bounds of each depth layer
#     - delta  : layer thicknesses

precompute_spline_structures <- function(dIn, lam = 0.1) {
  
  # Extract upper and lower bounds of each depth interval
  u <- dIn[-length(dIn)]
  v <- dIn[-1]
  n <- length(u)  # number of layers
  
  # Compute thicknesses of each layer
  delta <- v - u
  
  # Compute differences between subsequent layers for curvature constraint
  del <- c(u[-1], u[n]) - v
  
  nm1 <- n - 1  # one fewer than number of layers
  
  # Construct the r matrix for smoothness penalties
  r <- matrix(0, nm1, nm1)
  diag(r) <- 1
  r[row(r) == col(r) - 1] <- 1  # upper diagonal of 1s
  d2 <- diag(delta[-1], nrow = nm1, ncol = nm1)
  r <- d2 %*% r + t(d2 %*% r)                         # symmetric smoothing matrix
  r <- r + 2 * diag(delta[1:nm1]) + 6 * diag(del[1:nm1])
  
  # Construct the q matrix (difference matrix for first derivatives)
  q <- matrix(0, n, n)
  diag(q) <- -1
  q[row(q) == col(q) - 1] <- 1
  q <- q[1:nm1, , drop = FALSE]
  
  # Invert r (can fail for pathological depth configurations)
  rinv <- try(solve(r), silent = TRUE)
  if (!is.matrix(rinv)) return(NULL)  # return NULL if singular
  
  # Compute the z matrix used in spline solution
  z <- 6 * n * lam * t(q) %*% rinv %*% q + diag(n)
  
  # Return all structures needed for spline fitting
  return(list(
    z = z,
    rinv = rinv,
    q = q,
    u = u,
    v = v,
    delta = delta
  ))
}
