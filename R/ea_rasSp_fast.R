# Purpose        : Apply a mass-preserving spline to a multi-layer SpatRaster stack
# Maintainer     : Brendan Malone (brendan.malone@csiro.au)
# Contributions  : 
# Status         : working
# Note           : Efficient raster implementation using precomputed spline matrices
# Last Modified  : 2025-05-20

# Description:
# This function fits a mass-preserving spline to a raster stack representing soil values 
# at multiple depth intervals. The function uses precomputed spline matrix components 
# and applies the spline cell-wise using terra::app() for speed and scalability.

# Inputs:
#   obj         : a SpatRaster with layers corresponding to input depth intervals
#   lam         : smoothing parameter lambda (default = 0.1)
#   dIn         : numeric vector of input depth boundaries (e.g., c(0,5,15,...))
#   dOut        : numeric vector of output depth intervals (e.g., c(0,30,60))
#   vlow        : minimum value to clip predictions to (e.g., 0)
#   vhigh       : maximum value to clip predictions to (e.g., 100)
#   depth_res   : resolution for spline interpolation along depth (default = 1 cm)

# Output:
#   A SpatRaster with one layer per output interval, containing average values
#   from the fitted spline function over each depth range.

ea_rasSp_fast <- function(obj,
                          lam = 0.1,
                          dIn = c(0,5,15,30,60,100,200),
                          dOut = c(0,30,60),
                          vlow = 0,
                          vhigh = 100,
                          depth_res = 1) {
  
  # Precompute the spline matrix components to be reused across all cells
  spline_info <- precompute_spline_structures(dIn = dIn, lam = lam)
  
  # Stop if the system matrix could not be inverted
  if (is.null(spline_info)) {
    stop("Spline matrix could not be inverted. Check depth intervals or smoothing parameter.")
  }
  
  # Define a wrapper function that will be applied to each raster cell (profile)
  spline_apply <- function(x) {
    fit_mpspline_optimized(
      vals = x,
      spline_info = spline_info,
      dOut = dOut,
      vlow = vlow,
      vhigh = vhigh,
      depth_res = depth_res
    )
  }
  
  # Apply the spline function to each cell in the raster
  out_rast <- terra::app(obj, fun = spline_apply, filename = "")
  
  # Name output layers based on output depth intervals
  names(out_rast) <- paste0(dOut[-length(dOut)], "_", dOut[-1], "cm")
  
  # Return the final SpatRaster with spline-fitted averages
  return(out_rast)
}
