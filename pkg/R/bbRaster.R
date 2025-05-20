# Purpose        : Returns the bounding box coordinates of a SpatRaster object
# Maintainer     : Brendan Malone (brendan.malone@csiro.au)
# Contributions  : —
# Status         : Working
# Note           : Coordinates are returned as a 4x2 matrix.
# Last Modified  : 2025-05-20

#' @title Bounding box of a raster
#' @description Returns the bounding box of a SpatRaster as a 4x2 matrix of coordinate pairs.
#' @param obj SpatRaster object (from the terra package).
#' @return A 4x2 matrix where each row is a coordinate pair (xmin/ymin, xmin/ymax, xmax/ymin, xmax/ymax).
#' @export

bbRaster <- function(obj) {
  # Get extent from SpatRaster
  ext <- terra::ext(obj)
  
  # Create bounding box matrix (corners: bottom-left, top-left, bottom-right, top-right)
  bbox <- matrix(NA, nrow = 4, ncol = 2)
  bbox[1, ] <- c(ext[1], ext[3])  # xmin, ymin
  bbox[2, ] <- c(ext[1], ext[4])  # xmin, ymax
  bbox[3, ] <- c(ext[2], ext[3])  # xmax, ymin
  bbox[4, ] <- c(ext[2], ext[4])  # xmax, ymax
  
  return(bbox)
}
