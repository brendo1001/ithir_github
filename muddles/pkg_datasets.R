### Ithir R package
# package data 
# Author: Brendan Malone
# Email: brendan.malone@csiro.au
# Date Created: 16.8.23
# Date Modified: 16.8.23

# CODE PURPOSE
# Getting data organised for the package
# Using exisitng data some modifcations needed given shift to terra and sf packages
##

library(ithir)

# processing edgeTarget_C.rda
library(terra)
target <- rast(system.file("extdata/edgeTarget_C.tif", package="ithir"))
target
bbRaster(target)

## edgeGrids.rda to a spatRaster stack
data("edgeGrids")
edgeGrids
# save to geotiff
n<- nlayers(edgeGrids)
for (i in 1:n){
  r1<- rast(edgeGrids[[i]])
  nm<- names(r1)
  path<- "/home/brendo1001/mystuff/devs/ithir_github/pkg/inst/extdata/"
  terra::writeRaster(x = r1, filename = paste0(path,"edgeGrids_", nm, ".tif"),overwrite=T)}

## edgeroiCovariates to Geotiffs
data(edgeroiCovariates)
path<- "/home/brendo1001/mystuff/devs/ithir_github/pkg/inst/extdata/"
terra::writeRaster(x = elevation, filename = paste0(path,"edgeroiCovariates_elevation.tif"),overwrite=T)
terra::writeRaster(x = landsat_b3, filename = paste0(path,"edgeroiCovariates_landsat_b3.tif"),overwrite=T)
terra::writeRaster(x = landsat_b4, filename = paste0(path,"edgeroiCovariates_landsat_b4.tif"),overwrite=T)
terra::writeRaster(x = radK, filename = paste0(path,"edgeroiCovariates_radK.tif"),overwrite=T)
terra::writeRaster(x = twi, filename = paste0(path,"edgeroiCovariates_twi.tif"),overwrite=T)

## hunterCovariates_sub.rda to a spatRaster stack
data("hunterCovariates_sub")
hunterCovariates_sub
# save to geotiff
n<- nlayers(hunterCovariates_sub)
n
for (i in 1:n){
  r1<- rast(hunterCovariates_sub[[i]])
  nm<- names(r1)
  path<- "/home/brendo1001/mystuff/devs/ithir_github/pkg/inst/extdata/"
  terra::writeRaster(x = r1, filename = paste0(path,"hunterCovariates_sub_", nm, ".tif"),overwrite=T)}

## hunterCovariates_sub.rda to a spatRaster stack
data("hunterCovariates")
hunterCovariates
# save to geotiff
n<- nlayers(hunterCovariates)
n
for (i in 1:n){
  r1<- rast(hunterCovariates_sub[[i]])
  nm<- names(r1)
  path<- "/home/brendo1001/mystuff/devs/ithir_github/pkg/inst/extdata/"
  terra::writeRaster(x = r1, filename = paste0(path,"hunterCovariates_", nm, ".tif"),overwrite=T)}

## edgeroiCovariates to Geotiffs
data("hvGrid25m")
hvGrid25m
path<- "/home/brendo1001/mystuff/devs/ithir_github/pkg/inst/extdata/"
terra::writeRaster(x = hvGrid25m, filename = paste0(path,"hvGrid25m_grid.tif"),overwrite=T)









