pkgname <- "ithir"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('ithir')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("HV100")
### * HV100

flush(stderr()); flush(stdout())

### Name: HV100
### Title: Soil Point Data from the Hunter Valley, NSW, Australia
### Aliases: HV100
### Keywords: datasets

### ** Examples

# Load and inspect the dataset
data(HV100)

# Basic structure and summary
str(HV100)
summary(HV100)



cleanEx()
nameEx("HV_dem")
### * HV_dem

flush(stderr()); flush(stdout())

### Name: HV_dem
### Title: Digital Elevation Model of the Hunter Valley, NSW
### Aliases: HV_dem
### Keywords: datasets

### ** Examples

library(terra)
data(HV_dem)

# Convert to raster and plot
dem_rast <- terra::rast(x = HV_dem, type = "xyz")
plot(dem_rast, main = "Hunter Valley DEM")



cleanEx()
nameEx("HV_subsoilpH")
### * HV_subsoilpH

flush(stderr()); flush(stdout())

### Name: HV_subsoilpH
### Title: Hunter Valley Subsoil pH Data with Environmental Covariates
### Aliases: HV_subsoilpH
### Keywords: datasets

### ** Examples

library(ithir)
data(HV_subsoilpH)
summary(HV_subsoilpH)



cleanEx()
nameEx("USYD_dIndex")
### * USYD_dIndex

flush(stderr()); flush(stdout())

### Name: USYD_dIndex
### Title: Soil Drainage Index Observations and Predictions
### Aliases: USYD_dIndex
### Keywords: datasets

### ** Examples

data(USYD_dIndex)

# View summary statistics
summary(USYD_dIndex)

# Plot observed vs predicted
with(USYD_dIndex, plot(DI_observed, DI_predicted,
     main = "Observed vs Predicted Drainage Index",
     xlab = "Observed", ylab = "Predicted"))



cleanEx()
nameEx("USYD_soil1")
### * USYD_soil1

flush(stderr()); flush(stdout())

### Name: USYD_soil1
### Title: Random Selection of Soil Profile Data from New South Wales
### Aliases: USYD_soil1
### Keywords: datasets

### ** Examples

data(USYD_soil1)

# Overview
summary(USYD_soil1)

# Show profile IDs
table(USYD_soil1$id)



cleanEx()
nameEx("bbRaster")
### * bbRaster

flush(stderr()); flush(stdout())

### Name: bbRaster
### Title: Get Bounding Box Coordinates from a 'SpatRaster'
### Aliases: bbRaster
### Keywords: methods

### ** Examples

library(terra)

# Example raster
target <- rast(system.file("extdata/edgeTarget_C.tif", package = "ithir"))

# Get bounding box coordinates
bbRaster(target)



cleanEx()
nameEx("ea_rasSp_fast")
### * ea_rasSp_fast

flush(stderr()); flush(stdout())

### Name: ea_rasSp_fast
### Title: Fast Mass-Preserving Spline on Raster Soil Data
### Aliases: ea_rasSp_fast
### Keywords: methods

### ** Examples

# Not run on CRAN due to file size and download time
## Not run: 
##D library(terra)
##D 
##D # Define SLGA V2 clay layer URLs (0–200 cm depth range)
##D clay_urls <- c(
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif',
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_005_015_EV_N_P_AU_TRN_N_20210902.tif',
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_015_030_EV_N_P_AU_TRN_N_20210902.tif',
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_030_060_EV_N_P_AU_TRN_N_20210902.tif',
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_060_100_EV_N_P_AU_TRN_N_20210902.tif',
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_100_200_EV_N_P_AU_TRN_N_20210902.tif'
##D )
##D 
##D # Load and crop a small extent near Canberra
##D clay_stack <- rast(clay_urls)
##D aoi <- ext(149.00, 149.10, -36.00, -35.90)
##D clay_crop <- crop(clay_stack, aoi)
##D 
##D # Fit spline and generate interpolated output
##D out <- ea_rasSp_fast(
##D   obj = clay_crop,
##D   lam = 0.1,
##D   dIn = c(0, 5, 15, 30, 60, 100, 200),
##D   dOut = c(0, 30, 60),
##D   depth_res = 2
##D )
##D 
##D # Plot the result
##D plot(out)
## End(Not run)



cleanEx()
nameEx("ea_spline")
### * ea_spline

flush(stderr()); flush(stdout())

### Name: ea_spline
### Title: Fit a Mass-Preserving Spline to Soil Profile Data
### Aliases: ea_spline
### Keywords: methods

### ** Examples

# Example using a simple data.frame
data(oneProfile)
str(oneProfile)
sp_fit <- ea_spline(obj = oneProfile, var.name = "C.kg.m3.")

# Example using a SoilProfileCollection from the aqp package
# library(aqp)
# library(plyr)
# lon <- 3.90; lat <- 7.50; id <- "ISRIC:NG0017"
# top <- c(0, 18, 36, 65, 87, 127)
# bottom <- c(18, 36, 65, 87, 127, 181)
# ORCDRC <- c(18.4, 4.4, 3.6, 3.6, 3.2, 1.2)
# munsell <- c("7.5YR3/2", "7.5YR4/4", "2.5YR5/6", "5YR5/8", "5YR5/4", "10YR7/3")
# prof1 <- join(data.frame(id, top, bottom, ORCDRC, munsell),
#               data.frame(id, lon, lat), type = 'inner')
# depths(prof1) <- id ~ top + bottom
# site(prof1) <- ~ lon + lat
# ORCDRC.s <- ea_spline(prof1, var.name = "ORCDRC")
# str(ORCDRC.s)



cleanEx()
nameEx("edgeLandClass")
### * edgeLandClass

flush(stderr()); flush(stdout())

### Name: edgeLandClass
### Title: Land Classification Points from the Edgeroi District, NSW
### Aliases: edgeLandClass
### Keywords: datasets

### ** Examples

data(edgeLandClass)

# View land class summary
summary(edgeLandClass$LandClass)

# Plot the locations by class
plot(edgeLandClass$x, edgeLandClass$y, col = edgeLandClass$LandClass,
     pch = 20, main = "Edgeroi Land Classification Points")



cleanEx()
nameEx("edgeTarget_C")
### * edgeTarget_C

flush(stderr()); flush(stdout())

### Name: edgeTarget_C
### Title: 1 km Resolution Soil Carbon Stock Map (Subset) – Edgeroi
###   District, NSW
### Aliases: edgeTarget_C edgeTarget_C.tif
### Keywords: datasets

### ** Examples

library(terra)

# Load the raster from the package
soc_path <- system.file("extdata/edgeTarget_C.tif", package = "ithir")
soc_raster <- rast(soc_path)

# Plot the raster
plot(soc_raster, main = "Edgeroi SOC Stock (0–30 cm)", col = rev(terrain.colors(20)))



cleanEx()
nameEx("edgeroiCovariates")
### * edgeroiCovariates

flush(stderr()); flush(stdout())

### Name: edgeroiCovariates
### Title: Environmental Covariate Rasters for the Full Edgeroi District,
###   NSW
### Aliases: edgeroiCovariates edgeroiCovariates_elevation.tif
###   edgeroiCovariates_twi.tif edgeroiCovariates_radK.tif
###   edgeroiCovariates_landsat_b3.tif edgeroiCovariates_landsat_b4.tif
### Keywords: datasets

### ** Examples

library(terra)

# Load elevation layer
elev_path <- system.file("extdata/edgeroiCovariates_elevation.tif", package = "ithir")
elev_rast <- rast(elev_path)
plot(elev_rast, main = "Edgeroi Elevation Map")



cleanEx()
nameEx("edgeroi_covariates_subset")
### * edgeroi_covariates_subset

flush(stderr()); flush(stdout())

### Name: edgeroi_covariates_subset
### Title: Selected Subset of Environmental Covariates for the Edgeroi
###   District, NSW
### Aliases: edgeroi_covariates_subset edgeGrids_Doserate
###   edgeGrids_Elevation edgeGrids_Panchromat edgeGrids_Slope
###   edgeGrids_TWI
### Keywords: datasets

### ** Examples

library(ithir)
library(terra)

# Load and plot the elevation raster
elevation <- rast(system.file("extdata/edgeGrids_Elevation.tif", package = "ithir"))
plot(elevation, main = "Edgeroi Elevation Map")



cleanEx()
nameEx("edgeroi_splineCarbon")
### * edgeroi_splineCarbon

flush(stderr()); flush(stdout())

### Name: edgeroi_splineCarbon
### Title: Harmonised Soil Carbon Density Data from the Edgeroi District,
###   NSW
### Aliases: edgeroi_splineCarbon
### Keywords: datasets

### ** Examples

library(ithir)
library(sf)

# Load the dataset
data(edgeroi_splineCarbon)

# Convert to sf object and plot
spat_edgeroi <- st_as_sf(edgeroi_splineCarbon, coords = c("east", "north"), crs = 32755)
plot(spat_edgeroi, pch = 16, cex = 0.5, main = "Edgeroi Soil Carbon Point Locations")



cleanEx()
nameEx("fit_mpspline_optimized")
### * fit_mpspline_optimized

flush(stderr()); flush(stdout())

### Name: fit_mpspline_optimized
### Title: Fit Mass-Preserving Spline to a Single Soil Profile
### Aliases: fit_mpspline_optimized
### Keywords: methods

### ** Examples

# Not run on CRAN due to external raster data size
## Not run: 
##D library(terra)
##D 
##D # Define SLGA clay raster URLs
##D clay_urls <- c(
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif',
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_005_015_EV_N_P_AU_TRN_N_20210902.tif',
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_015_030_EV_N_P_AU_TRN_N_20210902.tif',
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_030_060_EV_N_P_AU_TRN_N_20210902.tif',
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_060_100_EV_N_P_AU_TRN_N_20210902.tif',
##D   '/vsicurl/https://esoil.io/TERNLandscapes/Public/Products/TERN/SLGA/CLY/CLY_100_200_EV_N_P_AU_TRN_N_20210902.tif'
##D )
##D 
##D # Load and crop raster stack
##D clay_stack <- rast(clay_urls)
##D aoi <- ext(149.00, 149.10, -36.00, -35.90)
##D clay_crop <- crop(clay_stack, aoi)
##D 
##D # Extract a single profile (pixel)
##D vals <- terra::extract(clay_crop, cbind(149.05, -35.95))[1, -1]
##D 
##D # Precompute spline structures
##D dIn <- c(0, 5, 15, 30, 60, 100, 200)
##D spline_info <- precompute_spline_structures(dIn, lam = 0.1)
##D 
##D # Fit spline to single profile
##D fit <- fit_mpspline_optimized(
##D   vals = vals,
##D   spline_info = spline_info,
##D   dOut = c(0, 30, 60),
##D   vlow = 0,
##D   vhigh = 100,
##D   depth_res = 1
##D )
##D 
##D fit
## End(Not run)



cleanEx()
nameEx("fuzzyEx")
### * fuzzyEx

flush(stderr()); flush(stdout())

### Name: fuzzyEx
### Title: Derivation of Fuzzy Memberships to Classes with Extragrades
### Aliases: fuzzyEx
### Keywords: methods

### ** Examples

## Not run: 
# Example (requires compatible centroid and covariance data)
# memberships <- fuzzyEx(data = input_data, 
#                        centroid = centroid_matrix, 
#                        cv = cov_matrix, 
#                        expon = 2, 
#                        alfa = 0.01)
## End(Not run)



cleanEx()
nameEx("goof")
### * goof

flush(stderr()); flush(stdout())

### Name: goof
### Title: Goodness of Fit Measures
### Aliases: goof
### Keywords: methods

### ** Examples

library(ithir)
library(MASS)

# Load sample soil data
data(USYD_soil1)

# Fit a linear model
mod.1 <- lm(CEC ~ clay, data = USYD_soil1, y = TRUE, x = TRUE)

# Calculate goodness-of-fit statistics and plot
goof(observed = mod.1$y, predicted = mod.1$fitted.values, plot.it = TRUE)



cleanEx()
nameEx("goofcat")
### * goofcat

flush(stderr()); flush(stdout())

### Name: goofcat
### Title: Goodness of Fit Measures for Categorical Models
### Aliases: goofcat
### Keywords: methods

### ** Examples

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



cleanEx()
nameEx("homosoil_globeDat")
### * homosoil_globeDat

flush(stderr()); flush(stdout())

### Name: homosoil_globeDat
### Title: Global Environmental Data
### Aliases: homosoil_globeDat
### Keywords: datasets

### ** Examples

library(ithir)
data(homosoil_globeDat)

# Structure of the dataset
str(homosoil_globeDat)

# Summary of the first few climatic variables
summary(homosoil_globeDat[, 1:6])



cleanEx()
nameEx("hunterCovariates")
### * hunterCovariates

flush(stderr()); flush(stdout())

### Name: hunterCovariates
### Title: Environmental Covariate Rasters for the Lower Hunter Valley, NSW
### Aliases: hunterCovariates hunterCovariates_A_AACN
###   hunterCovariates_A_elevation hunterCovariates_A_Hillshading
###   hunterCovariates_A_light_insolation hunterCovariates_A_MRVBF
###   hunterCovariates_A_Slope hunterCovariates_A_TRI
###   hunterCovariates_A_TWI
### Keywords: datasets

### ** Examples

library(ithir)
library(terra)

# Load and plot the slope raster layer
slope <- rast(system.file("extdata/hunterCovariates_A_Slope.tif", package = "ithir"))
plot(slope, main = "Hunter Valley - Slope")



cleanEx()
nameEx("hunterCovariates_sub")
### * hunterCovariates_sub

flush(stderr()); flush(stdout())

### Name: hunterCovariates_sub
### Title: Environmental Covariate Rasters for a Subset of the Lower Hunter
###   Valley, NSW
### Aliases: hunterCovariates_sub hunterCovariates_sub_Elevation
###   hunterCovariates_sub_Slope hunterCovariates_sub_TWI
###   hunterCovariates_sub_MRVBF hunterCovariates_sub_Hillshading
###   hunterCovariates_sub_Light_insolation
###   hunterCovariates_sub_Landsat_Band1 hunterCovariates_sub_NDVI
###   hunterCovariates_sub_AACN
###   hunterCovariates_sub_Terrain_Ruggedness_Index
###   hunterCovariates_sub_Mid_Slope_Position
### Keywords: datasets

### ** Examples

library(ithir)
library(terra)

# Load and plot the elevation raster
elevation <- rast(system.file("extdata/hunterCovariates_sub_Elevation.tif", package = "ithir"))
plot(elevation, main = "Hunter Valley Subset - Elevation")



cleanEx()
nameEx("hvGrid25m")
### * hvGrid25m

flush(stderr()); flush(stdout())

### Name: hvGrid25m
### Title: Raster Grid of the Lower Hunter Valley, NSW, Australia
### Aliases: hvGrid25m
### Keywords: datasets

### ** Examples

library(ithir)
library(terra)

# Load and plot the raster grid
hv.grid <- rast(system.file("extdata/hvGrid25m_grid.tif", package = "ithir"))
plot(hv.grid, main = "Hunter Valley 25m Grid Index")



cleanEx()
nameEx("hvPoints250")
### * hvPoints250

flush(stderr()); flush(stdout())

### Name: hvPoints250
### Title: Hunter Valley Soil Observation Points (n = 250)
### Aliases: hvPoints250
### Keywords: datasets

### ** Examples

data(hvPoints250)

# Summary
summary(hvPoints250)

# Simple plot
plot(hvPoints250, pch = 20, col = "darkgreen", main = "Hunter Valley Sample Points")



cleanEx()
nameEx("hvTerronDat")
### * hvTerronDat

flush(stderr()); flush(stdout())

### Name: hvTerronDat
### Title: Soil Point Data with Terron Class Labels from the Hunter Valley,
###   NSW, Australia
### Aliases: hvTerronDat
### Keywords: datasets

### ** Examples

library(ithir)
data(hvTerronDat)
head(hvTerronDat)
table(hvTerronDat$terron_class)



cleanEx()
nameEx("oneProfile")
### * oneProfile

flush(stderr()); flush(stdout())

### Name: oneProfile
### Title: Example Soil Profile for Soil Carbon Density
### Aliases: oneProfile
### Keywords: datasets

### ** Examples

data(oneProfile)

# Show profile structure
str(oneProfile)

# Fit spline (if using ithir)
# result <- ea_spline(oneProfile, var.name = "carbon_density")
# result$harmonised



cleanEx()
nameEx("plot_ea_spline")
### * plot_ea_spline

flush(stderr()); flush(stdout())

### Name: plot_ea_spline
### Title: Plot Soil Profile Outputs from 'ea_spline'
### Aliases: plot_ea_spline
### Keywords: methods

### ** Examples

library(ithir)
library(aqp)

# Load example data and convert to SoilProfileCollection
data(oneProfile)
depths(oneProfile) <- Soil.ID ~ Upper.Boundary + Lower.Boundary

# Fit spline to the profile
eaFit <- ea_spline(oneProfile, var.name = "C.kg.m3.", 
                   d = t(c(0, 5, 15, 30, 60, 100, 200)), 
                   lam = 0.1, vlow = 0, show.progress = FALSE)

# Plot the fitted spline
plot_ea_spline(splineOuts = eaFit, 
               d = t(c(0, 5, 15, 30, 60, 100, 200)), 
               maxd = 200, 
               type = 1, 
               label = "Carbon Density")



cleanEx()
nameEx("plot_soilProfile")
### * plot_soilProfile

flush(stderr()); flush(stdout())

### Name: plot_soilProfile
### Title: Plot Soil Profile Data
### Aliases: plot_soilProfile
### Keywords: methods

### ** Examples

library(ithir)

# Load example profile data
data(oneProfile)

# Visualize the profile
plot_soilProfile(
  data = oneProfile,
  vals = oneProfile$C.kg.m3.,
  depths = oneProfile[, 2:3],
  label = names(oneProfile)[4]
)



cleanEx()
nameEx("precompute_spline_structures")
### * precompute_spline_structures

flush(stderr()); flush(stdout())

### Name: precompute_spline_structures
### Title: Precompute Spline Matrix Structures for Soil Profile Modeling
### Aliases: precompute_spline_structures
### Keywords: methods

### ** Examples

# Standard SLGA input depths
dIn <- c(0,5,15,30,60,100,200)

# Compute the spline matrices
spline_info <- precompute_spline_structures(dIn = dIn, lam = 0.1)

# View structure
str(spline_info)



cleanEx()
nameEx("topo_dem")
### * topo_dem

flush(stderr()); flush(stdout())

### Name: topo_dem
### Title: Example Digital Elevation Model as a Matrix
### Aliases: topo_dem
### Keywords: datasets

### ** Examples

data(topo_dem)

# Basic inspection
str(topo_dem)
image(topo_dem, main = "Synthetic DEM", col = terrain.colors(20))



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
