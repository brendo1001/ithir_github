#library(aqp)
#library(plyr)
#library(sf)
#library(ithir)

#Fit spline 
#data(oneProfile)
#class(oneProfile)
#sp.fit<- ithir::ea_spline(obj = oneProfile, var.name="C.kg.m3.")

#Using a SoilProfileCollection
## sample profile from Nigeria:
#lon = 3.90; lat = 7.50; id = "ISRIC:NG0017"; FAO1988 = "LXp" 
#top = c(0, 18, 36, 65, 87, 127) 
#bottom = c(18, 36, 65, 87, 127, 181)
#ORCDRC = c(18.4, 4.4, 3.6, 3.6, 3.2, 1.2)
#munsell = c("7.5YR3/2", "7.5YR4/4", "2.5YR5/6", "5YR5/8", "5YR5/4", "10YR7/3")
## prepare a SoilProfileCollection:
#prof1 <- join(data.frame(id, top, bottom, ORCDRC, munsell), 
#         data.frame(id, lon, lat, FAO1988), type='inner')
#aqp::depths(prof1) <- id ~ top + bottom
#aqp::site(prof1) <- ~ lon + lat + FAO1988 

## fit spline:
#ORCDRC.s <- ea_spline(prof1, var.name="ORCDRC")
#str(ORCDRC.s)


# library(ithir)
# library(terra)

# load 
# elevation <- rast(system.file("extdata/edgeGrids_Elevation.tif", package="ithir"))

# simple plot
#plot(elevation, main= "Edgeroi Elevation Map")

# library(ithir)
# library(sf)

## load data
# data(edgeroi_splineCarbon)

## plot the point locations
# spat_edgeroi_splineCarbon<- sf::st_as_sf(x = edgeroi_splineCarbon,coords = c("east", "north"))
# plot(spat_edgeroi_splineCarbon, pch = 9, cex = 0.3)

## Goof

## NOT RUN
# library(ithir)
# library(MASS)
## some data
# data(USYD_soil1)
## fit a linear model
# mod.1 <- lm(CEC ~ clay, data = USYD_soil1 , y = TRUE, x = TRUE)
## Goodness of fit
# goof(observed = mod.1$y, predicted = mod.1$fitted.values, plot.it = TRUE)


# goofcat

## NOT RUN
## Using a pre-constructed confusion matrix
# con.mat <- matrix(c(5, 0, 1, 2, 0, 15, 0, 5, 0, 1, 31, 0, 0, 10, 2, 11),nrow = 4, ncol = 4)
# rownames(con.mat) <- c("DE", "VE", "CH", "KU")
# colnames(con.mat) <- c("DE", "VE", "CH", "KU")
# goofcat(conf.mat = con.mat, imp=TRUE)

## NOT RUN
## Using observations and corresponding predictions
## Using random intgers
# set.seed(123)
# observed<- sample(1:5,1000,replace=TRUE)
# set.seed(321)
# predicted<- sample(1:5,1000,replace=TRUE)
# goofcat(observed = observed, predicted = predicted)

 library(ithir)
 data(homosoil_globeDat)
 str(homosoil_globeDat)
 
 
 ## HV DEM
 # library(ithir)
 # library(terra)
 
 # data(HV_dem)
 # map<- terra::rast(x = HV_dem, type = "xyz")
 # plot(map, main = "Hunter Valley DEM") 
 
 # plot spline
library(ithir)
library(aqp)
 
 
##NOT RUN
data(oneProfile)
str(oneProfile)
##convert to SoilProfileCollection object
aqp::depths(oneProfile)<- Soil.ID ~ Upper.Boundary + Lower.Boundary
##fit spline
eaFit <- ea_spline(oneProfile, var.name="C.kg.m3.",d= t(c(0,5,15,30,60,100,200)),lam = 0.1, vlow=0, show.progress=FALSE )
##do plot
plot_ea_spline(splineOuts=eaFit, d= t(c(0,5,15,30,60,100,200)), maxd=200, type=1, label="carbon density") 


# plot soil profile

# library(ithir)

## NOT RUN
# data(oneProfile)
# str(oneProfile)

## do plot
# plot_soilProfile(data = oneProfile, vals = oneProfile$C.kg.m3., depths = oneProfile[,2:3], label= names(oneProfile)[4])
