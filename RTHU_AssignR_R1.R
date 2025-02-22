###########################################################################
##  AssignR script to place wintering grounds of RTHU

##  G Koehler, 2022, 2023, 2024, 2025

###########################################################################
library(terra, quietly = TRUE)
library(assignR, quietly = TRUE)
library(maps, quietly = TRUE)
library(raster, quietly = TRUE)
# library(fuzzySim, quietly=TRUE)  modOverlap function shamelessly stolen and modified


#library(rgdal)
# args = commandArgs(trailingOnly=TRUE)  # Uncomment if running from command line.

RTHUfunct <- function(x) {x1.17-22.15}  #transfer fucntion for D in RTHU

source("svwodds.R")

compare_spatrasters <- function(raster1, raster2) {
  #compare spatrasters using pearson correlations and Schoeners D
  # Check if the rasters have the same extent and resolution
  if (!compareGeom(raster1, raster2)) {
    stop("The rasters do not have the same extent and/or resolution.")
  }
  
  # Extract cell values
  values1 <- values(raster1, na.rm = TRUE)
  values2 <- values(raster2, na.rm = TRUE)
  
  # Perform paired t-test
  test_result <- t.test(values1, values2, paired = TRUE)
  corr_result <- cor.test(values1, values2, method = "pearson", paired = TRUE)
  
  #Calculate Schoeners D 
  
  D <-  1 - (0.5 * global(abs(raster1 - raster2), 'sum', na.rm=TRUE))
#  I <-  1 - global((sqrt(raster1) - sqrt(raster2))^2)/2), 'sum', na.rm=TRUE)

  
  return(corr=corr_result)
}


# Calculate Schoener's D-metric of spatial similarity between 
# two probability surfaces.  Modified from isocat to work with the Terra library.
schoenersD <- function(rast1, rast2){ 
  # Check if the rasters have the same extent and resolution
  if (!compareGeom(rast1, rast2)) {
    stop("The rasters do not have the same extent and/or resolution.")
  }
  
  
  D <-  1 - (0.5 * global(abs(rast1 - rast2), 'sum', na.rm=TRUE))
  return(D$sum)
  
} #end function

modOverlap <- function(pred1, pred2, na.rm = TRUE) {
  
  # version 2.0 (3 Oct 2024) shamelessly stolen from fuzzysim package.
  
  if (inherits(pred1, "SpatRaster"))
    pred1 <- terra::values(pred1, mat = FALSE, dataframe = FALSE)
  if (inherits(pred2, "SpatRaster"))
    pred2 <- terra::values(pred2, mat = FALSE, dataframe = FALSE)
  
  pred1 <- unlist(pred1)
  pred2 <- unlist(pred2)  # can't remember what this was for...
  
  stopifnot(length(pred1) == length(pred2),
            #pred1 >= 0 && pred1 <=1,
            #pred2 >= 0 && pred2 <=1
            min(c(pred1, pred2), na.rm = TRUE) >= 0,
            max(c(pred1, pred2), na.rm = TRUE) <= 1
  )
  
  p1 <- pred1 / sum(pred1, na.rm = na.rm)
  p2 <- pred2 / sum(pred2, na.rm = na.rm)
  
  SchoenerD <- 1 - 0.5 * sum(abs(p1 - p2), na.rm = na.rm)
  HellingerDist <- sqrt(sum((sqrt(p1) - sqrt(p2))^2, na.rm = na.rm))
  WarrenI <- 1 - ((HellingerDist^2) / 2)
  
  list(SchoenerD = SchoenerD,
       WarrenI = WarrenI,
       HellingerDist = HellingerDist
  )
}  # end modOverlap function

calibRTHU <- function(isodef) {
  #read in calibratin values from file
  RTHU.calib <- read.csv("RTHU-calib.csv")
  
  
  # Longitude and latitude values
  long <- RTHU.calib$Long
  lat <- RTHU.calib$Lat
  longlat <- cbind(long, lat)
  
  # CRS
  crspoints <- "+proj=longlat +datum=WGS84"
  
  # Attributes for points
  d <- data.frame( d2H = RTHU.calib$dDf, se = RTHU.calib$SD)
  
  # SpatVector object for calraster function
  pts <- vect(longlat, atts = d, crs = crspoints)
  
  #get global GS isoscape 
  
  if (isodef == 1) { # get RWCIP GS isoscape
   
    Hisoscape <- c(rast("RCWIP_grid_6_14_2H_GS.tif"), rast("RCWIP_grid_err_6_14_100_2H_GS.tif "))
    RTHU.isoscape = calRaster(known = pts, isoscape = Hisoscape)
    return(RTHU.isoscape )
  } #end if
  else {
    Gisoscape <- getIsoscapes(isoType = "GlobalPrecipGS", timeout = 1200)  #load the gloabal GS isoscape from waterisotopes.org
    Hisoscape <- subset(Gisoscape, 1:2) #extract hydrogen isosccape and se
    RTHU.isoscape = calRaster(known = pts, isoscape = Hisoscape)
    return(RTHU.isoscape)
  }  #end else
}

isoassign <- function(isodata, H2O.isoscape, range) {
  
  #  library(terra)  #load necessary libraries
  #  library(rgdal)
  
  
  crs(H2O.isoscape) <- crs(range)
  SD <-1.4 
  
  # isoscape<- calib(H2O.isoscape)  # transform H2O isoscape into tissue isoscape 
  
  ext<-ext(range) ### object to store spatial extent of the range map
  isoscape<-crop(H2O.isoscape,ext)  #crop isoscape to range extent
  isoscape<-mask(isoscape, range) # fill in map with NA outside range
  
  prior<-1  ### PLACEHOLDER for prior probabilities, current value is a 'uniform prior'
  
  
  ###### establish a list to hold the results of subsequent loops, objects in each list stores results separately so that each can be retrieved/examined individually
  r <- list()
  
  ######## Apply Normal Probability Density Function to Assess Likelihoods and Normalize to the sum #########
  for (i in 1:length(isodata)){
    r[i] <- (1/(sqrt(2*pi)*SD)*exp(-(1/(2*(SD)^2))*(isodata[i]-isoscape)^2)*prior)
    #    r[i] <- r[[i]] / global(r, fun ="sum")  ### normalizes data so whole map sums to one
  }
  
  
  
  
  ## to generate depiction of origins across all samples, first create a raster stack
  s <- rast(r)
  
  
  c.prob<-(app(s, fun=function(x) (cumprod(x[i]))))
  #   c.prob<- c.prob/ global(c.prob, fun ="sum")  ### normalize map values to the sum
  return(c.prob) #return cumulative probabilities object
  
  #####   Plot the cumulative probabilities  
  
  
}




#-------------------------------------------------------------
# isoscape definition for calibRTHU() function.  1=RCWIP GS hydrogen isoscape
#                                                2=Bowen MA isoscape
#                                                0=BowenGS isoscape  
isodef = 0

#-------------------------------------------------------------
#     Load the RWCIP 2H growing season isoscape and standard errors
#-------------------------------------------------------------

Hisoscape <- c(rast("RCWIP_grid_6_14_2H_GS.tif"), rast("RCWIP_grid_err_6_14_100_2H_GS.tif ")) # -> RWCIP (Terzer)
#Hisoscape <- c(rast("d2h_GS.tif"), rast("d2h_se_GS.tif"))     # -> Bowen GS
#-------------------------------------------------------------
#   Or load them from the AssignR database
#-------------------------------------------------------------

#Gisoscape <- getIsoscapes(isoType = "GlobalPrecipGS", timeout = 1200)  #load the gloabal GS isoscape from waterisotopes.org
#Hisoscape <- subset(Gisoscape, 1:2) # extract d2H and d2Hse

#-------------------------------------------------------------
# calibrate GS isoscape to RTHU isoscape

#iso_RTHU = c(Hisoscape[[1]]*1.17 -22.15, Hisoscape[[2]]*1.17) #calibrate RTHU isoscape
iso_RTHU = c(Hisoscape[[1]]*1.19 -14.98, Hisoscape[[2]]*1.19) #calibrate RTHU isoscape RWCIP
RTHUiso <- calibRTHU(isodef=isodef)  #use assignR calibration function - works but not RMA
RTHUiso2 <- subset(iso_RTHU, 1) # get hydrogen surface for isoassign functions


#plot(iso_RTHU, layer=1)

#range <- readOGR(dsn="./", layer="Archilochus_colubris_winter") #load Wintering ground range of RTHU

shapefile <- "Archilochus_colubris_winter.shp"
basename(shapefile)
range <- vect(shapefile) #load Wintering ground range of RTHU



crs(iso_RTHU) <- crs(range) #set crs to the same
#iso_RTHU<-mask(iso_RTHU, range) #mask the isoscape - dont need to do this as masked by pdraster
# rext<-ext(range) ### object to store spatial extent of the range map
# iso_RTHU<-crop(iso_RTHU,rext) 
# plot(iso_RTHU, 1)
# #plot(range, add=T)
# #map('world', add=T)

# dD values and standard deviations of RTHU for SK
RTHU_SK <- read.csv("RTHU_SK.csv")
#asn_SK2 = pdRaster(RTHUiso, RTHU_SK, mask = range, genplot = FALSE)  #assign them but dont plot
asn_SK = pdRaster(iso_RTHU, RTHU_SK, mask = range, genplot = FALSE)  #assign them but dont plot
SK_stack = stack(asn_SK)
SK_odds = rast(svwodds(SK_stack, 0.333333))
#asn_SK3 <- isoassign(RTHU_SK$dD_SK, RTHUiso2, range)
all_birds_SK <- unionP(asn_SK) #get union probability of all birds in SK
#all_birds_SK2 <- unionP(asn_SK2)

#Ontario

RTHU_ON <- read.csv("RTHU_ON.csv")
#asn_ON = pdRaster(RTHUiso, RTHU_ON, mask = range, genplot = FALSE)  #assign them but dont plot
asn_ON = pdRaster(iso_RTHU, RTHU_ON, mask = range, genplot = FALSE) 
ON_stack = stack(asn_ON)
ON_odds = rast(svwodds(ON_stack, 0.33333))
all_birds_ON <- unionP(asn_ON) #get union probability of all birds in SK

#Illinois
RTHU_IL <- read.csv("RTHU_IL.csv")
#asn_IL = pdRaster(RTHUiso, RTHU_IL, mask = range, genplot = FALSE)  #assign them but dont plot
asn_IL = pdRaster(iso_RTHU, RTHU_IL, mask = range, genplot = FALSE)  #assign them but dont plot
IL_stack = stack(asn_IL)
IL_odds = rast(svwodds(IL_stack, 0.3333333))
all_birds_IL <- unionP(asn_ON) #get union probability of all birds in SK


#all_asn = qtlRaster(all_birds, threshold = 0.67, thresholdType = "prob", genplot = FALSE) #get 2:1 odds ratio of all birds. 

#compare surfaces
#result2 <- layerCor(asn, "pearson", na.rm=TRUE) #get correlation coefficients between layers.
skon <- compare_spatrasters(all_birds_SK, all_birds_ON) #calculate p-values for cor cooeficcients between assignment surfaces
ilon <- compare_spatrasters(all_birds_IL, all_birds_ON) #along with Schoeners D value for surface correlation 
skil <- compare_spatrasters(all_birds_SK, all_birds_IL)

skon2 <- modOverlap(all_birds_SK, all_birds_ON) #modOverlap from fuzzysim package - provides quantification of similarity of surfaces
ilon2 <- modOverlap(all_birds_IL,all_birds_ON) # specifically, Persons correlation, Schoeners D, and Warrens I.
skil2 <- modOverlap(all_birds_SK, all_birds_IL)

#output to file stats.txt
sink("stats2.txt")
cat("comparison \t corr \t D \t I \t p-value\n")
cat("SK-IL \t", skil$estimate, "\t", skil2$SchoenerD, "\t", skil2$WarrenI, "\t", skil$p.value,"\n")
cat("SK-ON \t", skon$estimate, "\t", skon2$SchoenerD, "\t", skon2$WarrenI, "\t", skon$p.value,"\n")
cat("IL-ON \t", ilon$estimate, "\t", ilon2$SchoenerD, "\t", ilon2$WarrenI, "\t", ilon$p.value,"\n")
sink()

#output 
custom.palette <- colorRampPalette(c("white","gray92","yellow2","mediumseagreen","skyblue","blue4"),space="rgb") ### colourblind palette
custom.palette1 <- colorRampPalette(c("gray92","skyblue", "blue", "green", "yellow", "orange","red"),space="rgb") # I like this colorpalette
custom.palette2 <- colorRampPalette(c("gray92","yellow2","mediumseagreen","blue4"),space="rgb") #Colourblind palette 2

plot(SK_odds, col=custom.palette2(24), main = "SK", xlab="long(deg. W)", ylab="lat(deg)")  #plot SK birds
map('world', add=T)
plot(range, add=T)

plot(IL_odds, col=custom.palette2(24), main = "IL", xlab="long(deg. W)", ylab="lat(deg)")  #plot IL birds
map('world', add=T)
plot(range, add=T)

plot(ON_odds, col=custom.palette2(24), main = "ON", xlab="long(deg. W)", ylab="lat(deg)")  #plot ON birds
map('world', add=T)
plot(range, add=T)




if (args[1] == "genplot") { # Hardcopy
  
  
  setEPS()
  postscript("RTHU_SK.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 5, width = 5)
  plot(all_birds_SK, col=custom.palette2(24), main = "SK", xlab="long", ylab="lat")  #plot SK birds
  map('world', add=T)
  plot(range, add=T)
  
  postscript("RTHU_IL.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 5, width = 5)
  plot(all_birds_IL, col=custom.palette2(24), main = "IL", xlab="long", ylab="lat")  #plot IL birds
  map('world', add=T)
  plot(range, add=T)
  
  postscript("RTHU_ON.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 5, width = 5)
  plot(all_birds_ON, col=custom.palette2(24), main = "ON", xlab="long", ylab="lat")  #plot ON birds
  map('world', add=T)
  plot(range, add=T)
  
  postscript("RTHU_SK_odds.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 5, width = 5)
  plot(SK_odds, col=custom.palette2(24), main = "SK", xlab="long", ylab="lat")  #plot SK birds
  map('world', add=T)
  plot(range, add=T)
  
  postscript("RTHU_IL_odds.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 5, width = 5)
  plot(IL_odds, col=custom.palette2(24), main = "IL", xlab="long", ylab="lat")  #plot IL birds
  map('world', add=T)
  plot(range, add=T)
  
  postscript("RTHU_ON_odds.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 5, width = 5)
  plot(ON_odds, col=custom.palette2(24), main = "ON", xlab="long", ylab="lat")  #plot ON birds
  map('world', add=T)
  plot(range, add=T)
  
  
  dev.off()
}

