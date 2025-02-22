##  AssignR script to place the two Costa Rica birds

##  G Koehler, 2022, 2024  - use new Terra library

###############################################################################################################
library(terra)
library(assignR)
library(maps)
#library(rgdal)#load AssignR library

args = commandArgs(trailingOnly=TRUE)

RTHUfunct <- function(x) {x1.22-20.17}  #transfer fucntion for D in RTHU
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
  
}

modOverlap <- function(pred1, pred2, na.rm = TRUE) {
  
  # version 2.0 (3 Oct 2024)
  
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
#-------------------------------------------------------------
#     Load the RWCIP 2H growing season isoscape and standard errors
#-------------------------------------------------------------

Gisoscape <- c(rast("RCWIP_grid_6_14_2H_GS.tif"), rast("RCWIP_grid_err_6_14_100_2H_GS.tif "))

#-------------------------------------------------------------
#   Or load them from the AssignR database
#-------------------------------------------------------------

#Gisoscape <- getIsoscapes(isoType = "GlobalPrecipGS", timeout = 1200)  #load the gloabal GS isoscape from waterisotopes.org
Hisoscape <- subset(Gisoscape, 1:2) # extract d2H and d2Hse
iso_RTHU = c(Hisoscape[[1]]*1.19 -14.98, Hisoscape[[2]]*1.19) #calibrate RTHU isoscape RWCIP
#plot(iso_RTHU, layer=1)

shapefile <- "Archilochus_colubris.shp"
shapefile2 <- "Archilochus_colubris_winter.shp"
shapefile3 <- "RTHU_summer.shp"
basename(shapefile)
basename(shapefile2)  #strip out the path
basename(shapefile3)
range_t <- vect(shapefile) #load breeding ground range of RTHU
range_w <- vect(shapefile2)
range_s <- vect(shapefile3)

crs(iso_RTHU) <- crs(range_s) #set crs to the same

#iso_RTHU<-mask(iso_RTHU, range_s)  #dont need this because mask function is built into pdfraster()
#ext<-extent(range) ### object to store spatial extent of the range map
#iso_RTHU<-crop(iso_RTHU,ext)
#plot(iso_RTHU, layer=1)
#plot(range, add=T)
#map('world', add=T)

# two unknown CR captured RTHU

id = c("H1", "H2", "H3") #our two birds H1,H2, and those of Hutchinson et al 2010, "H3"
d2H = c(-110, -97, -56)
d2H.sd = c(2,2,2)
un = data.frame(id,d2H, d2H.sd)
custom.palette1 <- colorRampPalette(c("gray92","skyblue", "royalblue", "green", "yellow", "orange","red"),space="rgb") #colorpalete
custom.palette2 <- colorRampPalette(c("gray92","yellow2","mediumseagreen","blue4"),space="rgb") #Colourblind palette 2


asn = pdRaster(iso_RTHU, unknown = un, mask = range_s, genplot = FALSE)  #assign them and plot
#all_asn = qtlRaster(asn, threshold = 0.25, thresholdType = "prob", genplot = FALSE)  #assign most probable 25%

# Plot maps for export

plot(asn, 1, col=custom.palette2(24), xlab="long", ylab="lat")
map('world', add=T)
plot(range_s, add=T)

plot(asn, 2, col=custom.palette2(24), xlab="long", ylab="lat")
map('world', add=T)
plot(range_s, add=T)

plot(asn, 3, col=custom.palette2(24), xlab="long", ylab="lat")
map('world', add=T)
plot(range_s, add=T)

#compare rasters
# result <- all.equal(asn[[1]] - asn[[2]])  #doesnt work


result2 <- layerCor(asn, "pearson", na.rm=TRUE)
onA <- compare_spatrasters(asn[[1]], asn[[3]]) #calculate p-values for cor cooeficcients between assignment surfaces
onB <- compare_spatrasters(asn[[2]], asn[[3]]) #along with Schoeners D value for surface correlation 
onC <- compare_spatrasters(asn[[1]], asn[[2]])

onA2 <- modOverlap(asn[[1]], asn[[3]]) #modOverlap from fuzzysim package - provides quantification of similarity of surfaces
onB2 <- modOverlap(asn[[2]], asn[[3]]) # specifically, Persons correlation, Schoeners D, and Warrens I.
onC2 <- modOverlap(asn[[1]], asn[[2]])

#output to file stats.txt
sink("stats2.txt")
cat("comparison \t corr \t D \t I \t p-value\n")
cat("ON1-ON2 \t", onC$estimate, "\t", onC2$SchoenerD, "\t", onC2$WarrenI, "\t", onC$p.value,"\n")
cat("ON2-ON3 \t", onB$estimate, "\t", onB2$SchoenerD, "\t", onB2$WarrenI, "\t", onB$p.value,"\n")
cat("ON1-ON3 \t", onA$estimate, "\t", onA2$SchoenerD, "\t", onA2$WarrenI, "\t", onA$p.value,"\n")
sink()



if (args[1] == "genplot") {

# Hardcopy

setEPS()
postscript("RTHUCRfig1.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 4, width = 5)
plot(asn, 1, col=custom.palette2(24), xlab="long", ylab="lat")
map('world', add=T)
plot(range_s, add=T)

postscript("RTHUCRfig2.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 4, width = 5)
plot(asn, 2, col=custom.palette2(24), xlab="long", ylab="lat")
map('world', add=T)
plot(range_s, add=T)

postscript("RTHUCRfig3.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 4, width = 5)
plot(asn, 3, col=custom.palette2(24), xlab="long", ylab="lat")
map('world', add=T)
plot(range_s, add=T)

dev.off()
}
