###  Process results into easy graph for viewing
###  Part of the RTHU project
###  G. Koehler 2022, 2023, 2024
###  

######################################################################
#setwd("~/Dropbox/RTHU2020/Data")  #to files data location
#source("rma.R")  # load RMA fuction

RTHU <- read.csv("RTHU_byloc.csv", header=TRUE) #read data 
RTHU_bg <- read.csv("RTHU_boysgirls.csv") # read data for H and O isotopic compositions according to sex
#split dataframe into vectors for boxplot
CR <- na.omit(RTHU$CO)
SK <- na.omit(RTHU$SK)
IL <- na.omit(RTHU$IL)
ON <- na.omit(RTHU$ON)
boysH <- na.omit(RTHU_bg$boysH)
boysO <- na.omit(RTHU_bg$boysO)
girlsH <- na.omit(RTHU_bg$girlsH)
girlsO <- na.omit(RTHU_bg$girlsO)

#plot it
CRmean <- mean(CR)
CRsd <- sd(CR)
SKmean <- mean(SK)
SKsd <- sd(SK)
ILmean <- mean(IL)
ILsd <- sd(IL)
ONmean <- mean(ON)
ONsd <- sd(ON)
boysHmean <- mean(boysH)
girlsHmean <- mean(girlsH)
boysOmean <- mean(boysO)
girlsOmean <- mean(girlsO)
boysH_sd <- sd(boysH)
girlsH_sd <- sd(girlsH)
boysO_sd <- sd(boysO)
girlsO_sd <- sd(girlsO)


cat("SK - mean=", SKmean, "+/-", SKsd, "\n")
cat("IL - mean=", ILmean, "+/-", ILsd, "\n")
cat("ON- me an=", ONmean, "+/-", ONsd, "\n")
cat("CR - mean=", CRmean, "+/-", CRsd, "\n")
cat("male dD =", boysHmean, "+/-", boysH_sd, "\n")
cat("female dD =", girlsHmean, "+/-", girlsH_sd, "\n")
cat("male d18O =", boysOmean, "+/-", boysO_sd, "\n")
cat("female d18O =", girlsOmean, "+/-", girlsO_sd, "\n")
#plot it and create figure.
t.test(boysH, girlsH)
t.test(boysO, girlsO)
t.test(SK, IL)
t.test(SK, ON)
t.test(IL, ON)


setEPS() #output eps file 
postscript("RTHUfig2.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 4, width = 4)
boxplot(SK, IL, ON, CR,
        names = c("SK", "IL", "ON", "CR"),
#        boxwex = 0.2,
        ylab=expression(delta^2~H[f])
)
dev.off() 

setEPS()# output eps file for sex separation
postscript("boysgirls.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 5, width = 5)
par(mfrow = c(1, 2))
boxplot(boysH, girlsH,
        names = c("M", "F"),
        #        boxwex = 0.2,
        ylab=expression(delta^2~H[f])
)

boxplot(boysO, girlsO,
        names = c("M", "F"),
        #        boxwex = 0.2,
        ylab=expression(delta^{18}~O[f])
)
dev.off() 
