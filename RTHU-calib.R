########################################################
###  Routine to calculate RMA regression cooeficients ##
###  for the RTHU Project.  Can easily be adapted ######
###  to other projects.                             ####
###  G. Koehler, 2022                               ####
########################################################

library(lmodel2)
setwd("~/Dropbox/RTHU2020/Data")
source("rma.R")  #Load the RMA function from source
printf <- function(...) invisible(print(sprintf(...)))
compare_models <- function(model1, model2) {
  # Check if both inputs are linear models
  if (!inherits(model1, "lm") | !inherits(model2, "lm")) {
    stop("Both inputs must be linear model objects.")
  }
  
  # Compare the two models using ANOVA
  comparison <- anova(model1, model2)
  
  # Return the comparison result
  return(comparison)
}


###  Load RTHU-precip datasets #######


RTHUcalib <- read.csv("~/Dropbox/RTHU2020/Data/RTHU-calib.csv")  #load 2004 and 2007 data 
zenal2018data <- read.csv("~/Dropbox/RTHU2020/Data/zenal2018data.csv")  # read in ZenzaL 2018 data
#RTHU2004 <- read.csv("~/Dropbox/RTHU2020/Data/RTHUcalibration2004.csv") # read in 2004 data

#RTHUcal <- lmodel2(RTHU2007$dDp ~ RTHU2007$dDf)
#RTHU2004calMA <- lmodel2(RTHU2004$Feather.Mean ~ RTHU2004$Annual.RCWIP)  # lmodel2 regression
RTHU_MA_RCWIP <- rma(RTHUcalib$Annual.RCWIP, RTHUcalib$dDf)  #RMA regression vs RCWIP MA data
RTHU_MA_RCWIP2 = rma(RTHUcalib$Annual..RCWIP2., RTHUcalib$dDf) #RMA vs RCWIP2 MA values 
RTHU_GS = rma(RTHUcalib$Grwowing.Bowen, RTHUcalib$dDf) #RMA vs Growing season Values
RTHU_GS_RCWIP = rma(RTHUcalib$GS.RCWIP, RTHUcalib$dDf) #$RMA vs Growing season RCWIP
ZenzelRMA = rma(zenal2018data$Precip, zenal2018data$Tissue)

plot(RTHUcalib$Annual.RCWIP, RTHUcalib$dDf, col=RTHUcalib$dataset, pch = 16, cex = 2, main="MA RCWIP", ylab=expression(delta^2~H[f]), xlab = expression (delta^2~H[p]), cex.lab=2)
abline(RTHU_MA_RCWIP$rma_coefs)
legend("topleft", legend=c("2004", "2007"), fill = c("blue", "yellow"), title = "RTHU dataset")

plot(RTHUcalib$Annual..RCWIP2., RTHUcalib$dDf, col=RTHUcalib$dataset, pch = 16, cex = 2, main="MA RCWIP2", ylab=expression(delta^2~H[f]), xlab = expression (delta^2~H[p]),cex.lab=2)
abline(RTHU_MA_RCWIP2$rma_coefs)
legend("topleft", legend=c("2004", "2007"), fill = c("blue", "yellow"), title = "RTHU dataset")

plot(RTHUcalib$Grwowing.Bowen, RTHUcalib$dDf, col=RTHUcalib$dataset, pch = 16, cex = 2, main="Growing Season Bowen", ylab=expression(delta^2~H[f]), xlab = expression (delta^2~H[p]),cex.lab=2)
abline(RTHU_GS$rma_coefs)
legend("topleft", legend=c("2004", "2007"), fill = c("blue", "yellow"), title = "RTHU dataset")

plot(RTHUcalib$GS.RCWIP, RTHUcalib$dDf, col=RTHUcalib$dataset, pch = 16, cex = 2, main="Growing Season RWCIP", ylab=expression(delta^2~H[f]), xlab = expression (delta^2~H[p]),cex.lab=2)
abline(RTHU_GS$rma_coefs)
legend("topleft", legend=c("2004", "2007"), fill = c("blue", "yellow"), title = "RTHU dataset")


plot(zenal2018data$Precip, zenal2018data$Tissue, pch = 16, cex = 2, main="Zenzel, 2018", ylab=expression(delta^2~H[f]), xlab = expression (delta^2~H[p]),cex.lab=2)
abline(ZenzelRMA$rma_coefs)



#print out regression stats ##
cat("Model \t \t intercept \t slope \t r2\n")
cat("Growing Season", RTHU_GS$rma_coefs, RTHU_GS$correlation,"\n")
cat("RWCIP MA\t\t", RTHU_MA_RCWIP$rma_coefs, RTHU_MA_RCWIP$correlation,"\n")
cat("RCWIP GS\t\t", RTHU_GS_RCWIP $rma_coefs, RTHU_GS_RCWIP $correlation,"\n")
cat("RCWIP2 MA\t\t", RTHU_MA_RCWIP2$rma_coefs, RTHU_MA_RCWIP2$correlation,"\n")
cat("Zenzel, 2018", ZenzelRMA$rma_coefs, ZenzelRMA$correlation, "\n")

