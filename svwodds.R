#################################################################################
#
#
# svpodds function - function to calculate odds ratio from probability surface.  
#  Original Code from Steven Van Wilgenburg (2011)
# Modified G. Koehler 2025
#################################################################################

svwodds <- function(probsurface, odds){
  
#  prob_stack <- stack(probsurface) #stack probsurface to work with raster library
#  prob_stack <- probsurface
  
  ##################################################################
  #### EXTRACT the estimated probabilities
  
  ### FIRST, extract values to a data frame, excluding na (not applicable) values
  d = na.omit(data.frame(values(probsurface)))# ; summary(d)
  
  ###reshape the data frame, creating a variable called "layer" to act as a factor
  #### "varying" refers to the columns
  e<-reshape(d, direction="long", v.names="density",timevar = "layer", varying=c(1:ncol(d)),sep="")
  
  ################################################### INPUT REQUIRED HERE ################################################################
  ####                                    ENTER the desired cumulative probability bound (e.g. 0.33 to select upper 67%)                                                                  #####
  ########################################################################################################################################
  
 # odds<-0.33333333  ### MODIFY ONLY THE RIGHT HAND SIDE of arrow (i.e. change 0.33 to desired value)
  
  #####    Function to fit splines to estimate the probability densities at which the THRESHOLD CUMULATIVE probability is reached     ####
  fun <- function(x) predict(smooth.spline(cumsum(sort(x)),sort(x), spar=0.1),odds)  #### function for fitting a spline curve to cumulative probabilities
  
  ### below applies the spline function to the data, within each layer
  cutoffs<- by(e$density,e$layer,fun)
  cutoffs<-as.matrix(unlist(cutoffs))
  
  #### output from the spline function is an odd format, including both x (the cumulative probability of interest (e.g. 0.33333) and y, the following function removes, the unnecessary data
  cutoffs<-as.matrix(cutoffs[-which( cutoffs[,1] == odds)])
  colnames(cutoffs) <- "cutoff"
  ###view the cutoff for each layer
  # cutoffs
  
  
  ################################## RECLASSIFY each map to binary surface #################################################
  s<-list() ### storage object
  for (i in 1:nrow(cutoffs)){
    s[i] = reclassify(probsurface[[i]],c(-Inf,cutoffs[[i]],0,cutoffs[[i]],Inf,1))
  }
  
  ############ CALCULATE THE ASSIGNED ORIGINS FOR WHOLE SAMPLE
  ### change the list of individual assignments into a "raster stack"
  z<-stack(s)
  i<-maxValue
  origins <- overlay(z, fun=function(i){return(sum(i))} )
  
  return(origins)
  
  ###################################################################
  ### END CODE
  
  
  
  
}
