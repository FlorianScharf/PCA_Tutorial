############### STEP 5b: Analysis of Jackknife Samples ###############
# This script analyzes the results of the jackknife samples and provides inferential statistics
# for latency differences.
# Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
# Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig

# clear workspace
rm(list = ls())

############### Settings ###############

# Relative amplitude criterion (we used 0.8, i.e., 80% of peak amplitude)
crit <- 0.8
# Provide sampling rate
srate <- 500
# time of epoch beginning (in seconds)
xmin <- -0.2
# number of sampling points per epoch
pnts <- 500

# Please list jackknife fit files to be analyzed 
jackknifeFits <- c("results/05a_jackknifing_step_1/rotfit_ad23_jkpca.Rdata", 
                   "results/05a_jackknifing_step_1/rotfit_ch21_jkpca.Rdata")

# Please list the respective total sample efa fits 
# attention: you need to list as many files here 
# as you list jackknifeFits
totalFits <- c("results/02bc_rotation_score/rotfit_ad23_geomin0.01.Rdata", 
               "results/02bc_rotation_score/rotfit_ch21_geomin0.01.Rdata")

# For which components should the latency be extracted?
# (the index refers to loadings in totalFits)
facIdx <- c(3, 5)

#### Adults 
# iFacIdx <- 2 # Adults P2
# iFacIdx <- 5 # Adults early P3a
# iFacIdx <- 3 # Adults late P3a
# iFacIdx <- 1 # Adults LDN

#### Children 
# iFacIdx <- 2 # Children P2 
# iFacIdx <- 6 # Children early P3a
# iFacIdx <- 5 # Children late P3a
# iFacIdx <- 1 # Children LDN P3a
###############

############### Preparations ###############

# load packages
library(BayesFactor)
library(ggplot2)
library(reshape)

# a convenience function to extract the relevant information 
# from the jackknife fits
evalIndFits <- function(rotFit, rotFitAll, iFacIdx, crit){
  
  sapply(rotFitAll, function(iFit){

    # find factor in subsample by similarity with factor of interest
    factorIdx <- which.max(cor(rotFit$loadings[,iFacIdx], iFit$loadings))
    # keep the loadings of this factor and unstandardize
    tmpLoad <- iFit$loadings[,factorIdx] * iFit$varSD
    
    # determine peak amplitude
    peakAmp <- max(tmpLoad)
    # find the first sampling point above the crit threshold
    latIdx <- min(which(tmpLoad >= peakAmp * crit))

    # append peak latency and loadings in a vector
    c(latIdx = latIdx,tmpLoad)
    }, simplify = "array")
  
}
###############

############### Extract the individual peak latencie estimates ###############

# The following code loads all jackknifeFits and extracts individual 
# latency estimates for the requested components
# in addition, a diagnostic plot of the subsample loadings is provided
# so that problems with factor alingment or stability of the factor solution can be detected
indLatEst <- lapply(jackknifeFits, FUN = function(iFit){
  # get index of current fit 
  iFitIdx <- which(iFit == jackknifeFits)
  # get index of factor to be analyzed
  iFacIdx <- facIdx[iFitIdx] 
  # load current fit 
  load(iFit)
  # load respective total fit object
  load(totalFits[iFitIdx])
  
  # extract sample size
  N <- length(rotFitAll) # there are as many fits as participants
  
  # apply convenience function to extract loadings and peak latencies
  out <- evalIndFits(rotFit, rotFitAll, iFacIdx, crit)
  # convert indices to seconds
  lat <- (out["latIdx",] - 1) / srate + xmin
  # Estimate individual latencies; Smulders, 2010, eq. (1)
  indLat <- N * mean(lat) - (N - 1) * lat
  
  # for a more generea explanation of the underlying mathematics, see, e.g.,
  # Chapter 10 in:
  # Efron, B., & Hastie, T. (2016). 
  # Computer age statistical inference: Algorithms, evidence, and data science. 
  # https://doi.org/10.1017/CBO9781316576533
  
  
  # restructure the loadings to make plotting easier
  loadings = data.frame(timeAxis = (0:(pnts - 1)) / srate + xmin, out[-1,])
  loadings = melt(loadings, id.vars = "timeAxis", variable_name = "participant")
  diagPlot <- ggplot(data = loadings, aes(x = timeAxis, y = value, color = participant)) +
    geom_line()
  
  # name the latencies after the current fit
  # to avoid confusion
  indLat <- as.matrix(indLat)
  colnames(indLat) <- iFit
  
  # return the individual latencies, the plot, and the subsample loadings
  list(indLat = indLat, lat = lat, diagPlot = diagPlot, loadings = loadings)
})

# name list elements by the files inserted to exclude confusion from
# indexing
names(indLatEst) <- paste0(jackknifeFits, "_facIdx_", facIdx)
###############

############### Show diagnostic plots ###############
# ideally, the individual loadings should be highly similar
# any indication that the factor structure varies drastically between
# subsamples should be inspected closely!
indLatEst[[1]]$diagPlot
indLatEst[[2]]$diagPlot
###############

############### Subject latencies to a test ###############
# frequentist and bayesian t tests
t.test(indLatEst[[1]]$indLat, indLatEst[[2]]$indLat)
ttestBF(indLatEst[[1]]$indLat, indLatEst[[2]]$indLat)
###############

############### Notes on further use of indices ###############

# indLatEst is a so-called list object 
# it contains as many elements as fit objects were provided 
# each element can be accessed by simple indexing:
indLatEst[[1]]
# or altenatively by the names of the jackknife fit files:
indLatEst$rotfit_ad25_jkpca.Rdata_facIdx_5

## Each element contains three elements:
# individual latency estimates
indLatEst[[1]]$indLat
# the diagnostic plot of the subsample loadings
indLatEst[[1]]$diagPlot
# a data frame with all subsample loadings 
indLatEst[[1]]$loadings
# note that participant note the participant who was EXCLUDED
# in that specific jackknife subsample

###############
