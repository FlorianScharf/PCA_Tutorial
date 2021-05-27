############### STEP 5a: Jackknife resampling ###############
# This script conducts the jackknife resampling and safes the rotated EFA results for all samples.
# Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
# Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig

#### NOTE:
#### This script runs for several hours on typical hardware. Do only run when 
#### you are prepared to wait for the results this long.

## Load necessary packages
if(!require(foreach)) install.packages("foreach")
if(!require(doParallel)) install.packages("doParallel")
library(foreach)
library(doParallel)
library(psych)
library(GPArotation)

## Setup parallel environment
# The resampling is conducted in a parallel environment
# each single core estimates the EFA model for a different jackknife sample.
# The higher the number of cores the faster this will be.
# We recommend to use as many cores as your computer has available minus one 
# to prevent overload. If you use a computing cluster with many cores:
# The script will not benefit from more than 32 cores because this is the number
# of jackknife samples.

# Return the number of cores available
detectCores(logical = TRUE)

# Set number of cores
cores <-  32 # change manually

# Prepare parallel environment
cl <- makeCluster(cores[1], outfile = "")
registerDoParallel(cl)

# Load raw data
load("results/01_data_import/erpdata.Rdata")

# Group labels to iterate through
groups <- c("ad", "ch") 

# How many factors should be extracted per group?
# This should match the number of factors from the EFA over all participants.
nFac <- c(23,21)

# go through both groups
for (iGroup in groups){ 
  
  # set number of factors to be extracted accordingly
  iNfac <- nFac[which(groups == iGroup)]
  
  rotFitAll = foreach(subj = levels(droplevels(erpdata[erpdata$group == iGroup,]$subj)), .combine=cbind, .packages = c("psych", "GPArotation")) %dopar% {
    
    # The following code essentially replicates the estimation of the EFA model 
    # in step 2 for every jackknife subsample. 
    # We only increased the number of maximum iterations for the rotation as well
    # as the number of random starts to prevent suboptimal results due to 
    # local optima or non-convergence in single samples.
    
    source("tools/myFA.R")
    source("tools/geominQ_multstart.R")
    
    data = as.matrix(erpdata[erpdata$group == iGroup & erpdata$subj != subj, -c(1:4)])
    
    S = cov(data)
    Var = diag(S)
    varSD = sqrt(Var)
    
    efaFit = fa_simplified(data, nfactors = iNfac, rotate = "none", covar = TRUE)
    efaFit$loadings = efaFit$loadings / varSD # as in ERP PCA Toolkit
    
    rotFit <- geominQ_multstart(A = efaFit$loadings,  # unrotated loadings
                                delta = 0.01,     # rotation parameter (geomin epsilon)
                                # Note: We decided to name all parameters consistently with
                                # the GPArotation package despite its deviation from the
                                # conventional naming epsilon for this parameter.
                                normalize = F,     # No additional standardization
                                rand.start = T,    # Use multiple random starts
                                start.values = 100, # Number of random starts
                                maxit = 500000,     # Number of iterations 
                                # Note: After this number of iterations, the
                                # function stops trying to estimate the parameters
                                # from this random starting values.
                                eps = 1e-5)        # Level of accuracy to determine convergence
    # Note: This means that the rotation is declared sucessfully 
    # converged when the criterion changes less than eps between
    # two iterations. We recommend against using lower values 
    # (this is the GPArotation default) since it can prevent
    # the algorithm from converging.
    # You can choose a higher values (e.g., 1e-3) if the
    # rotation takes too long but we recommend using the default
    # value if you are seriously interested in the results.
    
    ## Transfer variances and standard deviations into the new fit object.
    rotFit$varSD = varSD
    rotFit$Var = Var
    rotFit$group = iGroup
    
    ##  Sort factors by variance explained
    # Compute factor variances
    facVar = apply(diag(rotFit$Var) %*% rotFit$loadings * rotFit$loadings %*% rotFit$Phi, MARGIN = 2, FUN = sum) / sum(rotFit$Var) # From Dien ep_doPCA.m
    # Return indices of the factors ordered by the variance explained
    alignment = order(facVar, decreasing = TRUE)
    # reorder columns of factor loadings matrix in descending order of variance
    # explained
    rotFit$loadings = rotFit$loadings[, alignment]
    # reorder factor correlation matrix as well
    rotFit$Phi = rotFit$Phi[alignment, alignment]
    
    ## Flip mostly negative factor loadings
    flip = sign(colSums(rotFit$loadings)) # flip now contains -1 and 1, we want to turn all the -1s
    rotFit$loadings = rotFit$loadings %*% diag(flip) # post multiplying with this turns the factor loadings
    rotFit$Phi = diag(flip) %*%  rotFit$Phi %*% diag(flip) # turning factor correlations
    
    list(rotFit)
    
  }
  
  save(rotFitAll, file = paste0("results/05a_jackknifing_step_1/rotfit_", iGroup, iNfac, "_jkpca.Rdata"))
  
}


stopCluster(cl)