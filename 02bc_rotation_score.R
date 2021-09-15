############### STEP 2: ESTIMATION OF THE FACTOR MODEL ###############
# This script conducts the factor rotation and factor scoring steps of the analysis.
# Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
# Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig

# empty workspace
rm(list=ls())

## Check if necessary packages are installed and if not
# install them
if(!require(GPArotation)) install.packages("GPArotation")

## Load necessary packages
library(GPArotation)
library(psych)
source("tools/geominQ_multstart.R") # load custom rotation function

## Load data
# Load the results file from the previous script.
# Please remember that this script needs to be run once for each group.
# load("results/02ab_pca/pcafit_ad23.Rdata")
load("results/02ab_pca/pcafit_ch21.Rdata")

############### 2b: Estimation of rotated factor loadings ###############

## Standardize factor loadings before rotation
pcaFit$loadings = pcaFit$loadings / pcaFit$varSD

## Use Geomin rotation with 30 random start values.
# Please be aware that this function can take a while. 
# This is a custom function which implements a procedure 
# available from the package GPArotation. 
# Please see tools/geominQ_multstart.R for details
rotFit <- geominQ_multstart(A = pcaFit$loadings,  # unrotated loadings
                            delta = 0.01,     # rotation parameter (geomin epsilon)
                            # Note: We decided to name all parameters consistently with
                            # the GPArotation package despite its deviation from the
                            # conventional naming epsilon for this parameter.
                            normalize = F,     # No additional standardization
                            rand.start = T,    # Use multiple random starts
                            start.values = 30, # Number of random starts
                            maxit = 50000,     # Number of iterations 
                            # Note: After this number of iterations, the
                            # function stops trying to estimate the parameters
                            # from this random starting values.
                            eps = 1e-5)        # Level of accuracy to determine convergence
                            # Note: This means that the rotation is declared successfully 
                            # converged when the criterion changes less than eps between
                            # two iterations. We recommend against using lower values 
                            # (this is the GPArotation default) since it can prevent
                            # the algorithm from converging.
                            # You can choose a higher values (e.g., 1e-3) if the
                            # rotation takes too long but we recommend using the default
                            # value if you are seriously interested in the results.

### IMPORTANT: Please increase maxit in case you encounter many warnings
### regarding non-converged solutions.

#### OPTIONAL: If you want to apply a Promax rotation.
# Promax does not rely on a general gradient projection algorithm
# and does not require random starts (but it tends to conflate factors
# earlier than Geomin):

# We set the rotation parameter to 3, mimicking the behavior of Dien's ERP
# PCA Toolkit and the commercial software SAS, 
# because this value is well-justified by simulation studies.
# Choosing a value of 2 will allow for solutions with higher temporal overlap,
# choosing higher values will increase the preference for too simple structures.

# rotFit <- Promax(pcaFit$loadings, m = 3, normalize = TRUE)

## Transfer variances and standard deviations into the new fit object.
rotFit$varSD = pcaFit$varSD
rotFit$Var = pcaFit$Var
rotFit$group = pcaFit$group

##  Sort factors by variance explained

# Compute unstandardized loadings
L <- rotFit$loadings * rotFit$varSD
# Compute proportion of variance explained by each factor
facVar <- diag(rotFit$Phi %*% t(L) %*% (L)) / sum(diag(pcaFit$S)) 


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

###############

############### 2c: Estimation of factor scores ###############

# Load raw data again (because they are necessary for factor scoring)
load("results/01_data_import/erpdata.Rdata")

# Again, need to delete columns not containing data from sampling points
age_group_data = as.matrix(erpdata[erpdata$group == pcaFit$group, -c(1:4)])


## Estimate factor scores
# The formula is taken from DiStefano et al. (2009; Appendix 2).
# Note that "age_group_data" contains the unstandardized ERP individual averages 
# of the respective group.
# Therefore, the factor scores are not centered.
FacScr = age_group_data %*% pcaFit$Sinv %*% (rotFit$loadings * rotFit$varSD) %*% rotFit$Phi

# Rinv ... generalized inverse of the correlation matrix of the sampling points
# rotFit$loadings ... standardized factor loadings
# rotFit$Phi ... factor correlation matrix

# rename columns because they are factors now
colnames(FacScr) = colnames(paste0("MR", 1:pcaFit$factors))

# Show descriptive statistics for the factor scores
psych::describe(FacScr)

## We save the factor scores in a separate object and "reunite" them
# with the indicator variables describing participant, condition and 
# electrode site
scores = data.frame(erpdata[erpdata$group == pcaFit$group, 1:4], FacScr)

## We save all estimated objects for later steps.
save(pcaFit, rotFit, scores, file = paste0("results/02bc_rotation_score/rotfit_", pcaFit$group, pcaFit$factors, "_geomin0.01.Rdata"))


## Export to MATLAB
# This is optional. If you prefer to work in MATLAB for visualization purposes,
# this is how you export the results to a mat-File.
library(R.matlab)
writeMat(paste0("results/02bc_rotation_score/rotfit_", pcaFit$group, pcaFit$factors, "_geomin0.01.mat"), loadings = unclass(rotFit$loadings), phi = unclass(rotFit$Phi), scores = FacScr, varSD = unclass(pcaFit$varSD))

###############


