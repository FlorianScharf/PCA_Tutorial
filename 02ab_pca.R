############### STEP 2: ESTIMATION OF THE FACTOR MODEL ###############
# This script conducts the initial estimation of the unrotated temporal pca model.
# Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
# Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig

## Check if necessary packages are installed and if not
# install them
if(!require(psych)) install.packages("psych")
if(!require(MASS)) install.packages("MASS")

## Load packages
library(psych)
library(MASS)

## empty work space
rm(list=ls())

## Load pre-processed data
load("results/01_data_import/erpdata.Rdata")

## Select group to be analyzed
# ad = adults; ch = children
group = "ch"
# Please remember that this script needs to be run once for each group.


## Prepare data set 
# Select a subset of the data based on the chosen group
# delete all variable except for the sampling points (necessary for the pca functions)
age_group_data = as.matrix(erpdata[erpdata$group == group, -c(1:4)])

############### 2a:Determination of the Number of Factors ############### 

## Compute correlation matrix of the sampling points 
# to be used for the Empirical Kaiser Criterion (EKC)
# Attention: R should ONLY contain correlations between the variables
# which shall enter the pca. Please make sure to remove any indicator variables
# such as participant or condition ID before computing R!
R = cor(age_group_data)

## Determine number of factors using 
# This function is taken from the package EFAtools. 
source('tools/EKC.R')
res_ekc = EKC(R, N = nrow(age_group_data))

## Save the resulting number of factors
nFac = res_ekc$n_factors

## Print the number of factors
nFac

## Plot the variance explained by each factor in the initial
# solution, we cut the solution at 40 factors just to make the
# cut point more visible.
plot(1:ncol(age_group_data), res_ekc$eigenvalues,
     xlab = "Factor", ylab = "Variance Explained",
     main = ifelse(group == "ch", "Children", "Adults"),
     xlim = c(0,40), pch = 16,
     col = (res_ekc$references <= res_ekc$eigenvalues) + 1)
lines(1:ncol(age_group_data), res_ekc$references, lty = 2, lwd = 3,
      col = "blue")

abline(v = nFac)
text(x = nFac, y = 100, pos = 2, cex = 0.8,  
     labels = paste0("Number of Factors\nto be extracted: ", nFac))

### Optional code for 
## Horn's parallel test as implemented by Dien's ERP PCA toolkit
# Eigen values of the data
S = cov(age_group_data)
eigenData = eigen(S)$values

# Generate uncorrelated random data with the same
# standard deviation as the original data
randData = matrix(rnorm(n = dim(age_group_data)[1] * dim(age_group_data)[2], 
                        sd = sqrt(diag(S))), 
                  nrow = dim(age_group_data)[1])

# Covariance matrix and eigen values of the random data
Srand = cov(randData)
eigenRandData = eigen(Srand)$values

# How many eigen values of the real data are
# above their random counterparts?
# sum(eigenRandData < eigenData) 
min(which(eigenRandData >= eigenData)) - 1
# Uncomment and replace nFac if you want to overwrite
# the EKC result and use a different number of factors
# nFac <- 9999

## A basic plot of parallel analysis
# we only plot the first 50 eigen values
# to keep the figure comprehensible
plot(eigenData[1:50],
     xlab = "Factors",
     ylab = "Eigen Values")
lines(eigenRandData, col = "red")

############### 

############### 2b: Estimation of unrotated factor loadings ############### 

## Estimate unrotated factor loadings
# This function is adapted from the package psych. We slightly modified it to 
# improve estimation speed. The original function tries to compute a number of
# fit indices which are not relevant for our purposes. However, this computation
# is very slow for data sets as large as ours and it often results in irrelevant 
# warnings about non-convergence of the fit index estimation. Otherwise, the
# function is not changed.
source("tools/fa_simplified.R")

# The argument covar = TRUE is set in order to produce an unstandardized 
# solution as intended.
pcaFit = fa_simplified(r = age_group_data, nfactors = nFac, rotate = "none", covar = TRUE)

# We add a number of relevant descriptive statistics to the fit object
# which we use again in later steps.

S = cov(age_group_data)
Sinv = ginv(S)
Rinv = ginv(R)
Var = diag(S)
varSD = sqrt(Var)

pcaFit$S = S # Sample Covariance Matrix
pcaFit$Sinv = Sinv # (Moore-Penrose) Inverse of the Sample Covariance Matrix
pcaFit$R = R # Sample Correlation Matrox
pcaFit$Rinv = Rinv # (Moore-Penrose) Inverse of the Sample Covariance Matrix
pcaFit$varSD = varSD # SDs of the Sampling Points
pcaFit$Var = Var # Variances of the Sampling Points
pcaFit$group = group

# Finally, we save the unrotated solutions in a separate file for the
# next analysis steps. 
# The file name is automatically adjusted to include group name
# and number of factors.
save(pcaFit, file = paste0("results/02ab_pca/pcafit_", group, nFac, ".Rdata"))
###############