############### STEP 1: Data import ###############
# Import data from MATLAB and save as Rdata
# Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
# Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig


## Empty workspace to start with a clean plate.
rm(list=ls())

## Check if necessary packages are installed and if not
# install them
if(!require(R.matlab)) install.packages("R.matlab")
if(!require(ggplot2)) install.packages("ggplot2")
library(R.matlab)
library(ggplot2)


### Read participant averages 
avrdata = readMat("data/avrdata.mat")

### Combine the data set into a "nice" and labeled R dataframe
colnames(avrdata$data) = paste0("erp_", 1:500)
colnames(avrdata$dataIdx) = c("group", "cond", "subj", "chan")

erpdata = data.frame(cbind(avrdata$dataIdx, avrdata$data))
erpdata$group = factor(erpdata$group, labels = c("ad", "ch"))
erpdata$cond = factor(erpdata$cond, labels = c("sta", "nov"))
erpdata$subj = factor(erpdata$subj)
erpdata$chan = factor(erpdata$chan, labels = c("Fp1", "Fz", "F3", "F7", "IO1", "FC5", "FC1", "C3", "T7", "M1", "CP5", "CP1", "Pz", "P3", "P7", "LO1", "Oz", "LO2", "P4", "P8", "M2", "CP6", "CP2", "Cz", "C4", "T8", "FC6", "FC2", "F4", "F8", "FP2"))
###

### We also save some data characteristics 
fs = 500 # sampling rate
xmin = -0.2 # baseline
pnts = 500 # epoch duration
lat = (1:pnts - 1) / fs + xmin
###

### Save as csv and as Rdata-File
write.csv(erpdata, file = "results/01_data_import/erpdata.csv")
save(erpdata, fs, xmin, pnts, lat, file = "results/01_data_import/erpdata.Rdata")
