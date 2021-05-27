############### STEP 3b: Visual Inspection of Topography ###############
# This script enables visual inspection of the factors.
# Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
# Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig

## Check if necessary packages are installed and if not
# install them
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(RColorBrewer)) install.packages("RColorBrewer")
if(!require(reshape)) install.packages("reshape")
if(!require(xlsx)) install.packages("xlsx")
if(!require(remotes)) install.packages("remotes")
if(!require(R.matlab)) install.packages("R.matlab")
if(!require(eegUtils)) remotes::install_github("craddm/eegUtils")
if(!require(dplyr)) install.packages("dplyr")
if(!require(grid)) install.packages("grid")
if(!require(gridExtra)) install.packages("gridExtra")

library(grid)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(xlsx)
library(eegUtils)
library(R.matlab)

## Load some convinience functions for data handling
source("tools/topoplot_functions.R")

###### Data preparations ########

## Load data
load("results/01_data_import/erpdata.Rdata")

## Compute grand average across participants per group, channel and condition
allGAVRS <- aggregate(. ~ cond + chan + group, data = erpdata, FUN = mean)

## Load the data from both EFAs and put the results into separate objects
load("results/02bc_rotation_score/rotfit_ad23_geomin0.01.Rdata",  temp_env <- new.env())
adEFA <- as.list.environment(temp_env)
load("results/02bc_rotation_score/rotfit_ch21_geomin0.01.Rdata",  temp_env <- new.env())
chEFA <- as.list.environment(temp_env)

## Rename the factor score matrix columns 
colnames(adEFA$rotFit$loadings) <- paste0("F", 1:adEFA$efaFit$factors)
colnames(chEFA$rotFit$loadings) <- paste0("F", 1:chEFA$efaFit$factors)
colnames(adEFA$scores)[-c(1:4)] <- paste0("F", 1:adEFA$efaFit$factors)
colnames(chEFA$scores)[-c(1:4)] <- paste0("F", 1:chEFA$efaFit$factors)


## average factor score per condition and channel
adEFA$average_scores <- aggregate(. ~ cond + chan, data = adEFA$scores, FUN = mean)
chEFA$average_scores <- aggregate(. ~ cond + chan, data = chEFA$scores, FUN = mean)

allEFAs <- list(ad = adEFA, ch = chEFA)

## load subject averages 
# we will use these to have the channel locations in a format
# easily understanbale for eegUtils
allAVRs <- list(
  ad = list(sta = eegUtils::import_set("data/ad-sta.set"), nov = eegUtils::import_set("data/ad-nov.set")),
  ch = list(sta = eegUtils::import_set("data/ch-sta.set"), nov = eegUtils::import_set("data/ch-nov.set"))
)


########


###### Make combined time course and topography figures ########
# Choose a group to be plotted 
iGroup <- "ch"

# Extract EFA results for the current group
iEFA <- allEFAs[[iGroup]]

# Choose a factor to be plotted 
iFactor <- "F2"

# Choose an electrode to be plotted
# here we automatically find maximum electrode
# but it would be better to use a pre-selected electrode
# based on the expected factor structure
iEl <- iEFA$average_scores$chan[which.max(abs(iEFA$average_scores[, iFactor]))]

# unstandardize rotated loadings
iEFA$rotFit$loadings <- iEFA$rotFit$loadings * iEFA$efaFit$varSD

# find peak loading
iPeakLoad <- max(iEFA$rotFit$loadings[, iFactor]) 

# Extract scores for current factor
iScores <- iEFA$average_scores[, c("cond", "chan", iFactor)]

# Extract the time codes of the epoch
times <- unique(allAVRs[[iGroup]][["sta"]]$timings$time)


### The following code restructures the extracted information 
# so that it fits the requirements of ggplot, most importantly,
# the data are reshaped into long format.
# The details can be found in tools/topoplot_functions.

# Get grand average to long format for plotting
iGAVR <- gavr2long(allGAVRS, iEl = iEl, iGroup = iGroup, times = times)


# Compute reconstructed ERP for each observation for this factor
# and reshape it to long-format.
iERPest <- makeERPest(iScores = iScores,
                      iEFA = iEFA, 
                      iFactor = iFactor, 
                      iEl = iEl,
                      times = times)


## 

## Plot grand averages separately for each condition
time_course <- ggplot(iGAVR,
       aes(x = times,
           y = value, 
           color = cond,
           group = cond)) +
  geom_line(size = 2, alpha = 0.3) + # plot grand average
  geom_line(data = iERPest, size = 2) + # plot reconstructed ERP
##### This part of the code only serves the purpose of making the plot look nicer. ####
  scale_color_manual(values =  c(brewer.pal(n = 11, "RdYlBu")[c(10,2)], "grey49")  ) +
  #theme_linedraw()  + # if you prefer a white background uncomment this
  scale_y_reverse() +  
  labs(x = "Time [s]" , y = "Voltage [µV]", 
       title = paste0("Group: ", iGroup,", Factor: ", iFactor,", ", iEl, ", Peak Latency: ", times[which.max(iEFA$rotFit$loadings[, iFactor])], " s"),
       color = "Condition") +
  theme(plot.title = element_text(color="black", size=12, hjust = 0.5, 
                                  face = "bold"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13, face = "bold"),
        legend.position = "bottom")
#####

time_course

## We use a trick to provide the reconstructed ERPs in a way compatible with 
## the eegUtils-package:
# We replace the original data with the reconstructed data.

signals <- efa2eeg(iEFA = iEFA, 
                   allAVRs = allAVRs,
                   iGroup = iGroup,
                   iFactor = iFactor)

# Calculate the limits of the topography color scale
topoLimit <- max(ceiling(sapply(signals, function(x) max(abs(x)))))

# Plot topographies for all conditions for the selected factor
allTopos <- lapply(names(signals), function(i){
  iSignals <- signals[[i]]
  tmp <- allAVRs[[iGroup]]$nov
  tmp$signals <- as_tibble(iSignals)
  
  topoplot(tmp, 
           interp_limit = "skirt",
           grid_res = 300,
           limits = c(-1, 1) * topoLimit, 
           scaling = 0.9,
           highlights = iEl,
           contour = TRUE
  ) +
    ggtitle(i) +
    theme(plot.title = element_text(color="black", size=18, hjust = 0.5,
                                    face = "bold"),
          legend.text = element_text(size = 12)) 



  
})

names(allTopos) <- names(signals)

# Look at the topographies within R
allTopos$sta 
allTopos$nov
allTopos$`nov - sta`

# Combine the time course plot and the topography plots into a single eps-file
cairo_ps(filename = paste0("results/03b_topoplot_selectFactor/topoplot_", iFactor,".eps"),
           width = 14, height = 4)

    grid.arrange(time_course,
                 allTopos$sta + guides(fill = F),
                 allTopos$nov + guides(fill = F),
                 allTopos$`nov - sta`,
                 widths = c(4, 2, 2, 2.85),
                 nrow = 1, ncol = 4
    )

dev.off()

############ Plot Grand Average ###########
# Choose electrode site:
iEl <- "Cz"

## Plot grand averages separately for each condition
gavr_plot_ad <- ggplot(gavr2long(allGAVRS, iEl = iEl, iGroup = "ad", times = times),
                       aes(x = times,
                           y = value, 
                           color = cond,
                           group = cond)) +
  geom_line(size = 1.5) + # plot grand average
  ##### This part of the code only serves the purpose of making the plot look nicer. ####
scale_color_manual(values =  c(brewer.pal(n = 11, "RdYlBu")[c(10,2)], "grey49") ) +
  #theme_linedraw()  + # if you prefer a white background uncomment this
  ylim(12,-8) +   
  labs(x = "Time [s]" , y = "Voltage [µV]", 
       title = paste0("Adults", ", Electrode: ", iEl),
       color = "Condition") +
  theme(plot.title = element_text(color="black", size=14, hjust = 0.5, 
                                  face = "bold"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13, face = "bold"),
        legend.position = "bottom")

gavr_plot_ch <- ggplot(gavr2long(allGAVRS, iEl = iEl, iGroup = "ch", times = times),
                      aes(x = times,
                          y = value, 
                          color = cond,
                          group = cond)) +
  geom_line(size = 1.5) + # plot grand average
  ##### This part of the code only serves the purpose of making the plot look nicer. ####
  scale_color_manual(values =  c(brewer.pal(n = 11, "RdYlBu")[c(10,2)], "grey49") ) +
  #theme_linedraw()  + # if you prefer a white background uncomment this
  ylim(12,-8) +  
  labs(x = "Time [s]" , y = "Voltage [µV]", 
       title = paste0("Children, Electrode: ", iEl),
       color = "Condition") +
  theme(plot.title = element_text(color="black", size=14, hjust = 0.5, 
                                  face = "bold"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 12),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13, face = "bold"),
        legend.position = "bottom")

gavr_plot_ad
gavr_plot_ch


cairo_ps(filename = paste0("results/03b_topoplot_selectFactor/Grand_averages.eps"),
         width = 10, height = 4)

grid.arrange(gavr_plot_ch,
            gavr_plot_ad,
             widths = c(5, 5),
             nrow = 1, ncol = 2
)

dev.off()
#####


