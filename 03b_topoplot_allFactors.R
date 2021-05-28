############### STEP 3b: Visual Inspection of Topography ###############
# This script enables visual inspection of the factors.
# Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
# Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig

# empty workspace
rm(list=ls())


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
# easily understandable for eegUtils
# Please note that these set-files do NOT contain epoched raw data.
# Instead, each "epoch" in these datasets is a participant average.
# see also: https://github.com/widmann/grandaverage/blob/master/pop_grandaverage.m
allAVRs <- list(
  ad = list(sta = eegUtils::import_set("data/ad-sta.set"), nov = eegUtils::import_set("data/ad-nov.set")),
  ch = list(sta = eegUtils::import_set("data/ch-sta.set"), nov = eegUtils::import_set("data/ch-nov.set"))
  )


########


###### Make combined time course and topography figures ########

### This generates the plots for all actors automatically and generates a PDF file.
allPlots <- lapply(c("ad", "ch"), function(iGroup){
  
  # Extract EFA results for the current group
  iEFA <- allEFAs[[iGroup]]
  
  # unstandardize rotated loadings
  iEFA$rotFit$loadings <- iEFA$rotFit$loadings * iEFA$efaFit$varSD
  
  # compute maximal factor score to scale all factor topographies equally
  allPeaks <- as.matrix(iEFA$average_scores[,-c(1:4)]) %*% diag(apply(iEFA$rotFit$loadings, 2, max))
  topoLimit <- ceiling(max(abs(allPeaks)))

  # compute ERP min and max for equal scaling across factors
  voltageLimits <- c(ceiling(max(allGAVRS[allGAVRS$group == iGroup, -c(1:4)])),
                     floor(min(allGAVRS[allGAVRS$group == iGroup, -c(1:4)])))
                     
  
  lapply(colnames(allEFAs[[iGroup]]$rotFit$loadings), function(iFactor){
  
    # Choose an electrode to be plotted
    # here we automatically find maximum electrode (as is done by the ERP PCA Toolkit)
    # but it would be better to use a pre-selected electrode
    # based on the expected factor structure
    iEl <- iEFA$average_scores$chan[which.max(abs(iEFA$average_scores[, iFactor]))]

  

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
      geom_line(size = 1.25, alpha = 0.3) + # plot grand average
      geom_line(data = iERPest, size = 1.25) + # plot reconstructed ERP
      ##### This part of the code only serves the purpose of making the plot look nicer. ####
     scale_color_manual(values =  c(brewer.pal(n = 11, "RdYlBu")[c(10,2)], "grey49")  ) +
      #theme_linedraw()  + # if you prefer a white background uncomment this
      ylim(voltageLimits) +
      #scale_y_reverse() +  
      labs(x = "Time [s]" , y = "Voltage [µV]", 
           title = paste0(ifelse(iGroup == "ch", "Children", "Adults"), 
                          ", Electrode: ", iEl,
                          ", Factor: ", iFactor,
                           ",\nPeak Latency: ", sprintf("%.3f",times[which.max(iEFA$rotFit$loadings[, iFactor])]), " s"),
           color = "Condition") +
      theme(plot.title = element_text(color="black", size=13, #hjust = 0.5, 
                                      face = "bold"),
            axis.text = element_text(color = "black", size = 12),
            axis.title = element_text(color = "black", size = 12),
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 12, face = "bold"),
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
    # topoLimit <- max(ceiling(sapply(signals, function(x) max(abs(x)))))
    
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
    
    list(
      time_course =  time_course,
      topo_sta = allTopos$sta + guides(fill = F),
      topo_nov = allTopos$nov + guides(fill = F),
      topo_diff = allTopos$`nov - sta`)
    
  })
})

## allPlots now contains the topography and time course plots
## for both groups, all factors and all conditions

##### We print these plots into convenient PDF files:
names(allPlots) <- c("ad", "ch")
# Determine order in which the plots are returned
layout_matrix = matrix(1:24, nrow = 6, ncol = 4, TRUE)

# adult EFA
pages <- marrangeGrob(unlist(allPlots$ad, recursive = F), # put all adult plots in one long vector
             ncol = 4, # 4 columns (do not change!)
             nrow = 6, # 6 rows (change if you want fewer factors per page)
             heights = rep(3,6), # heights of the rows
             widths = c(3.2, 2, 2, 3), # widths of the columns
             layout_matrix = layout_matrix)
# print on A3 (font size etc. fit better at this ratio)
ggsave("results/03b_topoplot_allFactors/adEFA.pdf", pages, width = 1.4*210, height = 1.4*297, units = "mm")

# child EFA
pages2 <- marrangeGrob(unlist(allPlots$ch, recursive = F), 
                      ncol = 4, 
                      nrow = 6,
                      heights = rep(3,6), # heights of the rows
                      widths = c(3.2, 2, 2, 3), # widths of the columns
                      layout_matrix = layout_matrix)

ggsave("results/03b_topoplot_allFactors/chEFA.pdf", pages2, width = 1.4*210, height = 1.4*297, units = "mm")





