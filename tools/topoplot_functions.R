
# Get Grand Averages to long format fpr plotting
gavr2long <- function(allGAVRs, iEl, iGroup, times){
  # Get grandaverages to plot 
  iGAVR <- allGAVRS[allGAVRS$chan == iEl & 
                      allGAVRS$group == iGroup, -c(4)]
  iGAVR$cond <- as.character(iGAVR$cond)
  iGAVR$chan <- as.character(iGAVR$cond)
  
  iGAVR[3,] <- c("nov - sta", iEl, iGroup, iGAVR[2,-c(1:3)] - iGAVR[1,-c(1:3)])
  
  # Reshape for plotting
  iGAVR <- melt(iGAVR, variable_name = "times")
  iGAVR$cond <- factor(iGAVR$cond, levels = c("sta", "nov", "nov - sta"))
  levels(iGAVR$times) <- times
  iGAVR$times <- as.numeric(as.character(iGAVR$times))
  
  return(iGAVR)
}


makeERPest  <- function(iScores, iPCA, iFactor, iEl, times){
  iERP_est <- t(apply(iScores, MARGIN = 1, function(x){
    iPCA$rotFit$loadings[, iFactor] * as.numeric(x[iFactor])
  }))
  
  iERP_est <- cbind(iScores[,1:2], iERP_est)
  
  iEst <- iERP_est[iERP_est$chan == iEl,]  
  iEst$cond <- as.character(iEst$cond)
  iEst$chan <- as.character(iEst$cond)
  iEst[3,] <- c("nov - sta", iEl, iEst[2,-c(1:2)] - iEst[1,-c(1:2)])
  
  # Reshape for plotting
  iEst <- melt(iEst, variable_name = "times")
  iEst$cond <- factor(iEst$cond, levels = c("sta", "nov", "nov - sta"))
  levels(iEst$times) <- times
  iEst$times <- as.numeric(as.character(iEst$times))
  iEst
  
}


pca2eeg <- function(iPCA, allAVRs, iGroup, iFactor){
  signals <- lapply(c("sta", "nov"), function(iCond) {
    iSignals <- matrix(iPCA$average_scores[iPCA$average_scores$cond == iCond, iFactor] * max(iPCA$rotFit$loadings[,iFactor]),
                       ncol = dim(allAVRs[[iGroup]][[iCond]]$signals)[2], 
                       nrow = dim(allAVRs[[iGroup]][[iCond]]$signals)[1], 
                       byrow = TRUE)
    
    
    # We need to convert this to the tibble_format used in eegUtils:
    colnames(iSignals) <- names(allAVRs[[iGroup]][[iCond]]$signals)
    
    iSignals
  })
  names(signals) <- c("sta", "nov")
  
  signals[["nov - sta"]] <- signals$nov - signals$sta
  
  signals
}
  