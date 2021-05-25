############### STEP 3: Visual Inspection of Time Courses and Topography ###############
# This script contains some convenience functions to make the correlation matrices easier to use
# Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
# Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig

sort_similarities <- function(similarities){
  # Report maximum correlation per row to see if there is 
  # a match for a factor.
  similarities_sorted <- t(apply(round(similarities,2), 1, sort, decreasing = TRUE))
  
  # Which child factor is most similar to each adult factor?
  whichFactor <- t(apply(similarities, 1, order, decreasing = TRUE))
  
  out <- matrix(paste0("chF", whichFactor,"_r:", similarities_sorted), nrow = nrow(similarities),
         ncol = ncol(similarities))
  
  rownames(out) <- rownames(similarities)
  
  out
}


