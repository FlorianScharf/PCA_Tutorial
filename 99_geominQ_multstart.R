## Define rotation function for Geomin rotation with multiple starts
# Author: Florian Scharf, florian.scharf@uni-muenster.de and Andreas Widmann, widmann@uni-leipzig.de
# Copyright (c) 2021 Florian Scharf, University of Münster and Andreas Widmann, University of Leipzig
## DEPENDENCY: GPArotation
# If you want to use other rotation criteria, you can easily adjust this function by
# changing the value of the method argument into any rotation available here:
# ?GPArotation::GPFoblq
# For orthogonal solutions, you could replace GPFoblq with GPForth 
# We did not implement this here because we strongly recommend to use oblique rotations
# for ERP data as outlined in the article.

geominQ_multstart = function(A, normalize = FALSE, delta = 0.01, Tmat = diag(ncol(A)), start.values = 30, rand.start = TRUE,
                             maxit = 10000, eps = 1e-5, method = "geomin"){
  
  
  if (rand.start){
    
    mrRot <- lapply(1:start.values, function(iStart){
      Tmat <- Random.Start(ncol(A))
      print(paste0("Random start number: ", iStart))
      rotFit <- GPFoblq(A = A, method = method, normalize = normalize, Tmat = Tmat,
                           eps = eps, maxit = maxit, methodArgs = list(delta = delta))
      
      list(criterion = ifelse(is.null(rotFit$Table), NA, min(rotFit$Table[,2])), rotFit = rotFit)
    })
    
    bestStart <- which.min(lapply(mrRot, function(x) x$criterion))
    bestStart <- ifelse(length(bestStart) == 0, 1, bestStart)
    rotFit <- mrRot[[bestStart]]$rotFit
    
  }
  else {
    rotFit = GPFoblq(A = A, method = method, normalize = normalize, Tmat = Tmat,
                        eps = eps, maxit = maxit, methodArgs = list(delta = delta))
  }
  
  rotFit$Th <- t(solve(rotFit$Th))
  return(rotFit)
}
