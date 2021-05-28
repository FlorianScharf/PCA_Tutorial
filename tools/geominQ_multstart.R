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
  
  
  if (rand.start){ # if you are supposed to use random starts
    
    # Use as many starting values as indicated
    mrRot <- lapply(1:start.values, function(iStart){
      # Generate a random starting rotation
      Tmat <- Random.Start(ncol(A))
      # Tell the user where you are in the process
      print(paste0("Random start number: ", iStart))
      
      # Apply gradient projection - parameters are passed from the function 
      rotFit <- GPFoblq(A = A, method = method, normalize = normalize, Tmat = Tmat,
                           eps = eps, maxit = maxit, methodArgs = list(delta = delta))
      # Remember the results in each round
      list(criterion = ifelse(is.null(rotFit$Table), NA, min(rotFit$Table[,2])), rotFit = rotFit)
    })
    
    # Which of the random starting values resulted in the best value
    # of the criterion?
    bestStart <- which.min(lapply(mrRot, function(x) x$criterion))
    
    # If bestStart is empty, notify the user about it.
    if (length(bestStart) == 0) stop("No best solution found, perhaps no rotation converged? Try increasing maxit or reducing eps!")
    # If the function received reasonable input, non-convergence is basically the only
    # thing that can go wrong. 
    
    # return the best rotation results as the result
    rotFit <- mrRot[[bestStart]]$rotFit
    
  }
  else { # else and rotate: start from a generic starting matrix
    # Apply gradient projection - parameters are passed from the function 
     rotFit = GPFoblq(A = A, method = method, normalize = normalize, Tmat = Tmat,
                        eps = eps, maxit = maxit, methodArgs = list(delta = delta))
  }
  
  return(rotFit)
}
