step.slow <- function(object, data, penalty=2, max.iter=50){
  ## object: A marima object
  ## data:   The same data as given as argument when fitting 'object'
  ## penalty: The penalty to be used
  ## max.iter: Number of iterations before evaluating the penalty
  
  # Init
  
  obj <- object # So that the original object is returned if no reduction is needed.
  
  if (any(obj$ar.fvalues[,,-1] >0)){
    ar.f <- obj$ar.fvalues[,,-1]
    ar.p <- obj$out.ar.pattern[,,-1]
    ar.min <- min(ar.f[ar.f>0])
  } else{
    ar.min <- Inf
    ar.p <- obj$out.ar.pattern
  }
  if (any(obj$ma.fvalues[,,-1] >0)){
    ma.f <- obj$ma.fvalues[,,-1]
    ma.p <- obj$out.ma.pattern[,,-1]
    ma.min <- min(ma.f[ma.f>0])
  } else {
    ma.min <- Inf
    ma.p <- obj$out.ma.pattern
  }
  print(c(ar.min, ma.min))
  # Now starting the actual model reduction
  while (min(ar.min, ma.min) < penalty){
    if (ar.min < ma.min){
      ar.p[ar.f ==ar.min] <- FALSE
    } else{
      ma.p[ma.f ==ma.min] <- FALSE
      
    } # else
    ## Now restimate
    obj <- marima(data, ar.pattern = check.one(ar.p), 
                  ma.pattern = check.one(ma.p), 
                  max.iter = max.iter,
                  Plot='log.det')
    
    if (any(obj$ar.fvalues[,,-1] >0)){
      ar.f <- obj$ar.fvalues[,,-1]
      ar.p <- obj$out.ar.pattern[,,-1]
      ar.min <- min(ar.f[ar.f>0])
    } else{
      ar.p <- obj$out.ar.pattern
      ar.min <- Inf
    }
    if (any(obj$ma.fvalues[,,-1] >0)){
      ma.f <- obj$ma.fvalues[,,-1]
      ma.p <- obj$out.ma.pattern[,,-1]
      ma.min <- min(ma.f[ma.f>0])
    } else {
      ma.p <- obj$out.ma.pattern
      ma.min <- Inf
    }
    
  } # while
  return(obj)
}

