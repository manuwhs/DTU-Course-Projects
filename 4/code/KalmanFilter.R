
KalmanFilterEstimate <- function(Y,       # Observations in form Y[Nsamples, Ndim]
                                 A,B,C,  # Dynamics matrices
                                 Xpred1, SigmaXXpred1, SigmaYYpred1){  # Prior predictions
  Matrix_tr = Y
  Ntr = dim(Y)[0]
  Ndim = dim(Y)[1]
  
  Xpred = Xpred1
  SigmaXXpred = SigmaXXpred1
  SigmaYYpred = SigmaYYpred1
  
  for(n in 1:Ntr){
    #############   ESTIMATING STEP  ###################
    # Estimation of theparameters ("hat" variables V(t|t))
    # Calculated from the predicted values and the new sample
    K = SigmaXXpred %*% t(C) %*% solve(SigmaYYpred)  # Kalman Gain
    Xhat = Xpred + K %*% (C %*% Matrix_tr[n,] - C %*% Xpred)
    SigmaXXhat <- SigmaXXpred - K%*%SigmaYYpred%*%t(K)
    
    #############   PREDICTION STEP  ###################
    ## Predict the variables for the next instance  V(t+1|t)
    Xpred <- A%*%Xhat  # We use the correlation bet
    SigmaXXpred <- A%*%SigmaXXhat%*%t(A) + Sigma1
    SigmaYYpred <- C%*%SigmaXXhat%*%t(C) + Sigma2
    SigmaXYpred <- SigmaXXpred%*%t(C)
    
    ################# Storing data  ################
    XhatList[,n] <-Xhat
    SigmaXXhatList[[n]] <- SigmaXXhat
    
    XpredList[,n] <-Xpred
    SigmaXXpredList[[n]] <- SigmaXXpred
    
  }
  results = list(XhatList,SigmaXXpredList)
  # Return the estimated X and its estimated Covariance Matrix
}

KalmanFilterPredict <- function(Ntst,       # Number of predictions in advance
                                 A,B,C,  # Dynamics matrices
                                 Xpred1, SigmaXXpred1, SigmaYYpred1){  # Prior predictions

  Xpredtest = Xpred1
  SigmaXXpredtest = SigmaXXpred1
  SigmaYYpredtest = SigmaYYpred1
  
  ## Variables to store the results
  XpredtestList <- matrix(0,nrow = 4, ncol = Ntst)
  SigmaXXpredtestList <- list()
  
  # The first prediction is the one calculated last in the previous loop
  XpredtestList[,1] <- Xpredtest
  SigmaXXpredtestList[[1]] <- SigmaXXpredtest
  # We calculate the rest
  for(n in 2:(Ntst)){
    # Perform future predictions
    Xpredtest <- A%*%Xpredtest
    SigmaXXpredtest <- A%*%SigmaXXpredtest%*%t(A) + Sigma1 
    XpredtestList[,n] <-Xpredtest
    SigmaXXpredtestList[[n]] <- SigmaXXpredtest
  }
  
  results = list(XpredtestList,SigmaXXpredtestList)
  # Return the predicted X and its estimated Covariance Matrix
}





