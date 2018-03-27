####################################################################################
############################ KALMAN FILTER !! #######################################
####################################################################################
# We create the data matrix by combining the positions and velocities in a single 
# matrix !! We will need only to process these values, because for the prediction
# of the last 50 samples, we do not use them to update the model.
vel_mlat_tr = c(0,diff_mlat_tr)  # Set initial velocities to 0
vel_mlong_tr = c(0,diff_mlong_tr)

aux = rbind(c(mlat_tr), 
            c(vel_mlat_tr), 
            c(mlong_tr), 
            c(vel_mlong_tr))

Matrix_tr = t(matrix(aux, nrow = 4))  # Nsam * 4

# System of equations A
#   Xt = A*Xt-1 + B*Ut + SystemNoise  
#   Yt = Cx_k  + MeasurementNoise
# In our case B*Ut is 0, we do not have acceleration modeled.
########################### AR matrix ###############################
# We asume that x and y coordinates are independent 
# and that the dependence between variables is:
# xt = x_(t-1) + vx_(t-1)
# vxt = vx_(t-1)
# This leads to the next Autocorrelation Matrix A
A <- matrix(rbind(c(1,1,0,0),
                  c(0,1,0,0),
                  c(0,0,1,1),
                  c(0,0,0,1)),
            nrow = 4)

######################### C Matrix ###################################
# Measurement matrix.
# What transformation of the real state we perfrom.
# In this case we just select the position variables xt and yt
C <- matrix(rbind(c(1,0,0,0),
                  c(0,0,1,0)),
            nrow = 2)

# Dynamic System error Covariance Matrix
# The error that ocurrs in the system due to variables that we are not taking
# into account in the model (external forces: wind, vibration)
varPos = 3*1e-4  # Variance of the position noise
varVel = 1e-5    # Variance of the velocity noise 
Sigma1 <- matrix(rbind( c(varPos,0,0,0),
                        c(0,varVel,0,0),
                        c(0,0,varPos,0),
                        c(0,0,0,varVel)),
                 nrow = 4)

# Measurement error Covariance Matrix
# The covariance matrix of the measurement error observed at a given time t.
# We are only measuring position and we believe the measurements are uncorrelated
varNoise = 1e-4
Sigma2 <- matrix(rbind(c(varNoise,0),
                       c(0,varNoise)),
                 nrow = 2)

####################### Initialize Kalman Filter ##########################
## Time to choose our priors of prediction !!
## Important: 
#      - Variables with "hat" are estimated variables (V(t|t))
#      - Variables with "pred" are predicted variables (V(t+1|t))
# In the begining we have to define Xpred(1|0).
# Observed samples start at time t = 1. t = 0 is our prior thingy.

## We just initialize the parameters of the kalman it with a good guess. 
# In this case we just use the initial obervation of Yo as mean of the prior. X0hat = Y0
Xpred <- matrix(Matrix_in[1,])  # Initilization Xpred for X1 given Y0

# Initilize the uncertainty of the first samples high enough if we are not sure
# or low enough if we are sure of out initial Xhat
# As more samples arise, this influence of this initial decision will be decrease.
# We choose a small uncertainty of the initializatoin
sigmapred0 = 1e-5
SigmaXXpred<- matrix(rbind( c(sigmapred0,0,0,0),   # Initilization SigmaXXpred for X1 given Y0
                            c(0,sigmapred0,0,0),
                            c(0,0,sigmapred0,0),
                            c(0,0,0,sigmapred0)),
                     nrow = 4)

## Initial uncertainty of the observed variables. We say, small.
# Maybe same as the measurement noise would be appropiate.
SigmaYYpred <-  C%*%SigmaXXpred%*%t(C) + Sigma2
## Calculate covariance matrix between observarionts and latent variables
SigmaXYpred <- SigmaXXpred %*% t(C)

######### Initialize variables just to keep intermediate results ##########
XhatList <- matrix(0,nrow = 4, ncol = Ntr) # Store here the estimations of X
SigmaXXhatList <- list()        # Store here the estimated Sigma of X

XpredList <- matrix(0,nrow = 4, ncol = Ntr) # Store here the predictions of X
SigmaXXpredList <- list()        # Store here the predicted Sigma of X

###########################################################################
################## RECONSTRUCTION OF THE STATE ############################
###########################################################################
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
  SigmaYYpred <- C%*%SigmaXXpred%*%t(C) + Sigma2
  SigmaXYpred <- SigmaXXpred%*%t(C)
  
  ################# Storing data  ################
  XhatList[,n] <-Xhat
  SigmaXXhatList[[n]] <- SigmaXXhat
  
  XpredList[,n] <-Xpred
  SigmaXXpredList[[n]] <- SigmaXXpred
}
###########################################################################
################## PREDICTION OF THE STATE ############################
###########################################################################

# Now, using the final state of training (last estimations), we will predict 
# the state and observations for the last 50 samples (test samples)

## The initial hat is the one predicted for the last training sample
# We initilize with last parameters calculated
Xhattest = matrix(XhatList[,Ntr], ncol = 1) 
SigmaXXhattest <- SigmaXXhatList[[Ntr]]
Xpredtest = matrix(XpredList[,Ntr], ncol = 1) 
SigmaXXpredtest <- SigmaXXpredList[[Ntr]]

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
