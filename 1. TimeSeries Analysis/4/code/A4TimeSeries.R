
source("./graphlib2.R")
source("./mystatlib.R")
source("./step.slow.marima.R")

################ Question 4-1 ######################
## Loading and preparing the data:
myData = read.csv("./A4_gps_log.csv", sep = ",", header = TRUE) # read the csv
t = as.matrix(myData[,1])   # Get the date
mlat = as.matrix(myData[,2])   # Get the hour
mlong = as.matrix(myData[,3])     # Get the NO concentration 

Nsam = dim(t)[1]             # Number of samples

### Split data into train and test
Ntst = 50
tr_indx = 1:(Nsam - Ntst)  # Indexes for estimating
Ntr = length(tr_indx)
tst_indx = (Ntr + 1):Nsam   # Indexes for estimating
Ntst = length(tst_indx)

# Obtain data
mlat_tr = mlat[tr_indx]
mlat_tst = mlat[tst_indx]

mlong_tr = mlong[tr_indx]
mlong_tst = mlong[tst_indx]

t_tr = t[tr_indx]
t_tst = t[tst_indx]

## PLOT BOTH SIGNALS IN TERMS OF TIME

plot_timeSeries(t_tr, mlong_tr, "Longitude and Latitude", "Time", "Space","l", 2, "red", NULL, c(-15,10),0)
lines(t_tr, mlat_tr, col="blue", lwd = 2)
legend("topleft",c("Longitude","Latitude"),lwd=c(2,2),col=c("red","blue"))

plot_timeSeries(mlat, mlong, "Spatial aspects", "Latitude (x)", "Longitude (y)", "l",2, "black", NULL, NULL, 0)


############## Transformations of the data !##############
diff_mlat_tr = diff(mlat_tr)
diff_mlong_tr = diff(mlong_tr)

diffdiff_mlat_tr = diff(diff_mlat_tr)
diffdiff_mlong_tr = diff(diff_mlong_tr)

############# Create Train Matrixes ###############
aux = rbind(c(mlat_tr), c(mlong_tr))
train_Matrix = t(matrix(aux, nrow = 2))

aux = rbind(c(diff_mlat_tr), c(diff_mlong_tr))
difftrain_Matrix = t(matrix(aux, nrow = 2)) 

aux = rbind(c(diff_mlat_tr), c(diff_mlong_tr))
diffdifftrain_Matrix = t(matrix(aux, nrow = 2)) 

##############################################
######  PLOTTING
par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(2,0.7,0))

# Velocity Properties !!
plot_timeSeries(t_tr[-Ntr], diff_mlong_tr, "mlong", "Time", 
                "Increase in longitude","l",1,"red",NULL,NULL,0)
plot_timeSeries(t_tr[-Ntr], diff_mlat_tr, "mlat", "Time", 
                "Increase in latitude","l",1,"blue",NULL,NULL,0)

par(mfrow=c(1,1))
# Joint changes
plot_timeSeries(diff_mlong_tr, diff_mlat_tr, "Correlation mlong mlat", "mlong", 
                "mlat","p",1,"black",NULL,NULL,0)
plot_timeSeries(mlong_tr, mlat_tr, "Correlation mlong mlat", "mlong", 
                "mlat","p",1,"black",NULL,NULL,0)
# Acceleration Poperties !!
# Velocity Properties !!
plot_timeSeries(t_tr[1:(Ntr-2)], diffdiff_mlong_tr, "mlong", "Time", 
                "Increase in longitud","l",1,"red",NULL,NULL,0)
plot_timeSeries(t_tr[1:(Ntr-2)], diffdiff_mlat_tr, "mlat", "Time", 
                "Increase in latitud","l",1,"blue",NULL,NULL,0)
# Joint changes
plot_timeSeries(diffdiff_mlong_tr, diffdiff_mlat_tr, "mlat", "mlong", 
                "Correlation mlong mlat","p",1,"blue",NULL,NULL,0)

##############################################
####### View Basic Stufff ##################
summary(myData[c("mlat","mlong")])  # Same stuff
summary(t(matrix(mlat, mlong)))

################ Question 2 ######################
# We just use the ACF and PACF functions to calculate the coefficients.
par(c(1,1))
## Original data !!
Nlag = 50
plot_acfpacf(mlat_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_acfpacf(mlong_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_ccfpccf(mlat_tr,mlong_tr, "NOtr_ACFPACF.png",Nlag, 0)

## Differenced data !
plot_acfpacf(diff_mlat_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_acfpacf(diff_mlong_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_ccfpccf(diff_mlat_tr,diff_mlong_tr, "NOtr_ACFPACF.png",Nlag, 0)

## Differenced data acceleration!
plot_acfpacf(diffdiff_mlat_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_acfpacf(diffdiff_mlong_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_ccfpccf(diffdiff_mlat_tr,diffdiff_mlong_tr, "NOtr_ACFPACF.png",Nlag, 0)

#############################################################################
################ Question 3 ################################################
#############################################################################
## Starting marima
library(marima)
##################### TRY DIFFERENT MODELS #####################################
## Let us start with a marima(1,0,1)
par(mfrow=c(1,1))

ar = c(1,2)
ma = c(1,2)
struct11 <- define.model(kvar=2, ar=ar, ma=ma, indep=NULL) # rem.var is to ignore the years
M1 <- marima(difftrain_Matrix, means=1, ar.pattern=struct11$ar.pattern,
             ma.pattern=struct11$ma.pattern, Check=FALSE, Plot="log.det", penalty=0)
M1
# print estimates
short.form(M1$ar.estimates, leading=FALSE) 
short.form(M1$ma.estimates, leading=FALSE)

Mfinal = step.slow (M1, difftrain_Matrix, penalty=2, max.iter=50)
short.form(Mfinal$ar.estimates, leading=FALSE) 
short.form(Mfinal$ma.estimates, leading=FALSE)


########## OTHER ADVANCES STUFF #############
## Using regression variables
Model3 <- define.model(kvar=7, ar=ar, ma=ma, rem.var=c(1), reg.var=c(6,7))
Marima3 <- marima(all.data,means=1, ar.pattern=Model3$ar.pattern,
                  ma.pattern=Model3$ma.pattern, Check=FALSE, Plot='log.det', penalty=0)
short.form(Marima3$ar.estimates, leading=FALSE) # print estimates
short.form(Marima3$ma.estimates, leading=FALSE)

## Using a penalty
Model4 <- define.model(kvar=7, ar=ar, ma=ma, rem.var=1, reg.var=c(6,7))
Marima4 <- marima(all.data[1:90,], means=1, ar.pattern=Model4$ar.pattern,
                  ma.pattern=Model4$ma.pattern, Check=FALSE, Plot='log.det', penalty=1)
short.form(Marima4$ar.estimates, leading=FALSE) # print estimates
short.form(Marima4$ma.estimates, leading=FALSE)
Marima4$ar.fvalues
Marima4$ma.fvalues


####################################################################################
############################ KALMAN SHIT !! ########################################
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
#   Xt = A*Xt-1 + B*Ut + SystemNoise  # In our case B*Ut is 0, we do not have acceleration modeled.
#   Yt = Cx_k  + MeasurementNoise

########################### AR matrix ###############################
# We asume that x and y coordinates are independent and that the dependence between variables is:
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
varPos = 3*10e-4  # Variance of the position noise
varVel = 10e-5    # Variance of the velocity noise 
Sigma1 <- matrix(rbind( c(varPos,0,0,0),
                        c(0,varVel,0,0),
                        c(0,0,varPos,0),
                        c(0,0,0,varVel)),
                        nrow = 4)

# Measurement error Covariance Matrix
# The covariance matrix of the measurement error observed at a given time t.
# We are only measuring position and we believe the measurements are uncorrelated
varNoise = 10e-1
Sigma2 <- matrix(rbind(c(varNoise,0),
                       c(0,varNoise)),
                      nrow = 2)

####################### Initialize Kalman Filter ##########################
## Important: 
#      - Variables with "hat" are estimated variables (V(t|t))
#      - Variables with "pred" are predicted variables (V(t+1|t))
# In the begining we have to define Xpred(1|0).
# Observed samples start at time t = 1. t = 0 is our prior thingy.

## We just initialize the parameters of the kalman it with a good guess. 
# In this case we just use the initial obervation of Yo as mean of the prior. X0hat = Y0
Xpred <- matrix(Matrix_tr[1,])  # Initilization Xpred for X1 given Y0

# Initilize the uncertainty of the first samples high enough if we are not sure
# or low enough if we are sure of out initial Xhat
# As more samples arise, this influence of this initial decision will be decrease.
# We choose a small uncertainty of the initializatoin
sigmapred0 = 10e-10
SigmaXXpred<- matrix(rbind( c(sigmapred0,0,0,0),   # Initilization SigmaXXpred for X1 given Y0
                           c(0,sigmapred0,0,0),
                           c(0,0,sigmapred0,0),
                           c(0,0,0,sigmapred0)),
                    nrow = 4)

## Initial uncertainty of the observed variables. We say, small.
# Maybe same as the measurement noise would be appropiate.
sigmapred0 = 10e-10
SigmaYYpred <- matrix(rbind(c(sigmapred0,0),  # Initilization SigmaYYpred for X1 given Y0
                            c(0,sigmapred0)),
                           nrow = 2)

## Calculate covariance matrix between observarionts and latent variables
SigmaXYpred <- SigmaXXpred %*% t(C)

######### Initialize variables just to keep intermediate results ##########
XhatList <- matrix(0,nrow = 4, ncol = Ntr) # Store here the estimations of X
SigmaXXhatList <- list()        # Store here the estimated Sigma of X

XpredList <- matrix(0,nrow = 4, ncol = Ntr) # Store here the predictions of X
SigmaXXpredList <- list()        # Store here the predicted Sigma of X

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


######################################################################
############ Compare the Y estimated with the measurements !!! #######
######################################################################
Yhat =   t(C %*% XhatList)
train_Matrix
diff = Yhat - train_Matrix

plot_timeSeries(t_tr[5:50], mlong_tr[5:50], "Estimations","space","time","l",2, "black", NULL, NULL, 0)
lines(t_tr[5:50], Yhat[5:50,2], col="green", lwd = 5)
legend("topleft",c("Estimations","Observations"),lwd=c(2,2),col=c("black","green"))


plot(Yhat, train_Matrix)

##################  Predictions  #########################################v
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

######################################################################
############ Compare the Y predicted with the measurements !!! #######
######################################################################
Ypredtest =   t(C %*% XpredtestList)
difftest = Ypredtest - test_Matrix
YpredSigmaYY = C %*% SigmaXXpredtestList[[1]] %*% t(C)

plot(t_tst, Ypredtest[,1])
plot(t_tst,mlong_tst)

plot(t_tst, Ypredtest[,2])
lines(t_tst,mlat_tst)

## Plot some elipses for the noise ?
library(mixtools)
library(mvtnorm) 

plot(t_tst, Ypredtest[,2])
lines(t_tst,mlong_tst)

plot(t_tst,mlat_tst, ylim = c(6.2,6.7))
lines(t_tst, Ypredtest[,1])

plot(Ypredtest[,1], Ypredtest[,2], xlim = c(6.2,6.7), ylim = c(-13,-11))
lines(mlat_tst, mlong_tst)


for (i in 1:5){
  j = i * 9
  ellipse(mu=Ypredtest[j,], sigma= 0.003*C %*% SigmaXXpredtestList[[j]]  %*% t(C), alpha = .05, npoints = 250, col="red") 
}






