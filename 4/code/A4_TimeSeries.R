source("./graphlib.R")
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
plot_timeSeries(t_tr, mlong_tr, "NO concentration", "Time", "Concentration", 2, 0)
lines(t_tr, mlat_tr, lwd = 5)

plot_timeSeries(mlat_tr, mlong_tr, "Title", "latitud (x)", "longitud (y)", 2, 0)

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

############# Create Test Matrixes ###############
aux = rbind(c(mlat_tst), c(mlong_tst))
test_Matrix = t(matrix(aux, nrow = 2))

############# Create Both Matrixes ###############
# Because of the stupid way the MARMA forecast works
aux = rbind(c(diff(mlat)), c(diff(mlong)))
dataMatrix = t(matrix(aux, nrow = 2))

##############################################
######  PLOTTING
par(mfrow=c(1,1))

# Velocity Properties !!
plot_timeSeries(t_tr[-Ntr], diff_mlong_tr, "mlong", "Time", 
                "Increase in longitud",2,0)
plot_timeSeries(t_tr[-Ntr], diff_mlat_tr, "mlat", "Time", 
                "Increase in latitud",2,0)

# Joint changes
plot_timeSeries(diff_mlong_tr, diff_mlat_tr, "mlat", "mlong", 
                "Correlation mlong mlat",2,0)

# Acceleration Poperties !!
# Velocity Properties !!
plot_timeSeries(t_tr[1:(Ntr-2)], diffdiff_mlong_tr, "mlong", "Time", 
                "Increase in longitud",2,0)
plot_timeSeries(t_tr[1:(Ntr-2)], diffdiff_mlat_tr, "mlat", "Time", 
                "Increase in latitud",2,0)
# Joint changes
plot_timeSeries(diffdiff_mlong_tr, diffdiff_mlat_tr, "mlat", "mlong", 
                "Correlation mlong mlat",2,0)

##############################################
####### View Basic Stufff ##################
summary(myData[c("mlat","mlong")])  # Same stuff
summary(t(matrix(mlat, mlong)))

################ Question 2 ######################
# We just use the ACF and PACF functions to calculate the coefficients.
par(c(1,2))
## Original data !!
Nlag = 50
source("./graphlib.R")
plot_acfpacf(mlat_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_acfpacf(mlong_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_ccfpccf(mlat_tr,mlong_tr, "NOtr_ACFPACF.png",Nlag, 0)

## Differenced data !
plot_acfpacf(diff_mlat_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_acfpacf(diff_mlong_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_ccfpccf(diff_mlat_tr,diff_mlong_tr, "NOtr_ACFPACF.png",Nlag, 0)
##mierda = acf(cbind(diff_mlat_tr,diff_mlong_tr), lag.max = Nlag)

## Differenced data acceleration!
plot_acfpacf(diffdiff_mlat_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_acfpacf(diffdiff_mlong_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_ccfpccf(diffdiff_mlat_tr,diffdiff_mlong_tr, "NOtr_ACFPACF.png",Nlag, 0)

library(tseries)
adf.test(diff_mlat_tr, alternative = "stationary",k = 24)
adf.test(diff_mlong_tr, alternative = "stationary",k = 24)
#############################################################################
################ Question 3 ################################################
#############################################################################
## Starting marima
library(marima)
##################### TRY DIFFERENT MODELS #####################################
## Let us start with a marima(1,0,1)
par(mfrow=c(1,1))

ar = c(1,2,3)
ma = c(1,1)
max.iter = 100
struct11 <- define.model(kvar=2, ar=ar, ma=ma, indep=NULL) # rem.var is to ignore the years
M1 <- marima(difftrain_Matrix, means=1, ar.pattern=struct11$ar.pattern, max.iter = max.iter,
             ma.pattern=struct11$ma.pattern, Check=FALSE, Plot="log.det", penalty=0)
M1
# print estimates
short.form(M1$ar.estimates, leading=FALSE) 
short.form(M1$ar.pvalues, leading=FALSE) 

short.form(M1$ma.estimates, leading=FALSE)
short.form(M1$ma.pvalues, leading=FALSE)

Mfinal = step.slow (M1, difftrain_Matrix, penalty=2, max.iter=100)
short.form(Mfinal$ar.estimates, leading=FALSE) 
short.form(Mfinal$ar.pvalues, leading=FALSE) 
short.form(Mfinal$ma.estimates, leading=FALSE)
short.form(Mfinal$ma.pvalues, leading=FALSE)

residuals = Mfinal$residuals[,4:(Ntr -1)]
residuals_mlat = residuals[1,]
residuals_mlong = residuals[2,]
plot_SigACFPACF(residuals_mlong, Nlag)
plot_SigACFPACF(residuals_mlat, Nlag)

plot_ccfpccf(residuals_mlat,residuals_mlong, "NOtr_ACFPACF.png",Nlag, 0)

## Error Sigma !!


SS <- sum(residuals_mlong^2)/Ntr  # Sum of square errors of the residual model 1
par(mfrow=c(1,2))

########################## Error Sigma ##############
# Histogram and draw gaussian
hist(residuals_mlong,probability=T,col='blue')
curve(dnorm(x,sd = sqrt(Mfinal$resid.cov[1,1])), col=2, lwd=2, add = TRUE)
# Do the samples quantities
qqnorm(residuals_mlong)
qqline(residuals_mlong)

########################## PERDIODOGRAM ##############
par(mfrow=c(1,2))
cpgram(diff_mlong_tr)
cpgram(residuals_mlong)

plot_ccfpccf(residuals_mlat,residuals_mlong, "NOtr_ACFPACF.png",Nlag, 0)
plot_ccfpccf(diff_mlat_tr,diff_mlong_tr, "NOtr_ACFPACF.png",Nlag, 0)

alltest(residuals_mlat)
### Check if the sum of squared acf values follow a chi-square distribution.
acfvals <- acf(residuals_mlat, type="correlation", plot=FALSE)$acf[1:15] # Get acf residuals
test.stat <- sum(acfvals^2) * (length(residuals)-1)
1 - pchisq(test.stat, length(acfvals)-6)
prob = pchisq(test.stat, length(acfvals)-6, lower.tail = FALSE) 

########################## VAR library tests ##############
var1l <- VAR(t(residuals), p=1)
summary(var1l)
vars.test(var1l)

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
################################## Prediction ############################
####################################################################################

# call the forecasting function using Marima and the prepared data:
Forecasts <-  arma.forecast(t(dataMatrix), nstart=949, nstep=50, marima=Mfinal)

##############################################################
###################### MLONG ###############################
##############################################################
### From here on the plot is constructed ###
## We start in 950, dataMatrix only has 999 samples corresponding to the las 999 samples, we
# cannot calculate the diff of the first sample.
# To generate the orginial signal for diff[i] we need to add the value of the signal signal[i-1]
# In this case, since the matrix diff is already desplazated to the future 1 positons, the arrays are
# already syncronized

########################### mlat ###########################

# For reconstruction x[i + 1] = x[i] + v[i]. Using the index yeah
# For the predictions, we are given v[i], v[i+1], v[i+2] ans so on,
# so to get the original signal we perfrom cumsum of these increments.
# Making gaussian assumption, the variance of the sum is equal to the sum
# of the variances 

Predict<- cumsum(Forecasts$forecasts[1,Ntr:(Nsam -1)]) + mlat[Ntr]
stdv<-sqrt(cumsum(Forecasts$pred.var[1,1,]))
upper.lim=Predict+stdv*1.96
lower.lim=Predict-stdv*1.96
par(mfrow=c(1,1))

## Plot the real signal !!! 
Nbefore = 10

## Check for correction :)
caca = dataMatrix[(Ntr - Nbefore - 1):(Nsam - 1),1] + mlat[(Ntr - Nbefore -1):(Nsam -1)]
caca2 = mlat[(Ntr - Nbefore):Nsam]
caca - caca2
plot(t[(Ntr - Nbefore):Nsam],caca ,type="l", xlab="time_index",
     ylab="mlat prediction", main="mlat prediction", col= "black",
     ylim = c(6.45, 6.7))

##lines(t[(Ntr - Nbefore):Nsam], caca2,type="l", col= "blue")

lines(t[(Ntr +1):(Nsam)], Predict,type="l", col= "blue")
lines(t[(Ntr+1 ):(Nsam)], upper.lim, type="l", lty = 2, col = "blue")
lines(t[(Ntr+1 ):(Nsam)], lower.lim, type="l", lty = 2, col = "blue")

########################### mlong ###########################
Predict<-cumsum(Forecasts$forecasts[2,Ntr:(Nsam -1)]) + mlong[Ntr]
stdv<-sqrt(cumsum(Forecasts$pred.var[2,2,]))
upper.lim=Predict+stdv*1.96
lower.lim=Predict-stdv*1.96
par(mfrow=c(1,1))

# plot results:
Nbefore = 10
plot(t[(Ntr - Nbefore):Nsam], dataMatrix[(Ntr - Nbefore - 1):(Nsam - 1),2] + mlong[(Ntr - Nbefore -1):(Nsam -1)],type="l", xlab="time_index",
     ylab="mlong prediction", main="mlong prediction", col= "black")

## Plot the real signal !!! 
lines(t[(Ntr +1):Nsam ], Predict,type="l", col= "blue")

lines(t[(Ntr +1):Nsam], upper.lim, type="l", lty = 2, col = "blue")
lines(t[(Ntr +1):Nsam], lower.lim, type="l", lty = 2, col = "blue")
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
sigmapred0 = 1e-4
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
######################################################################
############ Compare the Y estimated with the measurements !!! #######
######################################################################
Yhat =   t(C %*% XhatList)
Ypred = t(C %*% XpredList)
train_Matrix
diff = Yhat - train_Matrix

## Compute the variance:
mlat_stdXX = c()
mlat_stdYY = c()
mlong_stdXX = c()
mlong_stdYY = c()

for(n in 1:Ntr){
  mlat_stdXX = c(mlat_stdXX,  sqrt(SigmaXXhatList[[n]][1,1]))
  aux = C%*%SigmaXXhatList[[n]]%*%t(C) + Sigma2
  mlat_stdYY = c(mlat_stdYY, sqrt(aux[1,1])) 
  
  mlong_stdXX = c(mlong_stdXX,  sqrt(SigmaXXhatList[[n]][3,3]))
  aux = C%*%SigmaXXhatList[[n]]%*%t(C) + Sigma2
  mlong_stdYY = c(mlat_stdYY, sqrt(aux[2,2])) 
}


par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(2,0.7,0))
par(mfrow=c(1,1))
Nex = 50
plot_timeSeries(t_tr[1:Nex], mlat_tr[1:Nex], "Estimations", "t", "mlat_estimation", 2, 0)
lines(t_tr[1:Nex], Yhat[1:Nex,1], col = "blue", lwd = 2)
lines(t_tr[1:Nex], Yhat[1:Nex,1] + 1.96 * mlat_stdXX[1:Nex], col = "blue", lwd = 1,  lty = 2)
lines(t_tr[1:Nex], Yhat[1:Nex,1] - 1.96 * mlat_stdXX[1:Nex], col = "blue", lwd = 1,  lty = 2)
lines(t_tr[1:Nex], Yhat[1:Nex,1] + 1.96 * mlat_stdYY[1:Nex], col = "red", lwd = 1,  lty = 2)
lines(t_tr[1:Nex], Yhat[1:Nex,1] - 1.96 * mlat_stdYY[1:Nex], col = "red", lwd = 1,  lty = 2)

Nex = 50
plot(t_tr[1:Nex], mlong_tr[1:Nex],main = "Estimations", xlab = "t", 
     ylab = "mlong_estimation", type = "l")
lines(t_tr[1:Nex], Yhat[1:Nex,2], col = "blue", lwd = 2)
lines(t_tr[1:Nex], Yhat[1:Nex,2] + 1.96 * mlong_stdXX[1:Nex], col = "blue", lwd = 1,  lty = 2)
lines(t_tr[1:Nex], Yhat[1:Nex,2] - 1.96 * mlong_stdXX[1:Nex], col = "blue", lwd = 1,  lty = 2)
lines(t_tr[1:Nex], Yhat[1:Nex,2] + 1.96 * mlong_stdYY[1:Nex], col = "red", lwd = 1,  lty = 2)
lines(t_tr[1:Nex], Yhat[1:Nex,2] - 1.96 * mlong_stdYY[1:Nex], col = "red", lwd = 1,  lty = 2)


plot(Yhat, train_Matrix)

plot_timeSeries(t_tr[1:Nex], mlat_tr[1:Nex], "Estimations", "t", "mlat_estimation", 2, 0)
lines(t_tr[1:Nex], Ypred[1:Nex,1], col = "blue", lwd = 2)
plot(Ypred, train_Matrix)
plot(Ypred, Yhat)

Nlag = 30
Nlag2 = 30
residuals = t(Ypred[1:(Ntr-1),] - train_Matrix[2:(Ntr),]) # Right displacement

residuals_mlatMA = residuals_mlat
residuals_mlat = residuals[1,]
residuals_mlong = residuals[2,]
plot_SigACFPACF(residuals_mlong, Nlag)
tsdiag(fit1, gof.lag = Nlag2)  # You plot the residuals and their ACF and JunlgeBox

## Error Sigma !!

SS <- sum(residuals_mlat^2)/Ntr  # Sum of square errors of the residual model 1
par(mfrow=c(1,2))

########################## Error Sigma ##############
# Histogram and draw gaussian
hist(residuals_mlat,probability=T,col='blue')
curve(dnorm(x,sd = sqrt(Mfinal$resid.cov[1,1])), col=2, lwd=2, add = TRUE)
# Do the samples quantities
qqnorm(residuals_mlat)
qqline(residuals_mlat)

########################## PERDIODOGRAM ##############
par(mfrow=c(1,2))
cpgram(diff_mlat_tr)
cpgram(residuals_mlat)

alltest(residuals_mlat)
########################## VAR library tests ##############
var1l <- VAR(t(residuals), p=1)
summary(var1l)
vars.test(var1l)
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

######################################################################
############ Compare the Y predicted with the measurements !!! #######
######################################################################
Ypredtest =   t(C %*% XpredtestList)
difftest = Ypredtest - test_Matrix
YpredSigmaYY = C %*% SigmaXXpredtestList[[1]] %*% t(C)

## Compute the variance:
mlat_stdXXtest = c()
mlat_stdYYtest = c()
mlong_stdXXtest = c()
mlong_stdYYtest = c()

for(n in 1:Ntst){
  mlat_stdXXtest = c(mlat_stdXXtest,  sqrt(SigmaXXpredtestList[[n]][1,1]))
  aux = C%*%SigmaXXpredtestList[[n]]%*%t(C) + Sigma2
  mlat_stdYYtest = c(mlat_stdYYtest, sqrt(aux[1,1])) 
  
  mlong_stdXXtest = c(mlong_stdXXtest,  sqrt(SigmaXXpredtestList[[n]][3,3]))
  aux = C%*%SigmaXXpredtestList[[n]]%*%t(C) + Sigma2
  mlong_stdYYtest = c(mlong_stdYYtest, sqrt(aux[2,2])) 
}

par(mfrow=c(1,1))
Nbefore = 20
plot(t[(Ntr- Nbefore):Nsam], mlat[(Ntr- Nbefore):Nsam], type = "l",
     main = "Predictions", xlab = "t", ylab = "mlat_prediction", ylim  = c(5.0,8.0))

lines(t_tr[(Ntr-Nbefore):Ntr], mlat_tr[(Ntr-Nbefore):Ntr])
lines(t_tr[(Ntr-Nbefore):Ntr], Yhat[(Ntr-Nbefore):Ntr,1], col = "blue", lwd = 2)
lines(t_tr[(Ntr-Nbefore):Ntr], Yhat[(Ntr-Nbefore):Ntr,1] + 1.96 * mlat_stdXX[(Ntr-Nbefore):Ntr], col = "blue", lwd = 1,  lty = 2)
lines(t_tr[(Ntr-Nbefore):Ntr], Yhat[(Ntr-Nbefore):Ntr,1] - 1.96 * mlat_stdXX[(Ntr-Nbefore):Ntr], col = "blue", lwd = 1,  lty = 2)
lines(t_tr[(Ntr-Nbefore):Ntr], Yhat[(Ntr-Nbefore):Ntr,1] + 1.96 * mlat_stdYY[(Ntr-Nbefore):Ntr], col = "red", lwd = 1,  lty = 2)
lines(t_tr[(Ntr-Nbefore):Ntr], Yhat[(Ntr-Nbefore):Ntr,1] - 1.96 * mlat_stdYY[(Ntr-Nbefore):Ntr], col = "red", lwd = 1,  lty = 2)

lines(t_tst, Ypredtest[,1], col = "blue", lwd = 2, ylim  = c(6.0,8.0))
lines(t_tst, Ypredtest[,1] + 1.96 * mlat_stdXXtest, col = "blue", lwd = 1,  lty = 2)
lines(t_tst, Ypredtest[,1] - 1.96 * mlat_stdXXtest, col = "blue", lwd = 1,  lty = 2)
lines(t_tst, Ypredtest[,1] + 1.96 * mlat_stdYYtest, col = "red", lwd = 1,  lty = 2)
lines(t_tst, Ypredtest[,1] - 1.96 * mlat_stdYYtest, col = "red", lwd = 1,  lty = 2)

c(Ypredtest[1,1], Ypredtest[10,1], Ypredtest[25,1], Ypredtest[50,1])
c(Ypredtest[1,2], Ypredtest[10,2], Ypredtest[25,2], Ypredtest[50,2])

c(mlat_stdYYtest[1], mlat_stdYYtest[10], mlat_stdYYtest[25],mlat_stdYYtest[50])
c(mlong_stdYYtest[1], mlong_stdYYtest[10], mlong_stdYYtest[25],mlong_stdYYtest[50])

c(Ypredtest[1,1] + 1.96*mlat_stdYYtest[1], Ypredtest[10,1] + 1.96*mlat_stdYYtest[10], 
  Ypredtest[25,1] + 1.96*mlat_stdYYtest[25],Ypredtest[50,1] + 1.96*mlat_stdYYtest[50])
c(Ypredtest[1,1] - 1.96*mlat_stdYYtest[1], Ypredtest[10,1] - 1.96*mlat_stdYYtest[10], 
  Ypredtest[25,1] - 1.96*mlat_stdYYtest[25],Ypredtest[50,1] - 1.96*mlat_stdYYtest[50])

c(Ypredtest[1,2] + 1.96*mlong_stdYYtest[1], Ypredtest[10,2] + 1.96*mlong_stdYYtest[10], 
  Ypredtest[25,2] + 1.96*mlong_stdYYtest[25],Ypredtest[50,2] + 1.96*mlong_stdYYtest[50])

c(Ypredtest[1,2] - 1.96*mlong_stdYYtest[1], Ypredtest[10,2] - 1.96*mlong_stdYYtest[10], 
  Ypredtest[25,2] - 1.96*mlong_stdYYtest[25],Ypredtest[50,2] - 1.96*mlong_stdYYtest[50])

Nbefore = 20
plot(t[(Ntr- Nbefore):Nsam], mlong[(Ntr- Nbefore):Nsam], type = "l", ylim  = c(-13,-10),
     main = "Predictions", xlab = "t", ylab = "mlong_prediction")

lines(t_tr[(Ntr-Nbefore):Ntr], mlat_tr[(Ntr-Nbefore):Ntr])
lines(t_tr[(Ntr-Nbefore):Ntr], Yhat[(Ntr-Nbefore):Ntr,2], col = "blue", lwd = 2)
lines(t_tr[(Ntr-Nbefore):Ntr], Yhat[(Ntr-Nbefore):Ntr,2] + 1.96 * mlong_stdXX[(Ntr-Nbefore):Ntr], col = "blue", lwd = 1,  lty = 2)
lines(t_tr[(Ntr-Nbefore):Ntr], Yhat[(Ntr-Nbefore):Ntr,2] - 1.96 * mlong_stdXX[(Ntr-Nbefore):Ntr], col = "blue", lwd = 1,  lty = 2)
lines(t_tr[(Ntr-Nbefore):Ntr], Yhat[(Ntr-Nbefore):Ntr,2] + 1.96 * mlong_stdYY[(Ntr-Nbefore):Ntr], col = "red", lwd = 1,  lty = 2)
lines(t_tr[(Ntr-Nbefore):Ntr], Yhat[(Ntr-Nbefore):Ntr,2] - 1.96 * mlong_stdYY[(Ntr-Nbefore):Ntr], col = "red", lwd = 1,  lty = 2)

lines(t_tst, Ypredtest[,2], col = "blue", lwd = 2)
lines(t_tst, Ypredtest[,2] + 1.96 * mlong_stdXXtest, col = "blue", lwd = 1,  lty = 2)
lines(t_tst, Ypredtest[,2] - 1.96 * mlong_stdXXtest, col = "blue", lwd = 1,  lty = 2)
lines(t_tst, Ypredtest[,2] + 1.96 * mlong_stdYYtest, col = "red", lwd = 1,  lty = 2)
lines(t_tst, Ypredtest[,2] - 1.96 * mlong_stdYYtest, col = "red", lwd = 1,  lty = 2)

## Plot some elipses for the noise ?
library(mixtools)
library(mvtnorm) 

plot(t_tst, Ypredtest[,2])
lines(t_tst,mlong_tst)

plot(t_tst,mlat_tst, ylim = c(6.2,6.7))
lines(t_tst, Ypredtest[,1])

plot(Ypredtest[,1], Ypredtest[,2], xlim = c(5,8), ylim = c(-13,-11), 
     type = "l", xlab = "mlat", ylab = "mlong",  col = "blue")
lines(mlat_tst, mlong_tst, type = "l")


for (i in 1:5){
  j = (i-1) * 9 + 1
  ellipse(mu=Ypredtest[j,], sigma= 1*C %*% SigmaXXpredtestList[[j]]  %*% t(C), alpha = .05, npoints = 250, col="blue") 
  ellipse(mu=Ypredtest[j,], sigma= 1*C %*% SigmaXXpredtestList[[j]]  %*% t(C) + Sigma2, alpha = .05, npoints = 250, col="red") 
  }

########################################################################
####################### BONUS PART ########################
########################################################################

varPos = varPoses[1]
varVel = varVeles[4]
varNoise = varNoises[7]

source("./mystatlib.R") 

varPoses = exp(seq(log(1e-7),log(1e-5),length=10))
varVeles = exp(seq(log(1e-6),log(1e-5),length=10))
varNoises = exp(seq(log(1e-3),log(1e-1),length=10))

likies = array(0, dim=c(length(varPoses),length(varVeles),length(varNoises)))

for (i in 1:length(varPoses)){
  for (j in 1:length(varVeles)){
    for (k in 1:length(varNoises)){ 
      varPos = varPoses[i]
      varVel = varVeles[j]
      varNoise = varNoises[k]
      
Sigma1 <- matrix(rbind( c(varPos,0,0,0),
                        c(0,varVel,0,0),
                        c(0,0,varPos,0),
                        c(0,0,0,varVel)),
                 nrow = 4)

# Measurement error Covariance Matrix
# The covariance matrix of the measurement error observed at a given time t.
# We are only measuring position and we believe the measurements are uncorrelated

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
sigmapred0 = 10e-5
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

XpredList <- matrix(0,nrow = 4, ncol = Ntr + 1) # Store here the predictions of X
SigmaXXpredList <- list()        # Store here the predicted Sigma of X
SigmaYYpredList <- list() 
SigmaYYpredList[[1]] <- SigmaYYpred
XpredList[,0] = Xpred
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
  
  XpredList[,n+1] <-Xpred
  SigmaXXpredList[[n]] <- SigmaXXpred # We displace it 1 :) Necesary
  SigmaYYpredList[[n+1]] <- SigmaYYpred  # We displace it 1 :) Necesary
}
######################################################################
############ Compare the Y estimated with the measurements !!! #######
######################################################################
Ypred = t(C %*% XpredList)

like = get_likelihood(Matrix_tr %*% t(C), Ypred, SigmaYYpredList)
likies[i,j,k] = like

    }
  }
}


which(likies == max(likies), arr.ind = TRUE)


