a = as.matrix(cbind(c(1,0),c(1,1)))
B = matrix(c(2, 4, 3, 1, 5, 7),  nrow=3, ncol=2)
# Power of a matrix function
powA = function(a,n) {
  if (n < 0){
    n = -n
    a = solve(a)
  }
  if (n==1)  return (a)
  if (n==2)  return (a%*%a)
  if (n>2) return ( a%*%powA(a,n-1))
}

################ Question 1 ######################
## Loading and preparing the data:
myData = read.csv("./apartment_prices.csv", sep = " ") # read the csv
time = as.matrix(myData[,2])   # Get the time
price = as.matrix(myData[,3])  # Get the price
Nsam = dim(time)[1]             # Number of samples
# Vars for plotting
minPrice = min(price)
maxPrice = max(price)
rangePrice = maxPrice - minPrice
### Split data into train and test
tr_indx = as.matrix(which(time[,1] < 2015))  # Indexes for estimating
Ntr = dim(tr_indx)[1]

tst_indx = (Ntr+1):Nsam
Ntst = Nsam - Ntr

# Obtain data
time_tr = as.matrix(time[tr_indx,1])
time_tst = as.matrix(time[tst_indx,1])
price_tr = as.matrix(price[tr_indx,1])
price_tst = as.matrix(price[tst_indx,1])

## Do the plotting and saving it into an image
# Save the image in this file, with res =
png(file="Price_time.png",width=800, res=90) 
# Size and aspect ratio could be width=400,height=350,res=45
plot(time, price, 
     type = "l",                     # Draw as a line
     lwd= 3,                         # Line width
     main="Apartment prices", # Title of the graph
     xlab="Time (years)",            # x label
     ylab="Price (Kr.)")                  # y label
dev.off()

################ Question 2 ######################
# Obtain the mean and the std for estimation
mu = mean(price_tr)               # Calculate mean
Ytr = matrix(mu, Ntr)            # Output 
error_tr = Ytr - price_tr      # Obtain error estimate
# Obtain theoretical error std
std_error_tr = sqrt(t(error_tr)%*%error_tr / (Ntr - 1))  


# Obtain prediction and empirical std of prediction
Ytst = matrix(mu, Ntst)          # Obtain output prediction
error_tst = Ytst - price_tst      # Obtain error estimate
error_tst = Ytst - price_tst      # Obtain error estimate

# Obtain empirical error std
std_error_tst = sqrt(t(error_tst)%*%error_tst / (Ntst - 1))  

## Do the plotting and saving it into an image
png(file="Mean_Model.png",width=800, res=90)  # Save the image in this file, with res =
plot(time, price, type = "l",lwd= 3,main="CMM model for price estimation", 
     xlab="Time (years)", ylab="Price ",
     ylim = c(minPrice - rangePrice/2, maxPrice + rangePrice/3))

## Plot mean +- 2std lines
## Train prediction and std of tr error
lines(time_tr, Ytr,col = "blue", lwd=2) # Extend the values to create vector and plot
lines(time_tr, Ytr + 2 * matrix(std_error_tr, Ntr), col = "blue", type = "l", lty=2, lwd=1)
lines(time_tr, Ytr - 2 * matrix(std_error_tr, Ntr), col = "blue", type = "l", lty=2, lwd=1)

## Test prediction and std of tst empirical error 
lines(time_tst, Ytst,col = "red", lwd=2) # Extend the values to create vector and plot
lines(time_tst, Ytst + 2 * matrix(std_error_tst, Ntst), col = "red", type = "l", lty=2, lwd=1)
lines(time_tst, Ytst - 2 * matrix(std_error_tst, Ntst), col = "red", type = "l", lty=2, lwd=1)
dev.off()

################ Question 3 ######################
# Use a linear solver as y = w0 + w1 * t + e to predict the samples
Xtr <- cbind(matrix(1,Ntr),time_tr)   # Add a vector of 1s for the bias
thetahat <- solve(t(Xtr) %*% Xtr) %*% t(Xtr) %*% price_tr  # Calculate regression coefficients

Ytr = Xtr %*% thetahat        # Obtain output prediction
error_tr = Ytr - price_tr      # Obtain error estimate

# Obtain the unbiased estimator of tr error (estimation error). error_tr = [92,1]
NdimParam = dim(thetahat)[1]
std_error_tr = sqrt(t(error_tr)%*%error_tr / (Ntr - NdimParam))  # Obtain error std

# Obtain test estimates (prediction estimates):
Xtst = cbind(matrix(1,Ntst),time_tst)
Ytst = Xtst %*% thetahat        # Obtain output prediction
error_tst = Ytst - price_tst      # Obtain error estimate
error_tst = Ytst - price_tst      # Obtain error estimate

std_error_tst = sd(error_tst, na.rm = FALSE)   # Obtain empirical error std

# Plot the predicted and the real data.
## Do the plotting and saving it into an image
png(file="Linear_Model.png",width=800, res=90)  # Save the image in this file, with res =
plot(time, price, type = "l",lwd= 3,main=" SLM model for price estimation", xlab="Time (years)", 
     ylab="Price ", ylim = c(minPrice - rangePrice/5, maxPrice + rangePrice/5)) 

## Plot mean +- 2std lines
lines(time_tr, Ytr,col = "blue", lwd=2) # Extend the values to create vector and plot
lines(time_tr, Ytr + 2 * matrix(std_error_tr, Ntr), col = "blue", type = "l", lty=2, lwd=1)
lines(time_tr, Ytr - 2 * matrix(std_error_tr, Ntr), col = "blue", type = "l", lty=2, lwd=1)
## Test prediction and std of tst empirical error 
lines(time_tst, Ytst,col = "blue", lwd=2) # Extend the values to create vector and plot
lines(time_tst, Ytst + 2 * matrix(std_error_tst, Ntst), col = "red", type = "l", lty=2, lwd=1)
lines(time_tst, Ytst - 2 * matrix(std_error_tst, Ntst), col = "red", type = "l", lty=2, lwd=1)
dev.off()

################ Question 4 ######################

# Exponential smoothing: We predict the value in instant N + 1 from a
# weighted average of the previous samples. The weights follow a exponential
# decreasing.
# Y^_N+l+1|N+1 = (1-lambda)Y_N+1 + Y^_N+l|N 

lambda = 0.8
priceHatTr = matrix(price[1],1,1) #Initial value is the sample itself
for(i in 1:Ntr){
  yhat_i = (1-lambda)*price_tr[i]+lambda*priceHatTr[i] 
  priceHatTr = cbind(priceHatTr, yhat_i)
}

# Predictions for data Test
last_Tr = price[Ntr]    # Last known price 
priceHatTst = matrix(priceHatTr[Ntr],1,1)  # Price estimate for tst.
# First value will be removed after
for(i in 1:Ntst){
  pred = (1-lambda)*last_Tr+lambda*priceHatTst[i] 
  priceHatTst = cbind(priceHatTst, pred)
}

## Plot the initial

png(file="Exp_Smoothing.png",width=800, res=90)  # Save the image in this file, with res =
plot(time, price, type = "l",lwd= 3,main="Exponential Smoothing", xlab="Time (years)", ylab="Price ",
     ylim = c(min(price)*(1- 1/10), max(price) + 200)*(1+ 1/10))   
lines(time_tr, priceHatTr[1:Ntr], col="blue",lwd=2)

## Plot the predicted data !!
lines(time_tst, priceHatTst[-1],col="red",lwd=2)
dev.off()



##################################################################
################ Question 5.1 Unefficient way #####################
##################################################################

## Brute force, no updating.
# Calculate the estimation error in the training set
# Create "imaginary" time instances for ALL time instances
X_j = cbind(matrix(1,Nsam), matrix(-seq(Nsam-1,0,-1)))
lambda = 0.47

# Transition matrix
L = matrix(1,2,2) 
L[1,2] = 0
# f(0): Imaginary time in instant 0
f0 = matrix(1,2,1) #f(0)
f0[2] = 0
fj = L %*% f0

## ESTIMATION (GET Initial parameters)
## Create weight matrix
Ninit = 3
Nest = Ntr - Ninit   # Number of samples to use for updating
# Main variables where to store the predictions and theoretical std at every point.
Ypred = matrix(1,Nest,1)
Ypred_std = matrix(1,Nest,1)
for (Nini in (Ninit:Ntr)){
  X_j = cbind(matrix(1,Nini), matrix(-seq(Nini-1,0,-1)))
  CovInv = matrix(0,Nini,Nini)  # Inverse covariance Matrix of the samples. (ML model)
  for(i in 0:(Nini-1)){
    CovInv[Nini-i,Nini-i] = lambda^(i)
  }
  # Solve the system to obtain the GLS parameters.
  # N * Covariance of the dimensions of f(x)_j = [1 j] (Imaginary time)
  FMatrix = t(X_j[1:Nini,]) %*% CovInv %*% X_j[1:Nini,]
  H = t(X_j[1:Nini,]) %*% CovInv %*% price_tr[1:Nini,]      # Samples * Weight_Matrix * Y
  thetahat = solve(FMatrix) %*% H  # Same as GLS
  Ytr = X_j %*% thetahat
  
  # Get training variance
  erorr_tr = price_tr[1:Nini] - Ytr
  Sigma2 = (t((erorr_tr)) %*%CovInv %*% (erorr_tr)) 
  Sigma2 = Sigma2 / (Nini - Ndim)
  
  # Get next test and its std
  Ytst = t(powA(L,1)%*%f0)%*%thetahat 
  Ypred[Nini] = Ytst
  
  fl = powA(L,1) %*% f0
  var_error_tst = Sigma2 %*% (1 + t(fl)%*%solve(FMatrix)%*%fl)
  std_error_tst = sqrt(var_error_tst)
  
  Ypred_std[Nini] = std_error_tst
}

## Do the plotting and saving it into an image
png(file="LT_tr_est.png",width=800, res=90)  # Save the image in this file, with res =
plot(time, price, type = "l",lwd= 3,main="Linear Trend Estimation", xlab="Time (years)", ylab="Price ",
     ylim = c(minPrice - rangePrice/5, maxPrice + rangePrice/5)) 
# ylim = c(min(price)*(1- 1/10), max(price) + 200)*(1+ 1/10))   
#ylim = c(24000,28000)

## Plot mean +- 2std lines
lines(time_tr[Ninit:Ntr+1], Ypred[Ninit:Ntr+1], col = "blue", lwd=2) 
# Extend the values to create vector and plot
lines(time_tr[Ninit:Ntr+1], Ypred[Ninit:Ntr+1] - 2*Ypred_std[Ninit:Ntr+1] , 
      col = "blue", type = "l", lty=2, lwd=1)
lines(time_tr[Ninit:Ntr+1], Ypred[Ninit:Ntr+1] + 2*Ypred_std[Ninit:Ntr+1] , 
      col = "blue", type = "l", lty=2, lwd=1)
dev.off()

################ Question 5.2 Prediction ######################
# Create "imaginary" time instances 
Xtr_j = cbind(matrix(1,Ntr), matrix(-seq(Ntr-1,0,-1)))
lambda = 0.80

# Transition matrix
L = matrix(1,2,2) 
L[1,2] = 0
# f(0): Imaginary time in instant 0
f0 = matrix(1,2,1) #f(0)
f0[2] = 0

# Weight matrix for assigning weights to the samples.
# Older samples have decreasing weight using exponential decrease.
# It can be seen as WLS where Sigma inverse is a diagonal matrix.
# So it is also a ML for the Multivariate case of all samples, where
# all the samples are independent (diagonal matrix), but have differente
# variances (the more variance, the less weight )

## ESTIMATION (GET Initial parameters)
## Create weight matrix
CovInv = matrix(0,Ntr,Ntr)  # Inverse covariance Matrix of the samples. (ML model)
for(i in 0:(Ntr-1)){
  CovInv[Ntr-i,Ntr-i] = lambda^(i)
}
# Solve the system to obtain the GLS parameters.
# N * Covariance of the dimensions of f(x)_j = [1 j] (Imaginary time)
FMatrix = t(Xtr_j) %*% CovInv %*% Xtr_j 
H = t(Xtr_j) %*% CovInv %*% price_tr      # Samples * Weight_Matrix * Y
thetahat = solve(FMatrix) %*% H  # Same as GLS
Ytr = Xtr_j %*% thetahat # Prediction for training

# The unbias estimator for the variance of the noise in the training samples is:
Ndim = dim(fj)[1]
Sigma2 = (t((price_tr - Ytr)) %*%CovInv %*% (price_tr - Ytr)) %/% (Ntr - Ndim)
std_error_tr = sqrt(Sigma2)

## PREDICTION (GET Initial parameters)
# Now we predict the next samples, using the linear model with respecto to j (time transform)
Ytst = matrix(t(fj)%*%thetahat,1,1)
# For every new sample we get the new time and obtain the output.
fj = f0  # fj for calculating fj for different js
for(i in 1:Ntst){
  fj = L%*%fj     # Obtain new time instance (transformation of time j)
  Ytst_j = t(fj)%*%thetahat 
  Ytst = cbind(Ytst, Ytst_j)
}

## Std of the prediction error
# Tha variance for a prediction "l" samples away from the last training point is 
var_error_tst = matrix(0,1,1)  # Price estimate for tst.
for(l in 1:Ntst){
  fl = powA(L,l) %*% f0
  var_error_tst_l = Sigma2 %*% (1 + t(fl)%*%solve(FMatrix)%*%fl)
  var_error_tst = cbind(var_error_tst, var_error_tst_l)
}
std_error_tst = sqrt(var_error_tst)

## Do the plotting and saving it into an image
png(file="LT_test.png",width=800, res=90)  # Save the image in this file, with res =
plot(time, price, type = "l",lwd= 3,main="House Pricing Evolution", xlab="Time (years)", ylab="Price ",
     ylim = c(22000,28000), xlim = c(2014.5,2016))
# ylim = c(min(price)*(1- 1/10), max(price) + 200)*(1+ 1/10))   
#ylim = c(24000,28000)

## Plot mean +- 2std lines
lines(time_tr, Ytr,col = "blue", lwd=2) # Extend the values to create vector and plot
lines(time_tr, Ytr + 2 * matrix(std_error_tr, Ntr), col = "blue", type = "l", lty=2, lwd=1)
lines(time_tr, Ytr - 2 * matrix(std_error_tr, Ntr), col = "blue", type = "l", lty=2, lwd=1)
## Test prediction and std of tst empirical error 
lines(time_tst, Ytst[-1],col = "red", lwd=2) # Extend the values to create vector and plot
lines(time_tst, Ytst[-1] + 2 * matrix(std_error_tst[-1], Ntst), col = "red", type = "l", lty=2, lwd=1)
lines(time_tst, Ytst[-1] - 2 * matrix(std_error_tst[-1], Ntst), col = "red", type = "l", lty=2, lwd=1)
dev.off()



##################################################################
################ Question 6. Unefficient way #####################
##################################################################

## Brute force, no updating.
# Calculate the estimation error in the training set
# Transition matrix

L = matrix(1,2,2) 
L[1,2] = 0
# f(0): Imaginary time in instant 0
f0 = matrix(1,2,1) #f(0)
f0[2] = 0
fj = L %*% f0

# Search of lambda
Npoints = 20;
lambda_values = (1:Npoints)/(Npoints+1)
S_values = matrix(0,Npoints,1)   # Where to store the errors :)

Sindx = 0
for (lambda in (lambda_values)){
  Sindx = Sindx +1
  ## ESTIMATION (GET Initial parameters)
  ## Create weight matrix
  Ninit = 20
  Nest = Ntr - Ninit   # Number of samples to use for updating
  # Main variables where to store the predictions and theoretical std at every point.
  Ypred = matrix(1,Nest,1)
  Ypred_std = matrix(1,Nest,1)
  
  S_aux = 0 
  for (Nini in (Ninit:(Ntr-1))){
    X_j = cbind(matrix(1,Nini), matrix(-seq(Nini-1,0,-1)))
    CovInv = matrix(0,Nini,Nini)  # Inverse covariance Matrix of the samples. (ML model)
    for(i in 0:(Nini-1)){
      CovInv[Nini-i,Nini-i] = lambda^(i)
    }
    # Solve the system to obtain the GLS parameters.
    # N * Covariance of the dimensions of f(x)_j = [1 j] (Imaginary time)
    FMatrix = t(X_j[1:Nini,]) %*% CovInv %*% X_j[1:Nini,]
    H = t(X_j[1:Nini,]) %*% CovInv %*% price_tr[1:Nini,]      # Samples * Weight_Matrix * Y
    thetahat = solve(FMatrix) %*% H  # Same as GLS
    Ytr = X_j %*% thetahat
    
    # Get next test and its std
    Ytst = t(powA(L,1)%*%f0)%*%thetahat 
    
    S_aux = S_aux + (Ytst - price_tr[Nini+1])^(2)
    print (S_aux)
  }
  S_values[Sindx] = S_aux
}

best_lambda = lambda_values[which.min(S_values)]

## Do the plotting and saving it into an image
png(file="S_value.png",width=800, res=90)  # Save the image in this file, with res =
plot(lambda_values, S_values, type = "l",lwd= 3,main="Lambda value", xlab="Lambda", ylab="Square Error")

dev.off()

##################################################################
################ Question 5 Efficcient Update Misserably Failed ##
##################################################################

lambda = 0.99
# Transition matrix
L = matrix(1,2,2) 
L[1,2] = 0
# f(0): Imaginary time in instant 0
f0 = matrix(1,2,1) #f(0)
f0[2] = 0
fj = L %*% f0
# Weight matrix for assigning weights to the samples.
# Older samples have decreasing weight using exponential decrease.
# It can be seen as WLS where Sigma inverse is a diagonal matrix.
# So it is also a ML for the Multivariate case of all samples, where
# all the samples are independent (diagonal matrix), but have differente
# variances (the more variance, the less weight )

## ESTIMATION (GET Initial parameters)
## Create weight matrix
Nini = 20; # Initial number of samples to perform the initial GLS
# Create "imaginary" time instances for ALL time instances
X_j = cbind(matrix(1,Nsam), matrix(-seq(Nsam-1,0,-1)))
CovInv = matrix(0,Nini,Nini)  # Inverse covariance Matrix of the samples. (ML model)
for(i in 0:(Nini-1)){
  CovInv[Nini-i,Nini-i] = lambda^(i)
}

# Solve the system to obtain the GLS parameters.
# N * Covariance of the dimensions of f(x)_j = [1 j] (Imaginary time)
FMatrix = t(X_j[1:Nini,]) %*% CovInv %*% X_j[1:Nini,]
H = t(X_j[1:Nini,]) %*% CovInv %*% price_tr[1:Nini,]      # Samples * Weight_Matrix * Y
thetahat = solve(FMatrix) %*% H  # Same as GLS
Ytr = X_j[1:Nini,] %*% thetahat # Prediction of the previous samples

# The unbias estimator for the variance of the last samples is:
Ndim = dim(fj)[1]
Sigma2 = (t((price[1:Nini,] - Ytr)) %*%CovInv %*% (price[1:Nini,] - Ytr)) %/% (Nini - Ndim)
std_error_tr = sqrt(Sigma2)

## Estimation of the next sample.
# Now we predict only the next sample, using the linear model with respecto to j (time transform)
Ytst = matrix(t(fj)%*%thetahat,1,1)
# For every new sample we get the new time and obtain the output.
fj = L%*%f0   # fj for calculating fj for different js
Ytst = t(fj)%*%thetahat 

## Std of the prediction error
# Tha variance for a predictnion "l" samples away from the last training point is 
fl = powA(L,l) %*% f0
var_error_tst_l = Sigma2 %*% (1 + t(fl)%*%FMatrix%*%fl)
std_error_tst = sqrt(var_error_tst)

###################################################
##### UPDATING OF PARAMETERS WITH NEW SAMPLES #####
###################################################

Nest = Ntr - Nini   # Number of samples to use for updating
# Main variables where to store the predictions and theoretical std at every point.
Ypred = matrix(1,Nest,1)
initialthetahat = thetahat
initialFMatrix = FMatrix
initialH = H

### REAL UPDATION
thetahat = initialthetahat
FMatrix = initialFMatrix 
H = initialH 

for(i in 1:(Nest)){   # For every new sample
  #  print(thetahat)
  fi = powA(L,1) %*% f0  # We update the 0 point when we update this
  #### Get the estimate for the sample (from the past updating)
  Ypred[i] = t(fi)%*%thetahat 
  #  print (Ypred[i])
  #### UPDATE THE parameters with the new sample !
  # We update F, then H and then the parameters
  if (i < 5){
    print ("PUTA--------------------------")
    print(FMatrix)
    print (H)
    print (thetahat)
  }
  # Get the H and F
  fl = powA(L,-(Nini + i -1)) %*% f0  
  FMatrix = FMatrix + lambda^(Nini + i -1) * fl%*%t(fl)
  H = lambda * solve(L) %*% H + f0 * price[Nini + i, ]
  thetahat = solve(FMatrix) %*% H  # Same as GLS
  # thetahat = t(L)%*%thetahat + solve(FMatrix) %*% f0 *(price[Nini + i, ] - Ypred[i])  # Same as GLS
}

## Do the plotting and saving it into an image
png(file="Local_Trend.png",width=800, res=90)  # Save the image in this file, with res =
plot(time, price, type = "l",lwd= 3,main="House Pricing Evolution", xlab="Time (years)", ylab="Price ")
# ylim = c(min(price)*(1- 1/10), max(price) + 200)*(1+ 1/10))   
#ylim = c(24000,28000)

## Plot mean +- 2std lines
lines(time[Nini:(Ntr-1)], Ypred[,1] ,col = "blue", lwd=2) # Extend the values to create vector and plot
dev.off()

