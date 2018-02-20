
plot_timeSeries <- function(X,Y, name, xlabel, ylabel, lwd, gp){
  # gp: If set to 1, then it generates a physical image
  if (gp == 1){
    png(file = name,width=800, height = 600, res=130) 
  }
  # Size and aspect ratio could be width=400,height=350,res=45
  plot(X, Y, 
       type = "l",                     # Draw as a line
       lwd= lwd,                         # Line width
       main= name, # Title of the graph
       xlab=xlabel,            # x label
       ylab=ylabel
       
       )                  # y label
  
  if (gp == 1){
    dev.off()
  }
}

plot_acfpacf <- function(X, name, Nlag, gp){
  # gp: If set to 1, then it generates a physical image
  if (gp == 1){
    png(file = name,width=800, height = 600, res=130) 
  }
  par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(2,0.7,0))
  acf(X, lag.max = Nlag )
  pacf(X, lag.max = Nlag )
  
  #title(name); 
        
  if (gp == 1){
    dev.off()
  }
}

plot_SigACFPACF <- function(residuals, Nlag){
  par(mfrow=c(3,1), mar=c(3,3,1,1), mgp=c(2,0.7,0))
  plot(residuals, type="l")
  acf(residuals, lag.max = Nlag)
  pacf(residuals, lag.max = Nlag)
}

################ Question 3-1 ######################
## Loading and preparing the data:
myData = read.csv("./A3_jagt_NOx.csv", sep = ",") # read the csv
date = as.matrix(myData[,1])   # Get the date
hour = as.matrix(myData[,2])   # Get the hour
NO = as.matrix(myData[,3])     # Get the NO concentration 

Nsam = dim(NO)[1]             # Number of samples

# Convert date and hour to timestamp
# casa = as.matrix(strsplit(hour, "-"))[1,]    # Useless
hour <-substr(as.character(hour),2,3)
myTimeStamp <- as.POSIXlt(paste(as.character(date),hour),format="%d-%m %H")

## We reverse the data !!
myTimeStamp = rev(myTimeStamp)
NO = rev(NO)

NO = NO + 200   # Supermean


# Vars for plotting
minNO = min(NO)
maxNO = max(NO)
rangeNO = maxNO - minNO

### Split data into train and test
Ntst = 48
tr_indx = 1:(Nsam - Ntst)  # Indexes for estimating
Ntr = length(tr_indx)
tst_indx = (Ntr + 1):Nsam   # Indexes for estimating
Ntst = length(tst_indx)

# Obtain data
myTimeStamp_tr = myTimeStamp[tr_indx]
myTimeStamp_tst = myTimeStamp[tst_indx]

NO_tr = NO[tr_indx]
NO_tst = NO[tst_indx]

## Do the plotting and saving it into an image
plot_timeSeries(myTimeStamp, NO, "NO concentration", "Time", "Concentration", 2, 0)

Nmaxsam = Ntr # 
plot_timeSeries(myTimeStamp_tr[1:Nmaxsam], NO_tr[1:Nmaxsam], "NO concentration", 
                "Time", "Concentration", 2, 0)
## Do the plotting and saving it into an image
Nmaxsam = 24*7*2
plot_timeSeries(myTimeStamp_tr[(Ntr -Nmaxsam):Ntr], NO_tr[(Ntr -Nmaxsam):Ntr], 
                "NO concentration", "Time", "Concentration", 2, 0)

############## Transformations of the data !##############
diff_NO_tr = diff(NO_tr)
diff_NO_tr2 = diff(NO_tr, lag = 24)
log_NO_tr = log(NO_tr)

# The next two are the same
difflog_NO_tr = diff(log_NO_tr)  # -1 due to the return
normdiff_NO_tr = log( 1 + diff_NO_tr / NO_tr[-Ntr]) # -1 due to the return

## Difference twice!!
diffdiff_NO_tr = diff(diff_NO_tr, lag = 24)

## Get rid of part of the seasonal component !!
n = 24   # Differentiation
diff_NO_tr2 = NO_tr[-1] - 0.9*NO_tr[-Ntr]
Ntr_aux = Ntr - 1
diff_NO_tr2 = diff_NO_tr[(n+1):(Ntr_aux)] - 0.8*diff_NO_tr[1:(Ntr_aux-n)]

##############################################
######  PLOTTING
par(mfrow=c(1,1))

plot_timeSeries(myTimeStamp_tr, log_NO_tr, "logNOtr", "Time", 
                "Increase in NO concentration",2,0)

plot_timeSeries(myTimeStamp_tr[-Ntr], diff_NO_tr, "diffNOtr", "Time", 
                "Increase in NO concentration",2,0)

plot_timeSeries(myTimeStamp_tr[-Ntr], difflog_NO_tr, "difflogNOtr", "Time", 
                "Increase in NO concentration",2,0)
plot_timeSeries(myTimeStamp_tr[1:(Ntr-25)], diff_NO_tr2, "diff_NO_tr2", "Time", 
                "Increase in NO concentration",2,0)

plot_timeSeries(myTimeStamp_tr[1:(Ntr-25)], diffdiff_NO_tr, "diffdiff_NO_tr", 
                "Time", "Increase in NO concentration",2,0)
plot_timeSeries(myTimeStamp_tr[1:(Ntr-25)], diff_NO_tr2, "diff_NO_tr2", 
                "Time", "Increase in NO concentration",2,0)

#######################################
###### CHECK STATIONARITY
library(tseries)
adf.test(NO_tr, alternative = "stationary",k = 24)
adf.test(diff_NO_tr, alternative = "stationary",k = 24)
adf.test(log_NO_tr, alternative = "stationary",k = 24)
# Since the p-value is less than 0.01, and being the alrernatyive hyp 
# that the signal is stationary. The p-value, the probability of observing
# data more extreme that ours given the null-hypotesis 
# (opposite to alternative hypotesis)
# Since the probability of observing data that is more extreme is very low,
# Then the null hypotesis is rejected, because the data obtained is very unlikely
# under the null hypotesis.

## RANGE-MEAN


par(mfrow=c(2,2), mgp=c(2,0.7,0))
LenS = 12
Ranges = c()
Means = c()

lens = c(4,12,18,24)
for( LenS in lens){
  for (i in 1:(Ntr - LenS)){
    window = NO_tr[i:(i+LenS)]
    range = max(window) - min(window)
    mean = mean(window)
    Ranges = c(Ranges, range)
    Means = c(Means, mean)
  }
  plot(Means,
       Ranges, 
       lwd= lwd,                         # Line width
       main= paste("Mean-Range l = ", toString(LenS)), # Title of the graph
       xlab="Mean",            # x label
       ylab="Range")
}

################ Question 2 ######################
# We just use the ACF and PACF functions to calculate the coefficients.

## Original data !!
Nlag = 200
plot_acfpacf(NO_tr,"NOtr_ACFPACF.png",Nlag, 0)
plot_acfpacf(diff(log_NO_tr),"logNOtr_ACFPACF.png",Nlag, 0)
plot_acfpacf(diff_NO_tr,"diff_ACFPACF.png",Nlag, 0)
plot_acfpacf(difflog_NO_tr,"diff_ACFPACF.png",Nlag, 0)
plot_acfpacf(diff_NO_tr2, "diffdiff_ACFPACF.png",Nlag, 0)
plot_acfpacf(diffdiff_NO_tr, "diffdiff_ACFPACF.png",Nlag, 0)

dev.off()
################ Question 3 ######################
# Let us try to fit a model
# Model with ARIMA(2,0,2)(2,0,2) with seasonal component of order 1
# (0,1,1)(0,1,1)

fit1 <- arima(log_NO_tr, order=c(1,0,1), include.mean = FALSE, method="ML", 
              seasonal = list(order = c(1,1,1),period=24))

# require(forecast)
# ARIMAfit <- auto.arima(NO_tr, approximation=FALSE,trace=FALSE)
# Summary(ARIMAfit)

plot(fit1$coef)    # Plot the coefficientes obtained.
residuals = fit1$residuals
coef = fit1$coef
varMatrix_coef = fit1$var.coef

N_coef = length(coef)

# Calculate variance of the coefficients
vars_coef = c(1:N_coef)
for (i in 1:N_coef){
  vars_coef[i] = varMatrix_coef[i,i]
}
sigma_coef = sqrt(vars_coef)
T_value_coef = coef/sigma_coef
degrees_freedom = Ntr - 1 - N_coef
p_values = 2*pt(-abs(T_value_coef),df=degrees_freedom)

## More param
sigma_res = fit1$sigma2
loglik = fit1$loglik

## Error
SS <- sum(residuals^2)/Ntr  # Sum of square errors of the residual model 1

AIC =  Ntr * log(SS) + 2 * N_coef
BIC = Ntr * log(SS) +  N_coef * log(Ntr)

#### PREDICT
prediction = predict(fit1,48)
SSval = sum((prediction$pred - NO_tst)^2)/Ntst

par(mfrow=c(1,1))

#fit1 <- arima(diff_NO_tr, order=c(25,0,1), include.mean = FALSE, method="ML")
#plot(fit1$coef)    # Plot the coefficientes obtained.
#residuals = fit1$residuals
# This is not very good because it also obtaines valus for the inbetween ACF.

## Properties
## we want that collection of plots a lot so let's make a function:

view_residuals = 0
if (view_residuals == 1){
  Nlag2 = 30
  #### Graph1: Signal + ACF + PACF
  plot_SigACFPACF(residuals, Nlag)
  tsdiag(fit1, gof.lag = Nlag2)  # You plot the residuals and their ACF and JunlgeBox
  #  Ljung-Box Q null hypothesis is that there is no autocorrelation in the errors
  #where  \chi _{1-\alpha ,h}^{2}} \chi_{1-\alpha,h}^2 is the α-quantile of the chi-squared 
  # distribution with h degrees of freedom
  #### Graph2: Histogram and sample quantiles
  par(mfrow=c(1,2))
  
  # Histogram and draw gaussian
  hist(residuals,probability=T,col='blue')
  curve(dnorm(x,sd = sqrt(fit1$sigma2)), col=2, lwd=2, add = TRUE)
  # Do the samples quantities
  qqnorm(fit1$residuals)
  qqline(fit1$residuals)
  
  ## Graph 3: Simulate Residuals
  par(mfrow=c(2,1))
  ts.plot(residuals, ylim = c(-80,80))
  ts.plot(ts(rnorm(length(residuals))*sqrt(fit1$sigma2)),
          ylab='Simulated residuals',ylim = c(-80,80))
}

test_residuals = 0
if (test_residuals == 1){
  ### TESTS FOR THE RESIDUAL
  
  ### Binomial test
  # Confidence interval of the number of change of signs.
  n.residuals <- length(residuals)
  (n.residuals-1)/2
  sqrt((n.residuals-1)/4) 
  (n.residuals-1)/2 + 1.96 * sqrt((n.residuals-1)/4) * c(-1,1)
  
  ### Binomial test
  # P-value using binomial distribution that the number of sign changes is 1/2 probable
  (N.sign.changes <- sum( residuals[-1] * residuals[-n.residuals]<0 ))
  bt = binom.test(N.sign.changes, n.residuals-1)
  p_value = bt$p.value
  int_conf = bt$conf.int
  
  ### Check if the sum of squared acf values follow a chi-square distribution.
  acfvals <- acf(residuals, type="correlation", plot=FALSE)$acf[2:24] # Get acf residuals
  test.stat <- sum(acfvals^2) * (length(residuals)-1)
  ## Numerical issue ... next line is better
  1 - pchisq(test.stat, length(acfvals)-length(fit1$coef))
  prob = pchisq(test.stat, length(acfvals)-length(fit1$coef), lower.tail = FALSE) 

   ## Log likelihood of the data
  loglik = fit1$loglik
  # Informaiton creiterias
  res = residuals
  Nres = length(res)
  nparam = length(fit1$coef)  
  
  VarEstimate = sum(res^2)/(Nres - nparam) 
  VarEstimate2 = fit1$sigma2
  
  # Cumulative periodogram. The cumulative value of the frequency components of the signal.
  # If the signal is noise, then it should be a straing line since all the
  # frequency componentes should have the same value.
  

  par(mfrow=c(1,2))
  cpgram(NO_tr)
  cpgram(residuals)
  
}

compare_models = 0
if (compare_models == 1){
  
  # Likelihood ratio test
  prob_chisq = pchisq(-2* ( fit1$loglik - fit2$loglik ), df=1, lower.tail = FALSE)
  
  # F-test for lower order models
  s1 <- sum(residuals^2)  # Sum of square errors of the residual model 1
  s2 <- sum(residuals^2) # Sum of square errors of the residual model 2
  n1 <- 3    # Number of param for model 1
  n2 <- 4    # Number of param for model 2
  
  # Calculate the F-value of the models
  F_value = pf( (s1-s2)/(n2-n1) / (s2/(length(residuals)-n2)), 
                df1 = n2 - n1, df2 = (length(residuals)-n2), lower.tail = FALSE)
}

#### PREDICT

predic_shit = 1
if (predic_shit == 1){
  prediction = predict(fit1,48)
  VarEstimate = sum((prediction$pred - NO_tst)^2)
  
  Nbefore = 80
  
  lwd = 2
  
  plot(myTimeStamp[(Ntr -Nbefore):(Ntr+Ntst)],
       NO[(Ntr -Nbefore):(Ntr+Ntst)], 
       type = "l",                     # Draw as a line
       lwd= lwd,                         # Line width
       main= "NO concentration", # Title of the graph
       xlab="Time",            # x label
       ylab="Concentration",
       ylim = c(-200,300))
  
  lines(myTimeStamp_tst,  prediction$pred, col="blue")
  lines(myTimeStamp_tst,  prediction$pred + 1 * prediction$se, col="blue", lty=2)
  lines(myTimeStamp_tst,  prediction$pred - 1 * prediction$se, col="blue", lty=2)
  
  #legend(myTimeStamp[Ntr],200, legend = c("Actual values","Predicted values"))
  
  prse = c(prediction$se[1], prediction$se[24], prediction$se[48])
  pr = c(prediction$pred[1], prediction$pred[24], prediction$pred[48])
  
}

########################################## LOGGING #####################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
No = 1000
x = (1:No)/(No/10)
lnx = log(x)
dlnx = diff(lnx)*(No/10)
#dlnx = 1/x[-No]

par(mfrow=c(1,2))

plot(x,
     lnx, 
     type = "l",                     # Draw as a line
     lwd= lwd,                         # Line width
     main= "ln(x)", # Title of the graph
     xlab="x",            # x label
     ylab="Logarithmic transformation",
     xlim = c(0,4))

plot(x[-No],
     dlnx, 
     type = "l",                     # Draw as a line
     lwd= lwd,                         # Line width
     main= "derivative of ln(x)", # Title of the graph
     xlab="x",            # x label
     ylab="dln(x)",
     xlim = c(0,2),
     ylim = c(0,5))



##### IF WE ARE USING LOG OF THE INITIAL ISGNAL ############3

NO_tr_estimation = exp(-(fit1$residuals - log_NO_tr))
residuals = (NO_tr - NO_tr_estimation)
SS <- sum(residuals^2)/Ntr  # Sum of square errors of the residual model 1

plot_SigACFPACF(residuals, Nlag)
tsdiag(fit1, gof.lag = Nlag2)  # You plot the residuals and their ACF and JunlgeBox

par(mfrow=c(1,2))

# Histogram and draw gaussian
hist(residuals,probability=T,col='blue')
curve(dnorm(x,sd = sqrt(fit1$sigma2)), col=2, lwd=2, add = TRUE)
# Do the samples quantities
qqnorm(fit1$residuals)
qqline(fit1$residuals)

par(mfrow=c(2,1))
plot_timeSeries(myTimeStamp_tr, log_NO_tr, "logNOtr", "Time", 
                "Increase in NO concentration",2,0)

plot_timeSeries(myTimeStamp_tr, NO_tr_estimation, "logNOtr", "Time", 
                "Increase in NO concentration",2,0)



############### LOGSHIT ONLY !!!!! ######################
prediction = predict(fit1,48)
NO_tst_pred = exp(prediction$pred)
residuals_tst = (NO_tst - NO_tst_pred)

VarEstimate = sum((NO_tst - NO_tst_pred)^2)/(Ntst)


Nbefore = 80

lwd = 2
par(mfrow=c(1,1))
plot(myTimeStamp[(Ntr -Nbefore):(Ntr+Ntst)],
     NO[(Ntr -Nbefore):(Ntr+Ntst)], 
     type = "l",                     # Draw as a line
     lwd= lwd,                         # Line width
     main= "NO concentration", # Title of the graph
     xlab="Time",            # x label
     ylab="Concentration",
     ylim = c(0,300))

lines(myTimeStamp_tst,  NO_tst_pred, col="blue")
lines(myTimeStamp_tst,  exp((prediction$pred + 1 * prediction$se)), col="blue", lty=2)
lines(myTimeStamp_tst,  exp((prediction$pred - 1 * prediction$se)), col="blue", lty=2)

#legend(myTimeStamp[Ntr],200, legend = c("Actual values","Predicted values"))

prse = c(prediction$se[1], prediction$se[24], prediction$se[48])
pr = c(prediction$pred[1], prediction$pred[24], prediction$pred[48])

       