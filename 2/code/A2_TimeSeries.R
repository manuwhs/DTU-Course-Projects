
# 998 714 1009 1296 1083 751 950 1267 1134 838
mu = 1000
phi = c(0.8,-0.7,0.9,-0.72,0.63)
Yt1 = phi[1]*(838 - mu) +  phi[2]*(1134 - mu) + phi[3]*(950 - mu) + phi[4]*(751 - mu) +phi[5]*(1083 - mu)  + mu
Yt2 = phi[1]*(Yt1 - mu) +  phi[2]*(838 - mu) + phi[3]*(1267 - mu) + phi[4]*(950 - mu) +phi[5]*(751 - mu)  + mu

phi1 <- seq(from = -2.5, to = 2.5, length = 51) 
plot(phi1,1+phi1,lty="dashed",type="l",xlab="",ylab="",cex.axis=.8,ylim=c(-1.5,1.5))
abline(a = -1, b = 0, lty="dashed")
abline(a = 1, b = -1, lty="dashed")
title(ylab=expression(phi[2]),xlab=expression(phi[1]),cex.lab=.8)
polygon(x = phi1[6:46], y = 1-abs(phi1[6:46]), col="gray")
lines(phi1,-phi1^2/4)
text(0,-.5,expression(phi[2]<phi[1]^2/4),cex=.7)
text(1.2,.5,expression(phi[2]>1-phi[1]),cex=.7)
text(-1.75,.5,expression(phi[2]>1+phi[1]),cex=.7)


################ Question 1 ######################
Nts = 10    # Number of timeSeries
Nsa = 1000  # Number of samples for the timeSeries
par(mfrow=c(1,2)) 
Processes = matrix(0,Nts,Nsa)
time_i = 1:Nsa

mu = 1.0/4

for (i in 1:Nts){
    noise <- rnorm(Nsa)
    process = cumsum(noise)
    Processes[i,] = process + mu
}

## PLOTTING
png(file="Random Process Samples",width=800, height = 600, res=130) 
# Size and aspect ratio could be width=400,height=350,res=45
plot(time_i, Processes[i,], 
     type = "l",                     # Draw as a line
     lwd= 1,                         # Line width
     main="Random Process Samples.png", # Title of the graph
     xlab="t",            # x label
     ylab="Yt",
     ylim = c(min(Processes),max(Processes)))                  # y label

for (i in 2:Nts){
  lines(time_i, 
        Processes[i,], col=sample(rainbow(100)), type = "l", lwd=1)
}
dev.off()

### Autocorrelations ! 
Nlag = 80    # Maximum number of lag for autocorrelation
Autocorrelations = matrix(0,Nts,Nlag+1)
lag_i = 0:Nlag
for (i in 1:Nts){
  afc_obj = acf(Processes[i,],  lag.max = Nlag, plot = FALSE)
  Autocorrelations[i,] = afc_obj$acf
}
#sd = colMeans(Autocorrelations, na.rm = FALSE, dims = 1)
png(file="Autocorrelation.png",width=800, height = 600, res=130) 
plot(lag_i, Autocorrelations[1,], 
     type = "h",                     # Draw as a line
     lwd= 1,                         # Line width
     main="Autocorrelation", # Title of the graph
     xlab="lag",            # x label
     ylab="Autocorrelation of 10 realization",
     ylim = c(0,1.2))              # y label

for (i in 2:Nts){
  lines(lag_i, 
        Autocorrelations[i,], col=sample(rainbow(100)), type = "h", lwd=1)
}
dev.off()

################ Question 2.4 ######################

########## Model 1
# Simulating an AR(1) 
phi = c(0.9)
s = 12
Nlag = 40
Nsa = 1000

sim <- arima.sim(model = list(ar = phi, order=c(1,0,0)), n = Nsa)

png(file="SM1_Realizations.png",width=800, height = 600, res=130) 
par(mfrow=c(1,1)) 
plot(sim, lwd= 2)
dev.off()

png(file="SM1_ACFPACF.png",width=1000, height = 500, res=130) 
par(mfrow=c(1,2)) 
acf(sim, lag.max = Nlag )
pacf(sim)
dev.off()

########## Model 2
# Simulating an AR(1) 
Phi = c(0.7)
ar_coeff = c(matrix(0,1,s-1),Phi)
Nsa = 1000
sim <- arima.sim(model = list(ar = ar_coeff, order=c(12,0,0)), n = Nsa)

png(file="SM2_Realizations.png",width=800, height = 600, res=130) 
par(mfrow=c(1,1)) 
plot(sim)
dev.off()

png(file="SM2_ACFPACF.png",width=1000, height = 500, res=130) 
par(mfrow=c(1,2)) 
acf(sim, lag.max = Nlag , main = NULL, ylab = NULL)
pacf(sim, lag.max = Nlag, main = NULL , ylab = NULL)
dev.off()

########## Model 3
# Simulating an AR(1) 
Theta = c(0.4)
phi = c(0.9)
s = 12
Nsa = 100000
ar_coeff = c(phi)
ma_coeff = c(matrix(0,1,s-1),Theta)

sim <- arima.sim(model = list(ar = ar_coeff, order=c(1,0,12), ma = ma_coeff), n = Nsa)

png(file="SM3_Realizations.png",width=800, height = 600, res=130) 
par(mfrow=c(1,1)) 
plot(sim)
dev.off()

png(file="SM3_ACFPACF.png",width=1000, height = 500, res=130) 
par(mfrow=c(1,2)) 
acf(sim, lag.max = Nlag )
pacf(sim, lag.max = Nlag )
dev.off()

########## Model 4
# Simulating an AR(1) 
Phi = c(0.7)
s = 12
ar_coeff = c(matrix(0,1,s-1),Phi,-Phi*phi)
ar_coeff[1] = phi

sim <- arima.sim(model = list(ar = ar_coeff, order=c(13,0,0)), n = Nsa)

png(file="SM4_Realizations.png",width=800, height = 600, res=130) 
par(mfrow=c(1,1)) 
plot(sim)
dev.off()

png(file="SM4_ACFPACF.png",width=1000, height = 500, res=130) 
par(mfrow=c(1,2)) 
acf(sim, lag.max = Nlag )
pacf(sim, lag.max = Nlag )
dev.off()

########## Model 5
# Simulating an AR(1) 
Theta = c(0.3)
theta = c(0.4)
s = 12
Nsa = 10000
ma_coeff = c(matrix(0,1,s-1),Theta,-Theta*theta)
ma_coeff[1] = theta

sim <- arima.sim(model = list(order=c(0,0,13), ma = ma_coeff), n = Nsa)

png(file="SM5_Realizations.png",width=800, height = 600, res=130) 
par(mfrow=c(1,1)) 
plot(sim)
dev.off()

png(file="SM5_ACFPACF.png",width=1000, height = 500, res=130) 
par(mfrow=c(1,2)) 
acf(sim, lag.max = Nlag )
pacf(sim, lag.max = Nlag )
dev.off()

########## Model 6
# Simulating an AR(1) 
theta = c(0.4)
Phi = c(0.7)

ar_coeff = c(matrix(0,1,s-1),Phi)
ma_coeff = theta
s = 12
Nsa = 1000
sim <- arima.sim(model = list(order=c(12,0,1), ma = ma_coeff, ar = ar_coeff), n = Nsa)

png(file="SM6_Realizations.png",width=800, height = 600, res=130) 
par(mfrow=c(1,1)) 
plot(sim)
dev.off()

png(file="SM6_ACFPACF.png",width=1000, height = 500, res=130) 
par(mfrow=c(1,2)) 
acf(sim, lag.max = Nlag )
pacf(sim, lag.max = Nlag )
dev.off()


