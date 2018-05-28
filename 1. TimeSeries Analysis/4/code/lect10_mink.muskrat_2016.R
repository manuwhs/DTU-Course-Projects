## Example: Mink vs Muskrat
setwd("/home/laec/teaching/02417/2016E/lectures/week10")

## Data source:
##   Jones, J.W. (1914) ”Fur-farming in Canada”, Commission of
##   Conservation Canada, pp.209–214
##   Fur sales, 1850-1911
##   URL: http://robjhyndman.com/tsdldata/ecology1/

years <- 1850:1911
mink <- read.table("mink.dat",skip = 3)[[1]] # To get the only column of data
muskrat <- read.table("muskrat.dat",skip = 3)[[1]]
mink <- ts(mink, start=1850, frequency=1)
muskrat <- ts(muskrat, start=1850, frequency=1)
data <- cbind(mink,muskrat)
# write.table(cbind(years,mink,muskrat), file="MinkMuskrat.txt",sep="\t",row.names = FALSE)

plot(data)
summary(data)

# looks like increasing variance with increasing values so log transforming
ldata <- log(data)
summary(ldata)  # Ups ... the names should reflect the transformation.
colnames(ldata) <- c("lmink","lmuskrat") ## Better names
summary(ldata)
plot(ldata)
par(mgp=c(2,0.7,0), mar=c(3,3,3,1))
acf(ldata)
pacf(ldata)


var.auto <- ar(ldata, order.max = 10, method="mle")
var.auto <- ar(ldata, order.max = 10)
var.auto # VAR(5)
summary(var.auto) # Only provides an overview of what is there ...

var.auto$aic # Which model order is best and by how much!

## How about the residuals?
acf(var.auto$resid, na.action = na.omit)
pacf(var.auto$resid, na.action = na.omit)

ks.test(var.auto$resid[,1]/sd(var.auto$resid[,1],na.rm = TRUE),"pnorm")
ks.test(var.auto$resid[,2]/sd(var.auto$resid[,2],na.rm = TRUE),"pnorm")
## Cannot reject independence nor normality

###############
## This was using "ar" but there is an other package just for this ...
###############
install.packages("vars")
library(vars)

( var1 <- VAR(ldata, p=5) )
## Not quite the same estimates - due to different technique
summary(var1)

( vars.l <- VARselect(ldata, type="const", lag.max = 6) )
## SC is the same as BIC...
## http://sccn.ucsd.edu/wiki/Chapter_3.5._Model_order_selection

( var1l <- VAR(ldata, p=1) )
summary(var1l)



## Test for auto correlation
serial.test(x = var1l)

## Test for conditional heteroskedasticity
arch.test(var1l)

## Three tests for normality
normality.test(var1l)
vars.test(var1l)

par(mgp=c(2,0.7,0), mar=c(3,3,3,1))

plot(var1l)
summary(var1l)
## Note that the roots (At top of summary) have the same modulus.
## Here comes the actual roots:
roots(var1l, modulus=FALSE)
## [1] 0.7347993+0.3131011i 0.7347993-0.3131011i
## They are complex conjugated. So a two dimensional VAR(1) model can produce oscilations 
## which is not the case for an univariate AR(1) - An AR(2) is required and it can be
## written as a VAR(1) :-)

VARselect(ldata, type="const", lag.max = 6) # - 5.02
VARselect(ldata, type="trend", lag.max = 6) # -4.87
VARselect(ldata, type="both", lag.max = 6) # -5.12
VARselect(ldata, type="none", lag.max = 6) # -4.81

var1b <- VAR(ldata, p=1, type="both") 
summary(var1b)

## One has a significant constant and the other a signicant trend ... 
## but you cannot easily pick which parameters to keep.
