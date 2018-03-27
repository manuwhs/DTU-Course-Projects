# Marima.doc code

install.packages("marima")
library(marima)

k <- 3
p <- 2
dependent <- paste("y", c(1:k))
regressor <- paste("r", c(1:k))
order <- paste("o", c(0:p))
phi<-  array(data=NA, dim=c(k,k,(1+p)), dimnames=list(dependent,regressor,order))

## Example: Australian firearms ####
data(austr); 
all.data <- austr

old.data <- t(austr)[,1:83]
ar<-c(1)
ma<-c(1)
# Define the proper model:
Model1 <- define.model(kvar=7, ar=ar, ma=ma, rem.var=c(1,6,7), indep=c(2:5))
# Now call marima:
Marima1 <- marima(old.data,means=1,
                  ar.pattern=Model1$ar.pattern, ma.pattern=Model1$ma.pattern,
                  Check=FALSE, Plot='log.det', penalty=0.0)
short.form(Marima1$ar.estimates, leading=FALSE) # print estimates
short.form(Marima1$ma.estimates, leading=FALSE)

## Now a true multivariat model
Model2 <- define.model(kvar=7, ar=ar, ma=ma, rem.var=c(1,6,7), indep=NULL)
Marima2 <- marima(old.data, means=1, ar.pattern=Model2$ar.pattern,
                  ma.pattern=Model2$ma.pattern, Check=FALSE, Plot='log.det', penalty=0)
short.form(Marima2$ar.estimates, leading=FALSE) # print estimates
short.form(Marima2$ma.estimates, leading=FALSE)

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


# call the forecasting function using Marima and the prepared data:
Forecasts <-  arma.forecast(t(all.data[1:100,]), nstart=90, nstep=10, marima=Marima4)

### From here on the plot is constructed ###
Year<-t(all.data[91:100,1]);
Predict<-Forecasts$forecasts[2,91:100]
stdv<-sqrt(Forecasts$pred.var[2,2,])
upper.lim=Predict+stdv*1.645
lower.lim=Predict-stdv*1.645
Out<-rbind(Year,Predict,upper.lim,lower.lim)
print(Out)
# plot results:
plot(all.data[1:100,1], Forecasts$forecasts[2,],type="l", xlab="Year",
     ylab="Rate of armed suicides", main="Prediction of suicides by firearms",
     ylim=c(0.0,4.1))
lines(all.data[1:90,1], all.data[1:90,2], type="p")
grid(lty=3, lwd=1, col="grey")
Years<-2005:2014
lines(Years, Predict, type="l")
lines(Years, upper.lim, type="l")
lines(Years, lower.lim, type="l")
lines(c(2004,2004), c(0,2))

