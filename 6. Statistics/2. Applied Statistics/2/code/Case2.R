#Files
library(xtable)
library(MASS)  #atach package
library(smooth)
cam <- read.delim("case2regionsOnePerBatch.txt")

# my.dataframe[ , "new.col"] <- a.vector
#cam<-cam[,-(8:25)]
cam$rp<-(Camlob$pos/Camlob$total);
cam$rp1<-(Camlob$R1pos/Camlob$R1total);
cam$rp2<-(Camlob$R2pos/Camlob$R2total)
cam$rp3<-(Camlob$R3pos/Camlob$R3total);
cam$rp4<-(Camlob$R4pos/Camlob$R4total);
cam$rp5<-(Camlob$R5pos/Camlob$R5total)
cam$rp6<-(Camlob$R6pos/Camlob$R6total);
cam$rp7<-(Camlob$R7pos/Camlob$R7total);
cam$rp8<-(Camlob$R8pos/Camlob$R8total)

## Just exploratory analysis of how to use time here.
cam$sma_aveTemp2 = sma(cam$aveTemp, h = 2)$fitted
cam$sma_aveTemp3 = sma(cam$aveTemp, h = 3)$fitted
cam$sma_aveTemp4 = sma(cam$aveTemp, h = 4)$fitted

## Removing the nan variables
cam <- cam[!is.na(cam$rp),]
cam[cam$rp == 0]
## Adding 0.00001 just in case, so that we do not have 0 values and we can run
# the bowplot thingy
cam$rp <- cam$rp + 0.00001
cam$date <- as.Date(paste("1", cam$week, cam$year, sep = "-"), format = "%w-%W-%Y")

 #Strength
 str(cam)
 #Summary of Weather data
 summary(Camlob[,3:7])

 #Number of slaughtered floks
 colSums(Camlob[8])
 ############################################
 ## Box Cox transform
 par(mfrow=c(2,1),mgp=c(2,0.7,0),mar=c(3,3,1,1))

 boxcox(rp ~ aveTemp + maxTemp+ relHum+sunHours+precip,data=cam, plotit=TRUE)  
 #Explore a transformation on the response
 cam$rp_trans <-sqrt(cam$rp)
 boxcox(rp_trans ~ aveTemp + maxTemp+ relHum+sunHours+precip,data=cam, plotit=TRUE)  
 #Explore a transformation on the response
 # Could we transform the variables using the time information ? 
 ## Histogram
 hist(cam$rp)
 hist(cam$rp_trans)
 
 ## Some little trial of time analysis
 #Graphical description
 library(car)#relHum+sunHours+precip, 
 scatterplotMatrix( ~ rp_trans + aveTemp + sma_aveTemp2 + sma_aveTemp3 + sma_aveTemp4, 
                   diagonal= "boxplot", data = cam)
 #ratio seems to be more depending on temperature and not so much on humidity and sunshine

 #time series
 par(mfrow=c(3,2), mgp = c(2,0.7,0), mar = c(3,3,1,1))
 plot(cam$date,cam$relpos,type="l");plot(cam$date,cam$aveTemp,type="l");
 plot(cam$date,cam$maxTemp,type="l");
 plot(cam$date,cam$relHum,type="l")
 plot(cam$date,cam$sunHours,type="l");plot(cam$date,cam$precip,type="l")
 
 #We remove the 0 values
 cam$relHum[which(cam$relHum==0)]<-NA
 cam$sunHours[which(cam$sunHours==0)]<-NA

 common_columns = c("year","week","rp","rp_trans", "aveTemp", "maxTemp", 
                    "relHum", "sunHours","precip")
 ## Outliers removing 
 cam <- cam[-c(488),]
 cam <- cam[-c(437),]
 cam <- cam[-c(436),]
 cam <- cam[-c(382),]
 cam <- cam[-c(331),]
 
 cam <- cam[,common_columns]
 cam<-na.omit(cam)  # Omit nas !! Some func do not work with them
 
 #Making a Generalized Additive Model
 library(nlme)
 library(mgcv)
 par(mfrow = c(3,2), mgp = c(2,0.7,0), mar = c(3,3,1,1))
 
 cam$maxTemprelHum =   cam$maxTemp* cam$relHum
 gam_model = gam(rp_trans ~ s(aveTemp) + s(maxTemp)+ s(relHum)+
                   s(sunHours)+s(precip), data = cam)
 res = mean(residuals(gam_model)^2)
 plot(gam_model)
 caca = summary(gam_model)
 xtable(summary(gam_model)$p.table)
 xtable(summary(gam_model)$s.table)
 par(mfrow = c(3,2), mgp = c(2,0.7,0), mar = c(3,3,1,1))
 plot(gam(rp_trans ~ s(aveTemp) + s(maxTemp)+ s(relHum)+
            s(sunHours)+s(precip), data = cam))
 
 par(mfrow=c(1,3))
 plot(gam_model$fitted.values, gam_model$residuals)
 plot(gam_model$fitted.values, gam_model$residuals/sd(gam_model$residuals))
 qqnorm(gam_model$residuals)
 
 par(mfrow=c(2, 3))

 gam_model$residuals <- gam_model$residuals/sd(gam_model$residuals)
 plot(gam_model$residuals~cam$aveTemp, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)
 plot(gam_model$residuals~cam$maxTemp, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)
 plot(gam_model$residuals~cam$relHum, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)
 plot(gam_model$residuals~cam$sunHours, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)
 plot(gam_model$residuals~cam$precip, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)
 plot(gam_model$residuals~cam$week, data = cam, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)
 

 #Average temperature shows Piecewise linear
 # max temp shows some polynomial trend and .
 # realative humidity is almost linear
 #We try to describe the data with a piecewise linear

 #Making a model Linear Model with all interactions with removed NA weeks
 #cam<-cam[,1:8]
 #cam<-na.omit(cam)
 
 
 ######################### MODELS  ##############################

 ### Adding extra variables ! 

 ### Squaring variables
 cam$aveTemp_2 <-  cam$aveTemp^2 
 cam$maxTemp_2 <-  cam$maxTemp^2
 Model_init<-lm(rp_trans ~aveTemp + maxTemp+ relHum+sunHours+precip + 
                  aveTemp_2 + maxTemp_2,data=cam)
 #Test QQ and uniform Variance, CHECK!

 plot(Model_init,which=1:4)
 
 summary(Model_init)
 xtable(summary(Model_init))
 
 null <- lm(rp_trans~1, data = cam)
 Model_initRed<-step(Model_init, scope=list(lower=null, upper=Model_init), 
                     direction="both")
 
 par(mfrow=c(2,2))
 plot(Model_initRed,which=1:4)
 summary(Model_initRed)
 xtable(summary(Model_initRed))

 #Testing for linearity, Not so nessecary when we did the MAR

 par(mfrow=c(2, 3))
 plot(rstandard(Model_initRed)~aveTemp, data = cam, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)
 plot(rstandard(Model_initRed)~maxTemp, data = cam, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)
 plot(rstandard(Model_initRed)~relHum, data = cam, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)
 plot(rstandard(Model_initRed)~sunHours, data = cam, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)
 plot(rstandard(Model_initRed)~precip, data = cam, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)

 #Test for Independent Residuals, check!
 # par(mfrow=c(1,1))
 plot(rstandard(Model_initRed)~week, data = cam, ylim = c(-3,3), ylab="Standardized ...
residuals"); abline(h = 0); abline(h=c(-1, 1) * qt(.975, df = 419), lty = 2)

##### PIECEWISE LINEAR MODEL !!!!
 
 ## Piecewise linear model helper function
 pwl<-function(x,x0){
   ## x is data
   ## x0 is cut off
   ## The associated estimated parameter is for x > x0
   return( (x > x0) * (x-x0) )
 }
 
 ## This was not the optimal split point So finding it:
 #optimize optimizes the function (finds the minimum)
 #for values of zz (the split point) between 3 and 8
 optim<-optimize(function(aveTemp_sp){
   model = lm(rp_trans ~ aveTemp + pwl(aveTemp,aveTemp_sp) +  
                maxTemp +  # pwl(maxTemp, maxTemp_sp) + #+   maxTemp + 
                # aveTemp_2 + maxTemp_2 +
                relHum+sunHours+precip, data = cam)
   null <- lm(rp_trans~1, data = cam)
   model_red<-step(model, scope=list(lower=null, upper=model), direction="both")
   sum( residuals(model_red)^2 )
 },c(3,8))
 (x0.opt<-optim$minimum)
 
 
 aveTemp_pwl_value = 7.785
 maxTemp_pwl_value = 8.00
 cam$aveTemp_pwl<-(cam$aveTemp > aveTemp_pwl_value)*(cam$aveTemp-aveTemp_pwl_value)
 cam$maxTemp_pwl<-(cam$maxTemp > maxTemp_pwl_value)*(cam$maxTemp-maxTemp_pwl_value)
 
 cam$aveTemp_pwl_2<- cam$aveTemp_pwl^2
 cam$maxTemp_pwl_2<- cam$maxTemp_pwl ^2
 
 Model_init = lm(rp_trans ~  aveTemp + aveTemp_pwl  +  
                   maxTemp + maxTemp_pwl + #    maxTemp + 
              relHum + sunHours + precip +
                  aveTemp_2 + maxTemp_2 + aveTemp_pwl_2 + maxTemp_pwl_2 
              # aveTemp*aveTemp_pwl + relHum*aveTemp_pwl + precip*aveTemp_pwl
              ,data = cam)
 null <- lm(rp_trans~1, data = cam)
 Model_initRed<-step(Model_init, scope=list(lower=null, upper=Model_init), 
                     direction="both")
 
# drop1(Model_init)
#  Model_initRed <- Model_init
#  drop1(Model_init)
#  drop1(Model_init)
#  drop1(Model_init)
 summary(Model_init)
 
 par(mfrow=c(2,2))
 plot(Model_initRed,which=1:4)
 summary(Model_initRed)
 xtable(summary(Model_initRed))
 
 #Plotting for means model
 # Prediciton in terms od the aveTemp
 par(mfrow=c(2,3))
 
 ##################################  1  ##################################
 aveTemp_seg = seq(from=min(cam$aveTemp), to=max(cam$aveTemp), length.out=500)
 maxTemp_seg =seq(from=min(cam$maxTemp), to=max(cam$maxTemp), length.out=500)
 
 ave_Temps = c(4,6,8,9,10,12)
 for (ave_i in ave_Temps){
   newData <- data.frame(
         #aveTemp=aveTemp_seg
         maxTemp=maxTemp_seg,
       # maxTemp=mean(cam$maxTemp), 
        aveTemp=ave_i,precip=mean(cam$precip),
       relHum=mean(cam$relHum), sunHours=mean(cam$sunHours)) 
                         
   newData$maxTemp_pwl <-  (mean(newData$maxTemp) > maxTemp_pwl_value)*
     (mean(newData$maxTemp)-maxTemp_pwl_value)
   newData$aveTemp_pwl <-  (mean(newData$aveTemp) > aveTemp_pwl_value)*
     (mean(newData$aveTemp)-aveTemp_pwl_value)
   newData$maxTemp_pwl_2 <-  newData$maxTemp_pwl^2
   newData$aveTemp_pwl_2 <-  newData$aveTemp_pwl^2
   newData$maxTemp_2= newData$maxTemp^2
   newData$aveTemp_2 = newData$aveTemp^2
      
   Pred.ci <- predict(Model_initRed, newdata=newData, interval="confidence",level=.95)
   Pred.pi <- predict(Model_initRed, newdata=newData,interval="prediction",level=.95)
 
   #Plot Data and CI and PI
   subset = cam
   subset = cam[cam$aveTemp < ave_i +1,]
   subset = subset[subset$aveTemp > ave_i - 1,]
   
   ## We plot the originial rp and the square of the predictions
   plot(rp ~ maxTemp, data = subset, pch = 20, las = 1,ylim=c(0,1), 
        main = paste("aveTemp: ", as.character(ave_i), sep = "") )
   matlines(newData$maxTemp, Pred.ci^2, lty=c(1,2,2), col=c(1,3,3))
   matlines(newData$maxTemp, Pred.pi^2, lty=c(1,2,2), col=c(1,2,2))
 }
 
 # As we can see, around the average we have the best predictions
 ####################################################################
 #Looking at regions ratio distribution
