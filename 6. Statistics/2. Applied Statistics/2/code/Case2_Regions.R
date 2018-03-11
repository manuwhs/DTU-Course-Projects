#Open file 
Campylobacter <- read.delim("case2regionsOnePerBatch.txt")
#Ratios

cam<-Campylobacter;cam<-cam[,-(8:25)]
cam$relpos<-(Campylobacter$pos/Campylobacter$total);
cam$r1relpos<-(Campylobacter$R1pos/Campylobacter$R1total);
cam$r2relpos<-(Campylobacter$R2pos/Campylobacter$R2total)
cam$r3relpos<-(Campylobacter$R3pos/Campylobacter$R3total);
cam$r4relpos<-(Campylobacter$R4pos/Campylobacter$R4total);
cam$r5relpos<-(Campylobacter$R5pos/Campylobacter$R5total)
cam$r6relpos<-(Campylobacter$R6pos/Campylobacter$R6total);
cam$r7relpos<-(Campylobacter$R7pos/Campylobacter$R7total);
cam$r8relpos<-(Campylobacter$R8pos/Campylobacter$R8total)
cam$date <- as.Date(paste("1", cam$week, cam$year, sep = "-"), format = "%w-%W-%Y")

#Data transformation
cam$relpos = sqrt(cam$relpos)
camall$relpos = sqrt(camall$relpos)

#Tree model
library(tree)
TREE<-tree(relpos~+aveTemp+maxTemp+relHum+sunHours+precip,data=cam1)
plot(TREE)
text(TREE)
summary(TREE)
###############################################################
##Regions differences
#New dataset
camf<-cam[,9:17]
camf$year<-cam$year
camf$week<-cam$week

#Make regions a factor with 8 levels
library(reshape)
camf2<-melt(camf, id=c("year","week","date"))
camf2$aveTemp<-cam$aveTemp
camf2$maxTemp<-cam$maxTemp
camf2$relHum<-cam$relHum
camf2$sunHours<-cam$sunHours
camf2$precip<-cam$precip
camf2$relpos<-cam$relpos

names(camf2)<-c("year", "week", "date", "region", "ratio", "aveTemp",
                "maxTemp", "relHum", "sunHours", "precip", "relpos")
camf2$region = as.factor(camf2$region)

camf4<-Campylobacter[,1:2]
camf4$R1total<-Campylobacter$R1total
camf4$R2total<-Campylobacter$R2total 
camf4$R3total<-Campylobacter$R3total
camf4$R4total<-Campylobacter$R4total
camf4$R5total<-Campylobacter$R5total
camf4$R6total<-Campylobacter$R6total
camf4$R7total<-Campylobacter$R7total
camf4$R8total<-Campylobacter$R8total
camf4$date<-cam$date

camf3<-melt(camf4, id=c("year", "week", "date"))
names(camf3)<-c("year", "week", "date", "Region", "totvalue")
camf3$aveTemp <-camf2$aveTemp
camf3$maxTemp <-camf2$maxTemp
camf3$relHum <-camf2$relHum
camf3$sunHours <-camf2$sunHours
camf3$precip <-camf2$precip
camf3$region<-camf2$region
camf3$ratio <-camf2$ratio


## Omit the Nans and also the measurements with less than 5 flocks
# camf22<-na.omit(camf3) 
camf22 <- camf22[camf22$totvalue > 5,]

#Linear model with only 1 variable
model<-lm(sqrt(ratio)~region,data=camf2)
library(xtable)
xtable(summary(model))

## Boxplot 
par(mfrow = c(1,1), mgp = c(2,0.7,0), mar = c(3,3,1,1))
camf22$rp_trans <- sqrt(camf22$ratio)
boxplot(data=camf22, rp_trans~region, xlab = 'Region',ylab = 'sqrt(ratio)', las = 1)

##Testing variances
library(xtable)
xtable(anova(model))

#Linear model with all of the variables
model2<-lm(sqrt(ratio)~region+aveTemp+maxTemp+precip, data=camf22)
library(xtable)
xtable(summary(model2))
model2b <- step(lm(sqrt(ratio)~region+aveTemp+maxTemp+precip, data=camf22))
xtable(summary(model2b))

boxplot(residuals(model2b)~camf22$region, xlab = 'Region',ylab = 'Residual', las = 1)
summary(model2b)

## The final model
camf22$aveTemp_2<- camf22$aveTemp^2
camf22$maxTemp_2<- camf22$maxTemp^2
aveTemp_pwl_value = 7.785
maxTemp_pwl_value = 8.00
camf22$aveTemp_pwl<-(camf22$aveTemp > aveTemp_pwl_value)*
  (camf22$aveTemp-aveTemp_pwl_value)
camf22$maxTemp_pwl<-(camf22$maxTemp > maxTemp_pwl_value)*
  (camf22$maxTemp-maxTemp_pwl_value)

camf22$aveTemp_pwl_2<- camf22$aveTemp_pwl^2
camf22$maxTemp_pwl_2<- camf22$maxTemp_pwl ^2

Model_init = lm(sqrt(ratio) ~  aveTemp + aveTemp_pwl  +  
                  maxTemp + maxTemp_pwl + #    maxTemp + 
                  relHum  +
                  aveTemp_2 + maxTemp_2 + aveTemp_pwl_2 + maxTemp_pwl_2 
                  + region 
                # aveTemp*aveTemp_pwl + relHum*aveTemp_pwl + precip*aveTemp_pwl
                ,data = camf22)
null <- lm(sqrt(ratio)~1, data = camf22)
Model_initRed<-step(Model_init, scope=list(lower=null, upper=Model_init), direction="both")
summary(Model_initRed)

