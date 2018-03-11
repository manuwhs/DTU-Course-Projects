setwd("C:/Users/Paulina/Desktop")
spr <- read.delim("C:/Users/Paulina/Desktop/SPR.txt")

#Summary statistics
str(spr)
head(spr)
summary(spr)

library(MASS)  
#attach package
boxcox(Response~.,data=spr, plotit=TRUE)  
#Explore a transformation on the response
spr$Response = sqrt(spr$Response)

#Plot the distribution of calcium ions and detergents
par(mfrow=c(2,1),mgp=c(2,0.7,0),mar=c(3,3,1,1))
boxplot(data=spr, Response~CaStock*DetStock, xlab = 'Calcium ions and detergents',ylab = 'Protein [RU]', las = 1)
stripchart(Response~CaStock*DetStock,data=spr, vertical=TRUE, method="jitter", xlab='Calcium ions and detergents', ylab="Protein [RU]",cex=1.2,pch=16, las=1,)

#Make a model
Model1a<-lm(Response~CaStock*DetStock,data=spr)
summary(Model1a)
#Eliminate the least significant parameter
drop1(Model1a,test="F")

#Reduced model (interaction not significant)
Model1b<-lm(Response~CaStock+DetStock,data=spr)
summary(Model1b)
drop1(Model1b,test="F")

#Reduced model (calcium ions not significant)
Model1c<-lm(Response~DetStock,data=spr)
summary(Model1c)
#Model is as simple as it can be
#test for assumptions
par(mfrow=c(2,2))
plot(Model1b,which=1:4)

##################################################################################
#Enzyme concentration distribution
par(mfrow=c(2,1))
boxplot(Response~EnzymeConc,data=spr, xlab = 'Enzyme concentration',ylab = 'Catalytic activity [RU]', las = 1)

#Enzyme concentration+calcium ions+detergent distribution
stripchart(Response~CaStock*DetStock*EnzymeConc,data=spr, vertical=TRUE, method="jitter", 
           xlab='Calcium ions, detergent and enzyme concentration', ylab="Catalytic activity [RU]",cex=1.2,pch=16, 
           las=1,col=c(2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5))

#Add enzyme concentration to the previous model
null <- lm(Response~1, data = spr)
Model2a<-lm(Response~CaStock*DetStock*EnzymeConc,data=spr)
Model2b<-step(Model2a, scope=list(lower=null, upper=Model2a), direction="both")
summary(Model2a)
summary(Model2b)

##############################################################################################################
#Are there any differences in performance among the enzymes in this study regarding the factors mentioned above?
#atach package
library(MASS) 
#Explore a transformation on the response
boxcox(Response~.,data=spr, plotit=TRUE) 
spr$Response = sqrt(spr$Response)
par(mfrow=c(2,1),mgp=c(2,0.7,0),mar=c(3,3,1,1))

#plot to check how enzymes itself influeces the response
boxplot(Response~Enzyme,data=spr,ylab="Cataclytic Activity [RU]", xlab = "Enzyme") 

#Interactio plot to check interaction between enzymes and concentration
interaction.plot( spr$EnzymeConc,spr$Enzyme,spr$Response, type=c("b"), ylab = "Cataclytic Activity [RU]", las=1, trace.label = "Block", xlab = "Enzyme Concentration")

#Plot of the enzymes and concentrations
boxplot(Response~Enzyme + EnzymeConc, data=spr,ylab="Cataclytic Activity [RU]", xlab = "Enzyme and concentration", col=c(8,7,6,5,4,8,7,6,5,4,8,7,6,5,4,8,7,6,5,4))

#Linear models without interactions and t-tests
#to check significance of enzymes itself
ModelEE<-lm(Response~Enzyme, data = spr)
summary(ModelEE)
ModelEC<-lm(Response~Enzyme+EnzymeConc, data=spr)
summary(ModelEC)
ModelED<-lm(Response~Enzyme+DetStock, data=spr)
summary(ModelED)
ModelECal<-lm(Response~Enzyme+CaStock, data=spr)
summary(ModelECal)

#Linear models with interactions and t-tests
ModelE1<-lm(Response~Enzyme*EnzymeConc, data=spr)
summary(ModelE1)
ModelE2<-lm(Response~Enzyme*DetStock, data=spr)
summary(ModelE2)
ModelE3<-lm(Response~Enzyme*CaStock, data=spr)
summary(ModelE3)

## Checking a possible boxcox transformaiton of the data Response.
par(mfrow=c(2,1),mgp=c(2,0.7,0),mar=c(3,3,1,1))
library(MASS)  #atach package
boxcox(Response~.,data=spr, plotit=TRUE)  #Explore a transformation on the response
spr$Response = sqrt(spr$Response)
boxcox(Response~.,data=spr, plotit=TRUE)  #Explore a transformation on the response
## Checking for systematic errors
spr$RunDate = as.factor(spr$RunDate)
spr_0 = spr[spr[, "EnzymeConc"] == 0.0,]
boxplot(data=spr_0, Response~RunDate, xlab = 'Ca and Detergent',ylab = 'Protein 
        [RU]', las = 1,main="Distribution of Observations")
fit <- aov(Response ~ RunDate, data=spr_0)
fit2 <- lm(Response ~ RunDate, data=spr_0)
anova(fit2)

## Building a final model
spr <- read.delim("./SPR.txt")
hist(spr$Response)
spr$Response = sqrt(spr$Response)
hist(spr$Response)
spr <- read.delim("./SPR.txt")
spr$Response = sqrt(spr$Response)
#Initial model
Model_F1<-lm(Response~DetStock*EnzymeConc*Enzyme,data=spr)
summary(Model_F1)
Model_F1R<-step(Model_F1, scope=list(lower=null, upper=Model_F1), direction="both")
summary(Model_F1R)
AIC(Model_F1R)

par(mfrow=c(2,2))
plot(Model_F1R,which=1:4)

#Testing model assumptions
par(mfrow=c(1, 1))
plot(rstandard(Model_F1R) ~ EnzymeConc, data = spr, ylim = c(-3,3), ylab="Standardized ...
     residuals")
abline(h = 0)
abline(h=c(-1, 1) * qt(.975, df = 140), lty = 2)

# Remove outlier
spr <- spr[-c(160),]
spr <- spr[-c(147),]

#Sqrt transformation
Model_F1<-lm(Response~DetStock*sqrt(EnzymeConc)*Enzyme,data=spr)
summary(Model_F1)
Model_F1R<-step(Model_F1, scope=list(lower=null, upper=Model_F1), direction="both", trace=1)
summary(Model_F1R)
AIC(Model_F1R)
par(mfrow=c(2,2))
plot(Model_F1R,which=1:4)

#Testing model assumptions
par(mfrow=c(1, 2))
plot(rstandard(Model_F1R) ~ DetStock, data = spr, ylim = c(-3,3), ylab="Standardized ...
     residuals")
abline(h = 0)
abline(h=c(-1, 1) * qt(.975, df = 140), lty = 2)

plot(rstandard(Model_F1R) ~ Enzyme, data = spr, ylim = c(-3,3), ylab="Standardized ...
     residuals")
abline(h = 0)
abline(h=c(-1, 1) * qt(.975, df = 140), lty = 2)

par(mfrow=c(2,2))
plot(Model_F1R,which=1:4)
which(duplicated(spr))

