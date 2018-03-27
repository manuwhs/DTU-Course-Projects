#install.packages("psych")
library(psych)
#install.packages("ppcor")
library(ppcor)

plot_timeSeries <- function(X,Y, name, xlabel, ylabel, lwd, gp){
  # gp: If set to 1, then it generates a physical image
  if (gp == 1){
    png(file = name,width=800, height = 600, res=130) 
  }
  # Size and aspect ratio could be width=400,height=350,res=45
  plot(X, Y, 
       type = "l",               # Draw as a line
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

plot_ccfpccf <- function(X1,X2, name, Nlag, gp){
  # gp: If set to 1, then it generates a physical image
  if (gp == 1){
    png(file = name,width=800, height = 600, res=130) 
  }
  par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(2,0.7,0))
  ccf(X1 , X2, lag.max = Nlag )
  
if (0){
  ppcf_list = list()
  stat_list = list()
  lag_i = -Nlag:Nlag
  aux_i = 1:length(lag_i)
  for(i in aux_i){
    d = lag_i[i]
    if (d < 0 ){
      d = -d
      ppcfi= pcor(cbind(X1[1:(Ntr-d -1)],X2[(d+1):(Ntr-1)]))
    }
    else {
      ppcfi= pcor(cbind(X1[(d+1):(Ntr-1)],X2[1:(Ntr-d -1)]))
      
    }
    ppcf_list[[i]] <-ppcfi$estimate[1,2]
    stat_list[[i]] <-ppcfi$p.value[1,2]
  }
  plot(lag_i, ppcf_list, 
       type = "h",                     # Draw as a line
       lwd= 1,                         # Line width
       main="PCCF", # Title of the graph
       xlab="lag",            # x label
       ylab="PCCF",ylim = c(-0.2,0.9))
  
  lines(lag_i, stat_list, 
        type = "p",                     # Draw as a line
        lwd= 1)
}
  lag_i = -Nlag:Nlag
  mierda = pacf(cbind(X1,X2), plot = FALSE, lag.max = Nlag)
  ppcf = c(t(mierda$acf[,2,1]), t(mierda$acf[,1,2]))
  lag_i = c(t(mierda$lag[,2,1]), t(mierda$lag[,1,2]))
  plot(lag_i, ppcf, 
       type = "h",                     # Draw as a line
       lwd= 1,                         # Line width
       main="PCCF", # Title of the graph
       xlab="lag")
  
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

       