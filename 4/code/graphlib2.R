
plot_timeSeries <- function(X,Y, name, xlabel, ylabel, type ,lwd , col,xlim,ylim,gp ){
  # gp: If set to 1, then it generates a physical image
  if (gp == 1){
    png(file = name,width=800, height = 600, res=130) 
  }
  # Size and aspect ratio could be width=400,height=350,res=45
  plot(X, Y, 
       type = type,               # Draw as a line
       lwd= lwd,                         # Line width
       main= name, # Title of the graph
       xlab=xlabel,            # x label
       ylab=ylabel,
col=col,
xlim=xlim,
ylim=ylim
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
  pccf(X1, X2, lag.max = Nlag )
  
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

       