#install.packages("vars")
library(vars)

get_likelihood <- function(y, ypred, SigmaYYpredList){
  # This function gets the likelihood function of the model
  Sigma2inv = solve(Sigma2)
  
  logpred = 0
  Nt = dim(y)[1]
  for (i in 1:Nt){
    logpred = logpred + log(det(SigmaYYpredList[[i]])) + t(y[i,] - ypred[i,])%*%solve(SigmaYYpredList[[i]])%*%(y[i,] - ypred[i,])
  }
  logpred = - logpred/2
  return(logpred)
}

##signTest
sign.test = function(res) {
  binom.test(sum(res[-1]*res[-length(res)]<0), length(res)-1)
}

## Wrapping all tests in one function
vars.test <- function(x){
  norm <- sapply(normality.test(x)[[2]],function(xx) xx$p.value)
  arch <- arch.test(x)[[2]]$p.value
  ser <- serial.test(x)[[2]]$p.value
  ## Return vector of p values
  return(c(norm,"arch"=arch,"serial"=ser))
}

alltest <- function(residuals){
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
  p_value_bintest = bt$p.value
  int_conf = bt$conf.int
  
  return (list(p_value_bintest))
}