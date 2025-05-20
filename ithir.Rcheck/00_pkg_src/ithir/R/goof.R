# Purpose        : Goodness of fit 
# Maintainer     : Brendan Malone (brendan.malone@csiro.au); 
# Contributions  : Budiman Minasny, Daniel Saurette 
# Status         : working
# Note           : 
#



goof <- function(observed,predicted, plot.it = FALSE, type="DSM"){
  # Coefficient of determination
  # Step 1: calculate the mean of y variable, y_bar
  ybar<- mean(predicted) 
  # Step 2: calculate the regression coefficents of y on x
  xbar<- mean(observed) 
  slope<- sum((observed-xbar)*(predicted-ybar))/sum((observed-xbar)^2)
  intercept<- ybar-(slope*xbar)
  # Step 3: Now calculate y_hat
  yhat<- (observed*slope)+intercept
  # Step 4: Calculate the Sums of Squares
  SStot<- sum((predicted - ybar)^2) 
  SSres<- sum((yhat - predicted)^2) 
  #new_SSreg<- new_SStot-new_SSres
  
  R2<- 1-(SSres/SStot)


  # Standard error of prediction ^2
  SEP2 <- mean((observed - predicted)^2)
  
  # Standard error of prediction
  SEP <- sqrt(SEP2)
  
  #Bias
  bias <- mean(predicted) - mean(observed)
  
  # residual  variance
  SEP2c <- sum(((predicted - bias - observed)^2) / length(observed))
  SEPc <- sqrt(SEP2c)
  
  # ratio of performance to deviation
  RPD <- sd(observed) / SEP
  
  # Ratio of performance to interquartile distance
  IQ <- c(quantile(observed))[4] - c(quantile(observed))[2]
  RPIQ <- IQ / SEP
  
  # Concordance
  mx <- mean(observed)
  my <- mean(predicted)
  s2x <- var(observed)
  s2y <- var(predicted)
  sxy <- mean((observed-mx) * (predicted-my))
  ccc <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  
  if (plot.it==TRUE){MASS::eqscplot(observed, predicted)
  abline(a = 0, b = 1, col = "brown4")}
  
  if (type == "DSM"){ gf <- data.frame(R2=R2, concordance=ccc, MSE=SEP2, RMSE=SEP, bias=bias, row.names=NULL)}
  else if (type == "spec"){ gf <- data.frame(R2=R2, concordance=ccc, MSE=SEP2, RMSE=SEP, bias=bias, 
                                       MSEc=SEP2c,RMSEc=SEPc, RPD=RPD, RPIQ=RPIQ, row.names=NULL)}
  else {stop("ERROR: Revise the type of output you require. Select from either DSM or spec")} 
  
  return(gf)
}

#END SCRIPT