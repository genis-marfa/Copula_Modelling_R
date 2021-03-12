# Packages for the whole project
library(normtest)
library(qrmdata)
library(qrmtools)
library(QRM)
library(MASS)
library(xts)
library(copula)

# Importing data and converting into xts format 
data5Y  <- read.delim("E:/LSE/St 429/data5Y.txt")
dat     <- data.frame(data5Y)
rawdata <- dat[dat[,2]!=0&dat[,3]!=0&dat[,4]!=0&dat[,5]!=0&dat[,6]!=0&
               dat[,7]!=0&dat[,8]!=0&dat[,9]!=0&dat[,10]!=0&dat[,11]!=0,]
date    <- as.Date(as.character(rawdata[, 1]),"%d/%m/%Y")
stocks  <- xts(rawdata[, -1], date)

# Individual stocks
GE   <- stocks[, 1]
NEP  <- stocks[, 2]
FAN  <- stocks[, 3]
ENPH <- stocks[, 4]
FSLR <- stocks[, 5]
FTEK <- stocks[, 6]
RUN  <- stocks[, 7]
PZD  <- stocks[, 8]
FCEL <- stocks[, 9]
TSLA <- stocks[, 10]

# Plotting a stock return
plot(zoo(stocks), xlab = "Time", main = "Stock prices", col = c("darkblue", 
  "darkgreen", "red", "orange", "violet", "black", "gold", "grey", "brown", "pink"))

# Creating the portfolio and log-returns
portfolio    <- stocks
log_s        <- log(portfolio)
tt           <- nrow(stocks)
lr_portfolio <- diff(log_s)[2:tt, ] #Log returns on portfolio
lr_GE        <- lr_portfolio[, 1]
lr_NEP       <- lr_portfolio[, 2]
lr_FAN       <- lr_portfolio[, 3]
lr_ENPH      <- lr_portfolio[, 4]
lr_FSLR      <- lr_portfolio[, 5]
lr_FTEK      <- lr_portfolio[, 6]
lr_RUN       <- lr_portfolio[, 7]
lr_PZD       <- lr_portfolio[, 8]
lr_FCEL      <- lr_portfolio[, 9]
lr_TSLA      <- lr_portfolio[, 10]

Avg_lr            <- cbind(mean(lr_GE), mean(lr_NEP), mean(lr_FAN), mean(lr_ENPH), mean(lr_FSLR),
                   mean(lr_FTEK), mean(lr_RUN), mean(lr_PZD), mean(lr_FCEL), mean(lr_TSLA))
colnames(Avg_lr) <- c("GE", "NEP", "FAN", "ENPH", "FSLR", "FTEK", "RUN", "PZD", "FCEL", "TSLA")

plot(zoo(lr_portfolio), xlab = "Time", main = "Log returns", col = c("darkblue", 
  "darkgreen", "red", "orange", "violet", "black", "gold", "grey", "brown", "pink"))

init_investment <- rep(1000, 10) # initial investment
price_portfolio <- stocks[1, ] # initial price of stock
No_shares <- int_investment/price_portfolio # number of shares = Lambda
X <- diff(log_s)[2:tt, ] # risk factor changes
CV_portfolio <- sum(as.matrix(stocks[tt, ])*No_shares) # current value of portfolio {current time = end time}
weights_portfolio <- as.matrix(No_shares*as.matrix(stocks[tt, ])/CV_portfolio) # weights of each share in the portfolio


loss.sim <- function(Xval, proportion, value){ # Function to create the losses for the portfolio 
  # arguments:
  # Xval ... matrix or vector of d risk-factor changes 
  # proportion ... row vector of d weights representing the proportion invested in each stock 
  # value ... scalar representing the value of the portfolio
  
  if (is.matrix(Xval)){
    prod <- (exp(Xval)-1) %*% t(proportion) 
  } else {
    n <- length(Xval)
    prod <- proportion * (exp(Xval)-1)
  }
  loss <- -value * prod
  return(loss)
}

# Historical losses for the portfolio:
Loss_sim           <- loss.sim(as.matrix(X), weights_portfolio, CV_portfolio) 
colnames(Loss_sim) <- c("Historical losses")
summary(Loss_sim)
Exp_loss           <- mean(Loss_sim,na.rm=TRUE) # Expected losses
hist(Loss_sim, breaks = 50, xlab = "Historical losses", 
     main = "Histogram of Empirical losses",col = "light blue")

# Normalizing loss variable: Function to normalize the data
normalize <- function(x) { 
  return ((x - min(x)) / (max(x) - min(x)))
}
Loss_sim_standard <- normalize(Loss_sim)

# Lineralizing losses
losslin.sim <- function(Xval, proportion, value){
  if (is.matrix(Xval)){
    n <- dim(Xval)[1]
    prod <- (Xval) %*% t(proportion)
  } else {
    n <- length(Xval)
    prod <- proportion * Xval
  }
  loss <- -value * prod
  return(loss)
}

# Historical linearized losses for the portfolio
losslin_sim           <- losslin.sim(X, weights_portfolio, CV_portfolio) 
colnames(losslin_sim) <- c("Linearized  Losses")
summary(losslin_sim)

hist(losslin_sim, breaks = 50,xlab = "Linearlized losses",
     main = "Histogram of linearlized losses", col = "grey")

combinelin_loss <- xts(cbind(losslin_sim,Loss_sim), order.by = date[2:tt])
plot.xts(combinelin_loss, xlab = "Time", ylab = "Price", screens = factor(1, 1),
         auto.legend = TRUE, col = c("darkblue","gold"),
         main = "Plot of linearlized losses and Historical losses") 
addLegend(legend.loc = "bottomleft", bty = "n", y.intersp = 1.2,lty = rep(1,2),
          col = c("darkblue","gold"),
          c("Linearized losses", "Historical losses"))


# Mean and variance of linearnized loss (Theoretical calculation)
muX.hat     <- colMeans(X,na.rm = TRUE)
sigmaX.hat  <- var(X,na.rm = TRUE)
meanLosslin <- -value_portfolio*sum(weights_portfolio*muX.hat)
varlosslin  <- value_portfolio^2*(weights_portfolio %*% sigmaX.hat %*% t(weights_portfolio))

# Fitting a distribution- normal distribution
set.seed(16)
#Finding the mean and variance of the fitted normal distribution:
fit_normal <- fitdistr(Loss_sim, densfun = "normal") 
# Plot histogram for the loss function to check where most values are present:
hist(Loss_sim, breaks = 100, main = "Histogram of historical losses",
     xlab = "Historical losses", col = "green") 
# Create a random normal distribution with similar mean and variance
fit_normal_data         <- rnorm(nrow(Loss_sim), fit_normal$estimate[1], fit_normal$estimate[2]) 
fit_normal_combineddata <- xts(cbind(Loss_sim, rnorm(nrow(Loss_sim), fit_normal$estimate[1], 
                           fit_normal$estimate[2])), order.by = date[2:tt])
colnames(fit_normal_combineddata) <- c("Historical losses", "Normal distribution")

# Plotting historical losses and random normal distribution
plot.xts(fit_normal_combineddata, xlab = "Time", ylab = "Price",
         screens = factor(1, 1), auto.legend = TRUE,
         main = "Plot of historical losses and random normal distribution" ) 
addLegend(legend.loc = "bottomleft", bty = "n", y.intersp = 1.2,
          lty = rep(1, 2), col = c("red","black"),
          c("Normal distribution","Historical losses"))
hist(fit_normal_combineddata[, 2], breaks = 20, 
     xlab = "Random normal distribution", 
     main = "Histogram of Random normal distribution",
     col = "blue") # Plotting histogram for the random normal distribution

# Q-Q plot & Jarque-Bera test (To check the normality of the log returns(stocks),loss values, portfolio)
qqplot(Loss_sim,fit_normal_data,xlab = "Historical losses",ylab = "Normal distribution",col="Darkblue")
title("QQ plot") #QQ plot with the random generated normal distribution with same mean and SD
abline(0, 1, col = "grey",lwd = 1)
norm_test_loss <- jb.norm.test(Loss_sim) # Jarque beta test for historical losses
norm_test_lrportfolio <- jb.norm.test(lr_portfolio) # Jarque-Bera test for portfolio

# Normal test for log returns of all stocks
norm_test_lrGE   <- jb.norm.test(lr_GE)
norm_test_lrNEP  <- jb.norm.test(lr_NEP)
norm_test_lrFAN  <- jb.norm.test(lr_FAN)
norm_test_lrENPH <- jb.norm.test(lr_ENPH)
norm_test_lrFSLR <- jb.norm.test(lr_FSLR)
norm_test_lrFTEK <- jb.norm.test(lr_FTEK)
norm_test_lrRUN  <- jb.norm.test(lr_RUN)
norm_test_lrPZD  <- jb.norm.test(lr_PZD)
norm_test_lrFCEL <- jb.norm.test(lr_FCEL)
norm_test_lrTSLA <- jb.norm.test(lr_TSLA)


# VaR and ES for Losses -- different levels of alpha -- empirical distribution
alpha <- c(seq(0.1,0.8,0.1), seq(0.9,1,0.01))
VaR.hs <- quantile(Loss_sim, alpha, na.rm = TRUE) # Calculating VaR for empirical distribution
ES.hs <- rep(0, length(alpha)) # To calculate true ES for any distribution
for(i in 1:length(alpha)) {
  values <- Loss_sim[Loss_sim > VaR.hs[i]]
  ES.hs[i] <- mean(values, na.rm = TRUE)
}
ran <- range(VaR.hs, ES.hs, na.rm = TRUE) # creating range for y-axis
plot(alpha, VaR.hs, type = "l", ylim = ran, xlab = expression(alpha), lwd = 2,
     ylab = expression("Estimated VaR & ES"[alpha])) # true ES_alpha
lines(alpha, ES.hs , type = "l", col = "maroon3", lwd = 2) # ES_alpha estimate
legend("topleft", bty = "n", y.intersp = 1.2, lty = rep(1, 2),
       col = c("black", "maroon3"), legend = c(expression(VaR[alpha]),
       expression(ES[alpha])))

## Compare them with empirical distribution & normal distribution

# Var and ES for normal distribution

# Calculating VaR for normal distribution:
var_n <- qnorm(alpha,mean = fit_normal$estimate[1], sd = fit_normal$estimate[2]) 
# Calculating ES for normal distribution:
ES_n <- ESnorm(alpha,mu = fit_normal$estimate[1], sd = fit_normal$estimate[2]) 

# VaR plot for empirical and normal distribution
plot(alpha, VaR.hs, type = "l", col ="green", xlab = expression(alpha), 
     ylab = expression(Var[alpha]), lwd = 2)
lines(alpha, var_n, type = "l", col = "red", lwd = 2)
legend("topleft", bty = "n", y.intersp = 1.2, lty = rep(1, 2),
       col= c("green", "red"), legend = c(expression(VaR[alpha]), 
       expression(VaR_normal[alpha])))

# Expected shortfall(ES) plot for empirical and normal distribution
plot(alpha, ES.hs, type = "l", col = "blue", ylim = ran,
     xlab = expression(alpha), ylab = expression(ES[alpha]), lwd = 2)
lines(alpha, ES_n, type = "l",col = "gold2",lwd = 2)
legend("topleft", bty = "n", y.intersp = 1.2, lty = rep(1, 2),
       col = c("blue", "gold2"), legend = c(expression(ES[alpha]),
       expression(ES_normal[alpha])))

# Fixing alpha = 0.8 and then calculating ES for the empirical distribution
alpha2       <- 0.8
VaR.hs_fixed <- quantile(Loss_sim, alpha2, na.rm = TRUE)
ES.hs_fixed  <- rep(0, length(alpha2))

for(i in 1:length(alpha2)) {
  values         <- Loss_sim[Loss_sim > VaR.hs_fixed[i]]
  ES.hs_fixed[i] <- mean(values, na.rm = TRUE)
}

ran <- range(VaR.hs, ES.hs, na.rm = TRUE)
plot(alpha, VaR.hs, type = "l", ylim = ran, xlab = expression(alpha),
     col = "violet",lwd = 2, main = expression("VaR at ES"[alpha=0.8]),
     ylab = expression("VaR "[alpha])) # VaR_alpha
points(alpha2, ES.hs_fixed, pch = 17, col = "green", cex = 2)
abline(ES.hs_fixed, 0, col = "grey", lwd = 0.5)
points(0.93, ES.hs_fixed, pch = 18,col = "blue",cex = 2)
abline(v = 0.93, col = "darkgrey", lwd = 0.5)
legend("topleft", bty = "n", y.intersp = 1.2, pch = c(17, 18),
       col = c("green", "blue"),cex = 2, legend = c(expression(ES[alpha=0.8]),
       expression(VaR[alpha=0.93])))

# Extract dates and remove from data frame.
Dates     <- strptime(dat$DATE, "%d/%m/%Y") 
dat$DATE  <- NULL # Remove dates from dataframe. 

# Convert dataframe to matrix
Prices     <- data.matrix(dat)
LogPrices  <- log(Prices)
LogReturns <- diff(LogPrices)

# Remove First Entry in dates vector:
#   This is done because when we calculate the log returns as the difference
#   in log prices, we do not produce a value at the first data point.
Dates    <- Dates[2:length(Dates)]

# Copula fitting:
FitCopulas <- function(stock1, stock2) {
  Dat <- cbind(LogReturns[, stock1], LogReturns[, stock2])
  colnames(Dat) <- c(stock1, stock2)
  
  # Delete rows where both stocks have zero return: 
  #     These correspond to holidays  (e.g. 25 December)
  Dat <- Dat[Dat[, 1]!=0 & Dat[,2] !=0, ]
  CopulaX <- pobs(as.matrix(Dat))
  
  par(mfrow = c(1,1))
  plot(CopulaX, xlab = stock1, ylab = stock2)
  par(mfrow = c(2,2))
  
  print("Sample version of rank correlations:")
  cat("Spearman rho: ", Spearman(CopulaX)[1,2], '\n')
  cat("Kendal tau: ", Kendall(CopulaX)[1,2], '\n')
  cat('\n')
  
  # -------------------------
  #       Gauss Copula:
  # -------------------------
  cop_model <- normalCopula(dim = 2)
  normFit <- fitCopula(cop_model, CopulaX, method = 'ml')
  
  cat("Gaussian Fit: ", '\n')
  cat("Rho: ", coef(normFit), '\n')
  
  rho    <- coef(normFit)
  cat("Calibrated paramter (Spearman): ", 2* sin(pi/6 * Spearman(CopulaX))[1,2], '\n')
  cat("Calibrated paramter (Kendall): ",  sin(pi/2 * Kendall(CopulaX))[1,2], '\n')
  
  # Overlay fitted copula contours:
  normal <- normalCopula(param = rho, dim = 2)
  plot(CopulaX, xlab="", ylab="", col = "light blue", yaxt="n", main = "Gaussian fit")
  par(new=TRUE)
  contour(normal, dCopula, xlab = stock1, ylab = stock2)
  
  # Calculate corresponding loglikelihood
  normLL <- loglikCopula(param = coef(normFit), u=CopulaX, copula = normal)
  cat("LogLikelihood max: ", normLL, '\n', '\n')
  
  # -------------------------
  #       t Copula:
  # -------------------------
  cop_model <- tCopula(dim = 2)
  tFit <- fitCopula(cop_model, CopulaX, method = 'ml')
  cat("tStudent Fit: ", '\n')
  cat("(Rho, ndof): ", coef(tFit), '\n')
  
  rho <- coef(tFit)[1]
  ndof <- coef(tFit)[2]
  cat("Calibrated paramter (Kendall): ",  sin(pi/2 * Kendall(CopulaX))[1,2], '\n')
  
  # Overlay fitted copula contours:
  tstu <- tCopula(param = rho,  df = ndof, dim = 2)
  plot(CopulaX, xlab="", ylab="", col = "light blue", yaxt="n", main = "tStudent fit")
  par(new=TRUE)
  contour(tstu, dCopula, xlab = stock1, ylab = stock2)
  
  # Calculate corresponding loglikelihood
  tstuLL <- loglikCopula(param = coef(tFit), u=CopulaX, copula = tstu)
  cat("LogLikelihood max: ", tstuLL, '\n', '\n')
  
  # -------------------------
  #    Archimidean Copulas:
  # -------------------------
  # Gumbel:
  cop_model <- gumbelCopula(dim = 2)
  gumbFit <- fitCopula(cop_model, CopulaX, method = 'ml')
  cat("Gumbel Fit: ",  '\n')
  cat("Theta: ", coef(gumbFit), '\n')
  
  theta <- coef(gumbFit)
  cat("Calibrated paramter (Kendall): ",  1/(1-Kendall(CopulaX))[1,2], '\n')
  
  # Overlay fitted copula contours:
  gumb <- gumbelCopula(param = theta, dim = 2)
  plot(CopulaX, xlab="", ylab="", col = "light blue", yaxt="n", main = "Gumbel fit")
  par(new=TRUE)
  contour(gumb, dCopula, xlab = stock1, ylab = stock2)
  
  # Calculate corresponding loglikelihood
  gumbLL <- loglikCopula(param = coef(gumbFit), u=CopulaX, copula = gumb)
  cat("LogLikelihood max: ", gumbLL, '\n', '\n')
  
  # Clayton:
  cop_model <- claytonCopula(dim = 2)
  clayFit <- fitCopula(cop_model, CopulaX, method = 'itau')
  
  cat("Clayton Fit: ", '\n')
  cat("Theta: ", coef(clayFit), '\n')
  
  theta <- coef(clayFit)
  cat("Calibrated paramter (Kendall): ", (2*Kendall(CopulaX)/(1-Kendall(CopulaX)))[1,2] , '\n')
  
  # Overlay fitted copula contours:
  clay <- claytonCopula(param = theta, dim = 2)
  plot(CopulaX, xlab="", ylab="", col = "light blue", yaxt="n", main = "Clayton fit")
  par(new=TRUE)
  contour(clay, dCopula, xlab = stock1, ylab = stock2)
  
  # Calculate corresponding loglikelihood
  clayLL <- loglikCopula(param = coef(clayFit), u=CopulaX, copula = clay)
  cat("LogLikelihood max: ", clayLL, '\n', '\n')
}

FitCopulas("GE", "TSLA")
FitCopulas("GE", "RUN")
FitCopulas("GE", "NEP")
FitCopulas("GE", "PZD")
FitCopulas("GE", "FCEL")
FitCopulas("GE", "FTEK")
FitCopulas("GE", "PZD")
FitCopulas("GE", "FAN")
FitCopulas("GE", "FSLR")
FitCopulas("NEP", "FAN")
FitCopulas("NEP", "TSLA")
FitCopulas("NEP", "ENPH")
FitCopulas("NEP", "FSLR")
FitCopulas("NEP", "FTEK")
FitCopulas("NEP", "RUN")
FitCopulas("NEP", "PZD")
FitCopulas("NEP", "FCEL")
FitCopulas("FAN", "PZD")
FitCopulas("FAN", "FSLR")
FitCopulas("FAN", "ENPH")
FitCopulas("FAN", "TSLA")
FitCopulas("FAN", "FTEK")
FitCopulas("FAN", "FCEL")
FitCopulas("FAN", "RUN")
FitCopulas("ENPH", "FSLR")
FitCopulas("ENPH", "PZD")
FitCopulas("ENPH", "TSLA")
FitCopulas("ENPH", "RUN")
FitCopulas("ENPH", "FTEK")
FitCopulas("ENPH", "FCEL")
FitCopulas("FSLR", "FTEK")
FitCopulas("FSLR", "RUN")
FitCopulas("FSLR", "PZD")
FitCopulas("FSLR", "FCEL")
FitCopulas("FSLR", "TSLA")
FitCopulas("FTEK", "FCEL")
FitCopulas("FTEK", "RUN")
FitCopulas("FTEK", "PZD")
FitCopulas("FTEK", "TSLA")
FitCopulas("RUN", "TSLA")
FitCopulas("RUN", "FCEL")
FitCopulas("RUN", "PZD")
FitCopulas("PZD", "FCEL")
FitCopulas("PZD", "TSLA")
FitCopulas("FCEL", "TSLA")

## Calculate VaR (& ES) by Monte-Carlo Simulation Method

Loss <- as.matrix(-X)
alpha3 <- c(seq(0.8, 0.99, 0.01))

# Make loss(log return) arrays
gridm <- function(stock){
  
  grid <- cbind(Loss[, stock[1]], Loss[, stock[2]])
  grid <- grid[grid[, 1] != 0 & grid[, 2] != 0, ] # get rid of non-trading days
  return(grid)
  
}

# (L1, L2) ~ C
VaRcalC <- function(train, alpha) { 
  # function for calculating VaR when using t-student copula(C)
  
  CopulaX <- apply(train, 2, edf, adjust = 1)
  tFit <- fitCopula(tCopula(dim = 2), CopulaX, method = 'ml') # best copula models chosen from question 4 are all t copula
  
  sumVaRC_sim <- matrix(0, nrow = 5000, ncol = length(alpha))
  for (i in 1:5000) { # repeat 5000 times
    U <- rCopula(1500, tCopula(param = coef(tFit)[1],
                               df = coef(tFit)[2])) # generate two random arrays compelling copula we chose in question 4
    L1 <- quantile(train[, 1], U[, 1])
    L2 <- quantile(train[, 2], U[, 2]) # get two loss samples with the copula relation above
    L <- L1 + L2 # get sum of losses sample
    sumVaRC_sim[i, ] <- quantile(L, alpha) # calculate VaR from sum of losses sample
  }
  
  sumVaRC <- apply(sumVaRC_sim, 2, mean) # use mean of simulated VaRs to approach the real VaR 
  table <- rbind(alpha, sumVaRC) # record result
  return(table)
  
}

# Backtesting 
Backtest <- function(test, alpha, sumVaRC) { 
  # function for Unconditional Coverage Test Likelihood Ratio
  
  n <- nrow(test) # total number of test sample losses
  failuretimes <- rep(0, length(alpha))
  LR_uc <- rep(0, length(alpha))
  Pvalue <- rep(0, length(alpha))
  for(k in 1:length(alpha)) {
    x <- sum(apply(test, 1, sum) >= sumVaRC[2, k]) # the number of losses exceeding VaR under specific alpha
    p <- 1-alpha[k] # significance level
    failuretimes[k] <- x # the number of failure times
    LR_uc[k] <- -2*log((1-p)^(n-x)*p^x)+2*log((1-x/n)^(n-x)*(x/n)^x) # Likelihood Ratio of Unconditional Coverage Test
    Pvalue[k] <- 1-pchisq(LR_uc[k], df = 1) # LR_uc obeys chi-square(1) distribution
  }
  
  table <- rbind(failuretimes, LR_uc, Pvalue) # record result
  return(table)
}


# (L1, L2) ~ M
VaREScalM <- function(grid, alpha) {
  # function for calculating VaR and ES when using comonotone copula(M)
  
  sumVaRM_sim <- matrix(0, nrow = 5000, ncol = length(alpha))
  sumESM_sim <- matrix(0, nrow = 5000, ncol = length(alpha))
  for(i in 1:5000) { # repeat 5000 times
    U <- runif(1500) # generate random uniform distribution numbers
    L1 <- quantile(grid[, 1], U) # use all samples to fit because there is no need for backtesting
    L2 <- quantile(grid[, 2], U) # get two loss samples with the comonotone copula(M) relation
    L <- L1 + L2 # get sum of losses sample
    sumVaRM_sim[i, ] <- quantile(L, alpha) # calculate VaR from sum of losses sample
    
    for(j in 1:length(alpha)) {
      sumESM_sim[i, j] <- mean(L[L >= sumVaRM_sim[i, j]]) # ES equals to mean of loss values which are greater than or equals to corresponding VaR
    }
  }
  
  sumVaRM <- apply(sumVaRM_sim, 2, mean) # use mean of simulated VaRs to approach the real VaR 
  sumESM <- apply(sumESM_sim, 2, mean) # use mean of simulated ESs to approach the real ES 
  table <- rbind(sumVaRM, sumESM) # record result
  return(table)

}

result <- function(stock, alpha) {
  # get well-prepared data and combine all results from functions above
  
  grid <- gridm(stock) 
  train <- seq(1, nrow(grid)-99) 
  losstest <- grid[-train, ] # use last 99 samples for backtesting copula model
  losstrain <- grid[train, ] # use the rest of samples for training to get copula model
  
  C <- VaRcalC(losstrain, alpha)
  bt <- Backtest(losstest, alpha, C)
  M <- VaREScalM(grid, alpha)
  table <- rbind(C, bt, M)
  return(table)
  
}


## Chosen pairs for research from all combinations in the portfolio

chosen <- rbind(c("FAN", "PZD"), c("GE", "PZD"), c("FSLR", "PZD"), c("RUN", "PZD"))
list1 <- result(chosen[1, ], alpha3)
list2 <- result(chosen[2, ], alpha3)
list3 <- result(chosen[3, ], alpha3)
list4 <- result(chosen[4, ], alpha3)


## Plot for VaR and ES results of all chosen pairs

par(mfrow = c(2, 2))
# Pair 1
plot(list1[1, ], list1[2, ], type = "l", xlab = expression(alpha),
     xlim = c(0.8, 1), ylab = expression("risk measure of sum loss"),
     ylim = c(0, max(list1[7, ])), col = "red", main = "Pair 1") 
lines(list1[1, ], list1[6, ], type = "l", col = "green")
lines(list1[1, ], list1[7, ], type = "l", col = "darkblue")
legend("topleft", bty = "n", y.intersp = 1.2, lty = rep(1, 3),
       col = c("red", "green", "darkblue"),
       legend = c(expression(VaR[alpha]~C),
                  expression(VaR[alpha]~M),
                  expression(ES[alpha]~M)))

# Pair 2
plot(list2[1, ], list2[2, ], type = "l", xlab = expression(alpha),
     xlim = c(0.8, 1), ylab = expression("risk measure of sum loss"),
     ylim = c(0, max(list2[7, ])), col = "red", main = "Pair 2")  
lines(list2[1, ], list2[6, ], type = "l", col = "green")
lines(list2[1, ], list2[7, ], type = "l", col = "darkblue")
legend("topleft", bty = "n", y.intersp = 1.2, lty = rep(1, 3),
       col = c("red", "green", "darkblue"),
       legend = c(expression(VaR[alpha]~C),
                  expression(VaR[alpha]~M),
                  expression(ES[alpha]~M)))

# Pair 3
plot(list3[1, ], list3[2, ], type = "l", xlab = expression(alpha),
     xlim = c(0.8, 1), ylab = expression("risk measure of sum loss"),
     ylim = c(0, max(list3[7, ])), col = "red", main = "Pair 3") #
lines(list3[1, ], list3[6, ], type = "l", col = "green")
lines(list3[1, ], list3[7, ], type = "l", col = "darkblue")
legend("topleft", bty = "n", y.intersp = 1.2, lty = rep(1, 3),
       col = c("red", "green", "darkblue"),
       legend = c(expression(VaR[alpha]~C),
                  expression(VaR[alpha]~M),
                  expression(ES[alpha]~M)))

# Pair 4
plot(list4[1, ], list4[2, ], type = "l", xlab = expression(alpha),
     xlim = c(0.8, 1),ylab = expression("risk measure of sum loss"),
     ylim = c(0, max(list4[7, ])), col = "red", main = "Pair 4") 
lines(list4[1, ], list4[6, ], type = "l", col = "green") 
lines(list4[1, ], list4[7, ], type = "l", col = "darkblue")
legend("topleft", bty = "n", y.intersp = 1.2, lty = rep(1, 3),
       col = c("red", "green", "darkblue"),
       legend = c(expression(VaR[alpha]~C),
                  expression(VaR[alpha]~M),
                  expression(ES[alpha]~M)))

# PCA analysis:
PCA <- prcomp(LogReturns)
summary(PCA)

# prcomp returns a list with class "prcomp" containing the following components:
# stdev, rotation: matrix of whose cols are eigenvectors, x: centered data, center: mean.
G   <- PCA$rotation # each column is the eigenvector g_i
mu  <- PCA$center
Y   <- PCA$x
var <- (PCA$sdev)^2

# PCA rotation changes the sign of loadings, but the interpretation is the same:
# (see e.g. https://stats.stackexchange.com/questions/88880/does-the-sign-of-scores-or-of-loadings-in-pca-or-fa-have-a-meaning-may-i-revers)
# revert signs to have positive loadings on PCA1 (to create index later on). 
FirstComp  <- -G[, 1]
SecondComp <- -G[, 2]

# Check that sum of variances of PCs equals sum of variances of original data
S <- var(LogReturns) 
totalVar <- sum(diag(S))
c(totalVar, sum(var))

# Plot first 6 PCAs relative to total variance:
FirstSix <- head(var/totalVar, n=-4) # pick first 6.

# Choose ylim that'll leave sufficient space above the tallest bar
ylim     <- c(0, 1.2*max(FirstSix))     

# Produce bar plot:
plt <-  barplot(FirstSix, xlab = "Principal Components", ylim = ylim,
                names=c("PCA1", "PCA2", "PCA3", "PCA4", "PCA5", "PCA6"),
                ylab = "Proportion of Variance", col ="#69b3a2",
                main = "Principal Component Analysis")

# Add value on top of bar:
text(x = plt, y = FirstSix, label = round(FirstSix, 3), pos = 3, cex = 0.8, col = "black")

par(mfrow = c(1, 2)) 
StockNames <- colnames(LogReturns)
barplot(FirstComp, names = StockNames, horiz = TRUE, cex.names = 0.6, 
        col = "#69b3a2", main = "PCA1", xlab = "Factor Loadings")
barplot(SecondComp, names = StockNames, horiz = TRUE, cex.names = 0.6, 
        col = "#69b3a2", main = "PCA2", xlab = "Factor Loadings")

# Calculate stock volatilities in 5-year period:
Volatilities <- 100*apply(LogReturns, 2, sd)

sort(Volatilities, decreasing = TRUE)             # 5yr period volatilities
sort(Volatilities, decreasing = TRUE)/sqrt(5)     # Annual volatilites
sort(Volatilities, decreasing = TRUE)/sqrt(5*252) # Daily volatilities
sort(FirstComp, decreasing = TRUE)                # PCA1 Loadings

# Study of dependence structure:
# ------------------------------
# 1. Contruct index: 
# Simple Pie Chart of index weights
par(mfrow = c(1, 1)) 
RelWeights     <- FirstComp/sum(FirstComp)
OrderedWeights <- sort(RelWeights, decreasing=TRUE)
lbls           <- c(names(OrderedWeights)[1:6], "OTHER")
slices         <- c(OrderedWeights[1:6], 1-sum(OrderedWeights[1:6]))
colors         <- c("light blue", "pink", "maroon", "light green", 
                    "light yellow", "orange", "gray") 
pie(slices, labels = lbls, main = "Pie Chart of Index weights", col = colors)

# Construct index by matmult LogRet with FirstComp
Index <- LogReturns %*% FirstComp # matrix multiplication
plot(Dates, Index, type = "l", ylab = "Log Returns", main = "PCA1 Index")
grid(10, 10, col = "lightgray", lty = "dotted", equilogs = FALSE)

par(mfrow = c(1, 1))
plot(Dates, LogReturns[, "FCEL"], type = "s", col = "dark blue",
     ylab = "Log Returns", main = "FCEL")
plot(Dates, Index, type="s", col = "red",
     ylab = "Log Returns", main = "PCA1 Index")

# Copula fitting: 
# Chose 3 stocks with high (FCEL), medium (ENPH) and low (FAN) weights.
# Copula fitting: (use same function in question 4)

# Append PCA1 Index to LogReturns matrix:
LogReturns           <- cbind(LogReturns, Index) 
colnames(LogReturns) <- c(StockNames, "PCA1 Index")

FitCopulas("PCA1 Index", "FCEL")
FitCopulas("PCA1 Index", "ENPH")
FitCopulas("PCA1 Index", "FAN")
