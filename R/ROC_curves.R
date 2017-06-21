library(dplyr)
library(survminer)
library(survival)
library(ggplot2)
library(scales)
library(rstan)
library(rstanarm)

#Test dataset
data(wells)
wells$dist100 <- wells$dist / 100
bfit <- stan_glm(
  switch ~ dist100 + arsenic, 
  data = wells,
  family = binomial(link = "logit"),
  prior_intercept = normal(0, 10),
  QR = TRUE,
  chains = 2,
  iter = 2000
)
print(bfit)

#non-dependent variables only. 
wells_n = wells[,-1]
var = wells[,1] #this is the observed outcome

#use the model in generative mode with its own data
p_data <- posterior_predict(bfit, newdata=wells_n)

p_est <- round(apply(p_data, 2, sum)/2000, 3) #calculate p-values
hist(p_est)#possible outcome of function

tDF <- tbl_df(data.frame(var, p_est))

br <- seq(from=0, to=1, length.out=1000)
TPR <- c()
FPR <- c()
for (i in br){
  tDF1 <- filter(tDF, p_est<=i)
  tDF0 <- filter(tDF, p_est>i)
  
  size1 <- dim(tDF1)[1]
  size0 <- dim(tDF0)[1]
    
  true_pos =  sum(tDF1$var)
  false_pos = size1 - true_pos
  true_neg = size0 - sum(tDF0$var)
  false_neg = sum(tDF0$var)
  
  TPR <- c(TPR, true_pos/(true_pos+false_neg))
  FPR <- c(FPR, 1-(true_neg/(true_neg+false_pos)))

}

plot(TPR,FPR,type="s")
abline(a=0,b=1, col="yellow4")


library(zoo)

x <- 1:10
y <- 3*x+25
id <- order(x)

AUC <- sum(diff(x[id])*rollmean(y[id],2))

