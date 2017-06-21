# Functions for for calculating AUC from posterior distributions of stan models


#' Function that uses the raw data and a stanreg object of the model fitted with them.
#'
#' It first uses the model in a generative fashion to generate posterior distributions for each observation.
#' It uses the
#' @param data_frame the dataframe containing clinical data
#' @param response_name the column name of the data_frame that has the resonse variable
#' @param fitted_model the stanreg fitted model
#' @param h_gram do you want a histogram of all the estimated posterior p-values? Possible values TRUE/FALSE, defaults to FALSE
#'
#' @return numeric AUC
#' @import survminer
#' @import survival
#' @import dplyr
#' @import ggplot2
#' @import scales
#' @import rstan
#' @import rstanarm
#' @import zoo
#' @export
auc <- function(data_frame, response_name, fitted_model, h_gram=FALSE){
  var = data_frame[response_name][,1]  #observed outcomes
  ind_vars = data_frame[ , ! colnames(data_frame) %in% c(response_name) ] #only independent variables
  p_data <- posterior_predict(fitted_model, newdata=data_frame)
  p_est <- round(apply(p_data, 2, sum)/dim(p_data)[1], 3) #calculate p-values
  if (histogram==TRUE) {hist(p_est, xlab= "Estimated posterior p-values")}

  tDF <- tbl_df(data.frame(var=var, p_est=p_est))
  br <- seq(from=min(p_est), to=max(p_est), length.out=1000)
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

  AUC <- sum(diff(TPR[id])*rollmean(FPR[id],2))
  return(AUC)
}
