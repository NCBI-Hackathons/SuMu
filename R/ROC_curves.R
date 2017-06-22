# AUC from posterior distributions of stan models


#' Function that uses the raw data and a stanreg object of the model fitted with them.
#'
#' It first uses the model in a generative fashion to generate posterior distributions for each observation.
#' It then utilizes the generated distributions of the binary outcome to compute the probability that the draw will be "1" (eg. alive)
#' 
#' @param data_frame the dataframe containing clinical data
#' @param response_name the column name of the data_frame that has the resonse variable
#' @param fitted_model the stanreg fitted model
#' @param h_gram do you want a histogram of all the estimated posterior p-values? Possible values TRUE/FALSE, defaults to FALSE
#' @param roc_plot do you want a ROC curve produced at the end? Possible values TRUE/FALSE, defaults to FALSE
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
auc <- function(data_frame, response_name, fitted_model, h_gram=FALSE, roc_plot=FALSE){
  var = data_frame[response_name][,1]  #observed outcomes
  ind_vars = data_frame[ , ! colnames(data_frame) %in% c(response_name) ] #only independent variables
  p_data <- posterior_predict(fitted_model, newdata=data_frame)
  p_est <- round(apply(p_data, 2, sum)/dim(p_data)[1], 3) #calculate p-values
  if (h_gram==TRUE) {hist(p_est, xlab= "Estimated posterior p-values")}

  tDF <- tbl_df(data.frame(var=var, p_est=p_est))
  br <- seq(from=min(p_est), to=max(p_est), length.out=1000)
  TPR <- c()
  FPR <- c()

  for (i in br){
    tDF1 <- filter(tDF, p_est<=i)
    tDF0 <- filter(tDF, p_est>i)

    size1 <- dim(tDF1)[1]
    size0 <- dim(tDF0)[1]

    if (size1==0 | size0==0){next}
    true_pos =  sum(tDF1[,1])
    false_pos = size1 - true_pos
    false_neg = sum(tDF0[,1])
    true_neg = size0 - false_neg


    TPR <- c(TPR, true_pos/(true_pos+false_neg))
    FPR <- c(FPR, 1-(true_neg/(true_neg+false_pos)))
  }
  id <- order(TPR)
  AUC <- sum(diff(TPR[id])*rollmean(FPR[id],2))

  if (roc_plot==TRUE) {
    plot(TPR,FPR, type="s")
    abline(a=0, b=1, col="yellow4")
    mtext(paste("AUC : ",round(AUC,2)))
  }

  return(AUC)
}
