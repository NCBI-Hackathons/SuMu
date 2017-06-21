ROC_ <- function(data_frame, response_name, fitted_model){
  
 # '''This function takes as inputs the dataframe with the raw data, 
 # the column name of the response variable, and the fitted model (trained on a different),
 # calculates AUC for the model and spits out some pretty plots too (no plots yet)'''
  
  var = data_frame[response_name][,1]  #observed outcomes
  ind_vars = data_frame[ , ! colnames(data_frame) %in% c(response_name) ] #only independent variables
  p_data <- posterior_predict(fitted_model, newdata=ind_vars)
  p_est <- round(apply(p_data, 2, sum)/dim(p_data)[1], 3) #calculate p-values
  
  tDF <- tbl_df(data.frame(var=var, p_est=p_est))
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
  
  AUC <- sum(diff(TPR[id])*rollmean(FPR[id],2))
  return(AUC)
}
