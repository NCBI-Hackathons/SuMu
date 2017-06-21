#' Function for evaluation of model
#' @param data stanreg object with fit data
#' comparison of modeled vs observed data

#' Function for textual summary of data
#'
#' @param data stanreg object with fit data
#' @param features list object with specific features of interest (if not specified, all features are plotted)
#' @return summary dataframe with 1) genetic feature 2) posterior mean 3) posterior SD/SE 4) associated clinical/genetic variants
#' @export
#'
#'
feature_table <- function (data,features)
{
  #Order data
  if(missing(features))
    ordered_data=data_summary(data)
  else
    ordered_data=data_summary(data,features)

  #Report features correlated with genetic variants
  cov_table=data$covmat #Pull covariance matrix
  cor_table=cov2cor(cov_table)[-1,-1] #Generate correlation matrix

  n=nrow(cor_table)
  for (i in 1:n)
  {
    f=colnames(cor_table)[i] #select column to analyze
    col=abs(cor_table[,i]) #select values to check
    names=rownames(cor_table)[which(col>0.5)] #select names of rows with correlation coeff <-.5 or >.5
    ordered_data[f,4]=paste(paste(names,' '),collapse='') #convert to string
  }

  return(ordered_data)
}




#' Function for graphical summary of data
#' Plots mean +/- sd of genetic variants or clinical characteristic posteriors ordered by mean
#'
#' @param data stanreg object with fit data
#' @param features list object with specific genetic/clinical features of interest (if not specified, all features are plotted)
#' @export
#'
#'
#'
#'
feature_graph = function (data,features)
{
  #Order data
  if(missing(features))
    ordered_data=data_summary(data)
  else
    ordered_data=data_summary(data,features)

  #Prepare x/y variables
  avg=ordered_data$mean
  sd=ordered_data$sd
  variants=1:nrow(ordered_data)

  #Plot means & SDs
  plot(
    avg, variants,
    xlim=range(c(avg-sd, avg+sd)),
    yaxt="n",
    ylab="Gene Variants", xlab="Mean +/- SD",
    pch=19,
    main="Effects of Gene Variants on Survival"
  )
  arrows(avg-sd, variants, avg+sd, variants, length=0, angle=90, code=3) #Plot error bars
  axis(2,at=variants,labels=rownames(ordered_data),las=2) #Label axis with gene names
}

#' Helper function to summarize & order data
#'
#' @param data stanreg object with fit data
#' @param features list object with specific genetic/clinical features of interest (if not specified, all features are plotted)
#' @return mean_sd dataframe with mean/se/sd, ordered on mean
#'
data_summary = function(data,features)
{
  data_table=data$stan_summary
  data_table=data.frame(data_table[2:(nrow(data_table)-3),1:3])
  ordered_data=data_table[order(data_table$mean),]

  #If feature list is specified, subset for specified features
  if(!missing(features))
  {
    ordered_data=ordered_data[features,]
  }
  return(ordered_data)
}

#' Function for exploration of specific genetic feature
#'
#' @param mutations matrix with samples as rows, mutations as columns
#' @param clinical matrix of TCGA clinical data
#' @param survfit dataframe of predicted survival model
#' @param feature genetic variant (or clinical characteristic?) of interest
#'
#' @import survminer
#' @export
#'
#' Plot survival of population vs survival of subjects w/ genetic variant overlay w/ acutal data
#' Plot bar graph of covariance w/ genetic variant of clinical data
#'
view_feature = function(mutations, clinical, feature)
{
  #Plot observed data
  #Create survival table
  survival_table_cols=c("sampleID","OS","OS_IND")
  survival_table=clinical[survival_table_cols]

    #Identify samples with feature of interest & annotate samples w/o sequencing data
    mutated_samples=mutations$sample[mutations[feature]!=0]#sample names with feature
    not_sequenced=clinical$sampleID[match(clinical$sampleID,mutations$sample,nomatch=0)==0]

    for (i in 1:nrow(survival_table))
    {
      if(is.element(survival_table[i,1],mutated_samples))
      {  survival_table[i,4]="Mutated"
      } else if(is.element(survival_table[i,1],not_sequenced))
      {    survival_table[i,4]="Not_sequenced"
      } else
      {    survival_table[i,4]="No_Mutation_Detected"}
    }

    colnames(survival_table)[4]="Status"

  # Visualize with survival curve
  fit <- survfit(Surv(OS, OS_IND) ~ Status,
                 data = survival_table)
  survminer::ggsurvplot(fit, legend = "right", title = feature)

  #Plot simulated data
#  predicted_survival=posterior_survfit(survfit,)

#  ggplot(.,
#         aes(x = obstime, y = survpred, group = Status, colour = Status)) +
#    geom_line() +
#    geom_ribbon(aes(ymin = ci_lb, ymax = ci_ub, colour = NULL, fill = drug), alpha = 0.2) +
#    facet_wrap(~prevOI) +
#    theme_minimal() +
#    scale_y_continuous('Posterior-predicted survival', labels = percent) +
#    scale_x_continuous('Months')

  #Plot sorted correlations

}
