#' Standardize numerical covariates using parameters of the training data
#'
#' @param train [list]
#' List with training data:
#' "surv": two-column matrix with columns named "time" and "status"
#' "genes": p-column matrix with gene expression data
#' "group": factor vector with subgroup levels
#' "clin": data.frame with clinical covariates, only for real data, not in simulations
#' for real data there are also list elements with gene indices for gene filters
#' @param test [list]
#' List with test data analogously to training data
#' @param mwei [character(1)]
#' Type of Cox model and type of weight, can be one of 
#' c("all", "sub", "lasso", "ridge", "rF", seq(0.1,0.9,0.1)).
#' "all": Standard combined model (unweighted).
#' "sub": Standard subgroup model (unweighted).
#' "lasso", "ridge", "rF": Proposed model with weights estimated by multinomial 
#' logistic regression with lasso or ridge penalty, or by random forest (rF).
#' "0.1",...,"0.9": Weighted model with fixed weights in (0,1).
#' @param covars [character(1)]
#' Type of covariates, can be one of c("clin", "genes", "clin+genes")
#' "clin": only clinical variables (in real data, not simulation)
#' "genes": only genes (in real data and simulation)
#' "clin+genes": combination of clinical and genetic variables (only in real data), where only 
#' genes are penalized and clinical variables are included as mandatory (unpenalized)
#' @param gfilter [character(1)]
#' Name of gene filter for real data application
#'  
#' @return [list]
#' List of length 3 with
#' "train": standardized training data
#' "test": standardized test data
#' "pen.index": penalty value for each covariate

DataPreparation <- function(train, test, mwei, covars, gfilter){
  library("dummy")
  
  # apply gene filter for real data
  if(!is.null(gfilter)){ 
    train$genes = as.matrix(train$genes[,train[[gfilter]]])
    test$genes = as.matrix(test$genes[,test[[gfilter]]])
  }else{ # no gene filter for simulated data
    train$genes = as.matrix(train$genes) 
    test$genes = as.matrix(test$genes) 
  }
  
  # standardize genes with parameters of training data 
  if(covars %in% c("genes", "clin+genes") ){
    
    mu = apply(train$genes, 2, mean)
    sd = apply(train$genes, 2, sd)
    
    train$genes = scale(train$genes, scale=TRUE)
    test$genes = scale(as.matrix(test$genes), center = mu, scale = sd)
  }
  
  # standardize survival time for weights estimation, since time is included as covariate in 
  # classification
  if(mwei %in% c("lasso", "ridge", "rF") ){ 
    time = as.vector( scale(train$surv[,1L], center = mean(train$surv[,1L]),  scale = sd(train$surv[,1L])) )
    time. = as.vector( scale(test$surv[,1L], center = mean(train$surv[,1L]),  scale = sd(train$surv[,1L])) )
    train$surv2 = train$surv; train$surv2[, 1L] = time
    test$surv2 = test$surv; test$surv2[, 1L] = time.
  }
  
  # standardize numerical clinical variables (only in real data, not simulation)
  if(covars %in% c("clin", "clin+genes") ){  
    fac.id = sapply(1:ncol(train$clin), function(i) is.factor(train$clin[,i]) )
    covar2 = as.data.frame(train$clin[,!fac.id])
    covar2. = as.data.frame(test$clin[,!fac.id])
    if(ncol(covar2)==0){
      covar = as.matrix( dummy( train$clin, int=TRUE ) )
      covar. = as.matrix( dummy( test$clin, int=TRUE ) )
    }else{
      mu = apply(as.data.frame(covar2), 2, mean)
      sd = apply(as.data.frame(covar2), 2, sd)
      
      covar2 = scale(covar2, center = mu, scale = sd)
      colnames(covar2) = colnames(train$clin)[!fac.id]
      
      covar2. = scale(covar2., center = mu, scale = sd)
      colnames(covar2.) = colnames(test$clin)[!fac.id]
      
      covar = cbind(covar2, as.matrix( dummy( train$clin, int=TRUE ) ))
      covar. = cbind(covar2., as.matrix( dummy( test$clin, int=TRUE ) ))
    }
  }
  
  if(covars == "clin") {
    train$covar = covar 
    test$covar = covar. 
    pen.index = rep(1L, ncol(train$covar))
  }
  if(covars == "genes") {
    train$covar = train$genes
    test$covar = test$genes
    pen.index = rep(1L, ncol(train$covar))
  }
  if(covars == "clin+genes") {
    n.unpen = ncol(covar)
    # penalize only genes, not clinical variables:
    pen.index = c(rep(0L, n.unpen), rep(1L, ncol(train$genes)))    
    train$covar = cbind(covar, train$genes)
    test$covar = cbind(covar., test$genes)
  }
  
  if(mwei == "all"){
    # for combined Cox model, subgroup indicator variable is included as covariate:
    group = as.matrix( dummy( data.frame(group = train$group), int = TRUE ) )
    train$covar = cbind(train$covar, group) 
    test$covar = cbind(test$covar, 
                       as.matrix( dummy( data.frame(group = test$group), int = TRUE ) ))   
    
    # change penalty 
    if(covars == "clin"){
      pen.index = rep(1L, ncol(train$covar))  
    }
    else{
      pen.index = c(pen.index, rep(0L, ncol(group)))
    }
  }
  
  return(list("train" = train, "test" = test, "pen.index" = pen.index))
}
