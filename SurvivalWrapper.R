#' @title Building a Cox model based on (weighted) l1-regularized regression
#'
#' @description
#' First, training and test data are standardized.
#' Next, for weighted Cox models fixed weights or estimated weights are defined based 
#' on the training data.
#' Finally, the respective penalized Cox model is fitted and evaluated with the C-index.
#' 
#' @param job [Job]
#' Formal argument needed for algorithm function in 'batchtools' registry. 
#' Contains information about the job.
#' @param data [any]
#' Formal argument needed for algorithm function in 'batchtools' registry. 
#' @param instance [any]
#' Formal argument needed for algorithm function in 'batchtools' registry.
#' Passes the return value of the evaluated 'problem' function.
#' @param method.weights [character(1)]
#' Type of Cox model and type of weight, can be one of 
#' c("all", "sub", "lasso", "ridge", "rF", seq(0.1,0.9,0.1)).
#' "all": Standard combined model (unweighted).
#' "sub": Standard subgroup model (unweighted).
#' "lasso", "ridge", "rF": Proposed model with weights estimated by multinomial 
#' logistic regression with lasso or ridge penalty, or by random forest (rF).
#' "0.1",...,"0.9": Weighted model with fixed weights in (0,1).
#' @param covariates [character(1)]
#' Type of covariates, can be one of c("clin", "genes", "clin+genes")
#' "clin": only clinical variables (in real data, not simulation)
#' "genes": only genes (in real data and simulation)
#' "clin+genes": combination of clinical and genetic variables (only in real data), where only 
#' genes are penalized and clinical variables are included as mandatory (unpenalized)
#' @param gene.filter [character(1)]
#' Name of gene filter for real data application
#' 
#' @return [list]
#'  named list of length 5:
#'  "cox.Lambda": list of length S (S = number of subgroups); for each subgroup
#'  value of penalty parameter lambda.min (value of lambda giving minimum cvm).
#'  "cox.Beta": list of length S (S = number of subgroups); for each subgroup 
#'  named numeric vector with estimated regression coefficients.
#'  "cox.CI": list of length S (S = number of subgroups); C-index for each subgroup.
#'  "weights.fit": weights matrix with estimated or fixed weights.
#'  "subsampling": list of length 2 with results of object instance 
#'  (training and test data).

SurvivalWrapper <- function(job, data, instance, method.weights, 
                            covariates = "genes", gene.filter = NULL) {

  source("DataPreparation.R")
  source("EstimateWeights.R")
  source("CoxModelBuilder.R")
  
  library("mlr")  
  
  Traindata <- instance$Traindata  
  Testdata  <- instance$Testdata  

  # Combined model
  if(method.weights == "all"){
    weights.fit = NULL
    
    # Standardize training and test data 
    dprep = DataPreparation(train = Traindata, test = Testdata, mwei = method.weights, 
                            covars = covariates, gfilter = gene.filter)

    # Fit penalized (unweighted Cox model)
    all <- CoxModelBuilder(data = dprep$train, 
                           weights = rep(1L, nrow(dprep$train$covar)),  
                           pen.index = dprep$pen.index) 
    
    cox.CI <- cox.Beta <- cox.Lambda <- vector("list", nlevels(dprep$train$group))
    names(cox.CI) <- names(cox.Beta) <- names(cox.Lambda) <- levels(dprep$train$group)
    
    # Get Cox model parameters
    if(!is.null(all$modellfit)){
        beta <- as.numeric( predict(all$modellfit$learner.model$glmnet.fit,
                                    s = all$modellfit$learner.model$lambda.min,
                                    type = "coefficients") )
        names(beta) <- all$modellfit$features
        lambda <- all$modellfit$learner.model$lambda.min
     
      # Compute C-Index based on test data of each subgroup
      for (i in 1:nlevels(dprep$test$group)) {
        cox.Beta[[i]] <- beta
        cox.Lambda[[i]] <- lambda
        
        test.sub <- which(dprep$test$group == (levels(dprep$test$group)[i]))
        
        task.data           <- data.frame(dprep$test$surv, dprep$test$covar)
        colnames(task.data) <- make.names(colnames(task.data))
        tsk                 <- makeSurvTask("coxmodel", task.data, target = c("time", "status"))
        
        cox.CI[[i]] <- performance(predict(all$modellfit, task=tsk, subset=test.sub), measures = cindex )
      }
    }
  }
  else{
    # Standard subgroup model
    if(method.weights == "sub"){
      weights.fit <- NULL
      
      cox.CI <- cox.Beta <- cox.Lambda <- vector("list", nlevels(Traindata$group))
      names(cox.CI) <- names(cox.Beta) <- names(cox.Lambda) <- levels(Traindata$group)
      
      # For each subgroup i, fit Cox model based on training data of subgroup i
      # and compute C-index based on test data of subgroup i
      for (i in 1:nlevels(Traindata$group)) {
        sub.ids = which(Traindata$group == (levels(Traindata$group)[i])) 
        sub.ids. = which(Testdata$group == (levels(Testdata$group)[i])) 
        train.sub = Traindata; test.sub = Testdata;
        
        if(covariates != "genes"){
          train.sub[c("surv","clin","genes")] = lapply(train.sub[c("surv","clin","genes")], function(x) x[sub.ids,])
          test.sub[c("surv","clin","genes")] = lapply(test.sub[c("surv","clin","genes")], function(x) x[sub.ids.,])
        } else {
          train.sub[c("surv","genes")] = lapply(train.sub[c("surv","genes")], function(x) x[sub.ids,])
          test.sub[c("surv","genes")] = lapply(test.sub[c("surv","genes")], function(x) x[sub.ids.,])
        }
        
        dprep = DataPreparation(train = train.sub, test = test.sub, mwei = method.weights, 
                                covars = covariates, gfilter = gene.filter)
 
        sub <- CoxModelBuilder(data = dprep$train,
                               weights = rep(1L, nrow(dprep$train$covar)),  
                               pen.index = dprep$pen.index) 
        
        if(!is.null(sub$modellfit)){
            cox.Beta[[i]] <- as.numeric( predict(sub$modellfit$learner.model$glmnet.fit,
                                            s = sub$modellfit$learner.model$lambda.min,
                                            type = "coefficients") )
            names(cox.Beta[[i]]) <- sub$modellfit$features
            cox.Lambda[[i]] <- sub$modellfit$learner.model$lambda.min

            task.data           <- data.frame(dprep$test$surv, dprep$test$covar)
            colnames(task.data) <- make.names(colnames(task.data))
            tsk                 <- makeSurvTask("coxmodel", task.data, target = c("time", "status"))
          
            cox.CI[[i]] <- performance(predict(sub$modellfit, task=tsk), measures = cindex )
        }
    }
    }
    # Weighted model
    else{
      dprep = DataPreparation(train = Traindata, test = Testdata, mwei = method.weights,
                              covars = covariates, gfilter = gene.filter)
 
      # Estimate weights and generate weights matrix
      weights.fit <- EstimateWeights(train.set = dprep$train, 
                                     test.set = dprep$test,
                                     method  = method.weights)
      
      cox.CI <- cox.Beta <- cox.Lambda <- vector("list", nlevels(dprep$train$group))
      names(cox.CI) <- names(cox.Beta) <- names(cox.Lambda) <- levels(dprep$train$group)
      
      # For each subgroup i, fit weighted Cox model based on complete training data 
      # and compute C-index based on test data of subgroup i
      for (i in 1:nlevels(dprep$test$group)) {
        test.sub <- which(dprep$test$group == (levels(dprep$test$group)[i]))
        
        wei <- CoxModelBuilder(data = dprep$train, 
                               weights = weights.fit$weightmat[, i],  
                               pen.index = dprep$pen.index)     
        
        if(!is.null(wei$modellfit)){
            cox.Beta[[i]]        <- as.numeric( predict(wei$modellfit$learner.model$glmnet.fit,
                                                 s = wei$modellfit$learner.model$lambda.min,
                                                 type = "coefficients") )
            names(cox.Beta[[i]]) <- wei$modellfit$features
            cox.Lambda[[i]]      <- wei$modellfit$learner.model$lambda.min
          
            task.data           <- data.frame(dprep$test$surv, dprep$test$covar)
            colnames(task.data) <- make.names(colnames(task.data))
            tsk                 <- makeSurvTask("coxmodel", task.data, target = c("time", "status"))
            
            cox.CI[[i]] <- performance(predict(wei$modellfit, task=tsk, subset=test.sub), measures = cindex )
        }
      }
    }
  }
  
  print(sessionInfo())
  
  # Return results
  return(list(cox.Lambda = cox.Lambda,
              cox.Beta = cox.Beta,
              cox.CI = cox.CI,
              weights.fit = weights.fit,
              subsampling = instance 
             )) 
}