#' @title Estimate subgroup weights
#' 
#' @description
#' Function to determine individual weights for each patient in the training set and in each 
#' subgroup model (either estimated weights as proposed by us with different classification 
#' methods or fixed weights as proposed by Weyer and Binder, 2015).
#' 
#' @param train.set [list]
#' list with training data
#' "covar": design matrix
#' "surv": two-column matrix with columns named 'time' and 'status'
#' "surv2": two-column matrix with columns named 'time' and 'status', with 'time'
#' standardized to be included as covariate in weights estimation.
#' "group": factor vector with subgroup levels
#' @param test.set [list]
#' list with test data, analogous to train.set
#' @param method [character(1)]
#'  can either be one of c("lasso", "ridge", "rF"): proposed model with weights 
#'  estimated by multinomial logistic regression with lasso or ridge penalty or by 
#'  random forest ("rF");
#'  or any numeric value in the interval (0,1)
#'  
#' @return [list]
#' named list containing the following objects:
#' "weightmat": matrix with weights for each patient (rows) and each subgroup (column)
#' if method is one of c("lasso", "ridge", "rF") additional list elements are:
#' "pred.train": prediction based on cross-validated training data  
#' "var.imp": for random forest mean decrease in gini index; for multinomial logistic
#' regression: estimated regression coefficients and penalty parameter lambda.min 
#' (value of lambda giving minimum cvm).
#' "pred.test": prediction based on test data

EstimateWeights <- function(train.set, test.set, method){
  
  if(method %in% c("lasso", "ridge", "rF")){   # estimate weights
    library(mlr)

    # combine task data with covariates and categorical outcome y (= subgroup variable)
    task.data      <- data.frame(train.set$surv2, train.set$covar, y = train.set$group )
    task.data.test <- data.frame(test.set$surv2, test.set$covar, y = test.set$group )
  
    tsk <- makeClassifTask(id = "weights", data = task.data, target = "y")
    
    # 10-fold CV stratified according to outcome y 
    rsmpl <- makeResampleDesc(method = "CV", iters = 10, stratify = TRUE) 
    
    if(method == "rF"){  # random forest
      lrn <- makeLearner("classif.ranger", predict.type = "prob", par.vals = list(importance="impurity"))  
      
      res <- resample(learner = lrn, task = tsk, resampling = rsmpl, 
                     extract = function(x) x$learner.model$variable.importance)
      
      # variable importance (mean decrease in gini index):
      var.imp <- res$extract
    } 
    else{ # penalized multinomial logistic regression
      if (method == "lasso") alpha = 1L
      if (method == "ridge") alpha = 0L   
      
        lrn <- makeLearner("classif.cvglmnet", predict.type = "prob", alpha = alpha, 
                          s = "lambda.min", standardize = FALSE)

        res <- resample(learner = lrn, task = tsk, resampling = rsmpl, 
                       measures = list(multiclass.au1p, acc, mmce), models = TRUE)

        # For each CV step, one value of lamda.min and one named vector of estimated
        # regression coefficients is obtained.
        lambda  <- lapply(res$models, function(x) x$learner.model$lambda.min )
        betas   <- lapply(res$models, function(x){s = which(x$learner.model$glmnet.fit$lambda == x$learner.model$lambda.min);
                                               bet = lapply(x$learner.model$glmnet.fit$beta, function(xx) xx[,s]);
                                               return(bet) })                      
        var.imp <- list(beta = betas, lambda.min = lambda)
    }

    # resample function changes order of training indices.
    # Thus, re-ordering of training indices:
    pred.prob <- res$pred$data[order(res$pred$data$id,decreasing=F), grep("prob.", colnames(res$pred$data))] 
    
    # compute likelihood ratios (weights)
    frequencies <- table(task.data$y) / length(task.data$y)   
    weights.mat <- t(t(pred.prob) / as.numeric(frequencies)) 
    
    # classifier fit based on training data
    m.fit <- train(lrn, tsk)
    
    # prediction based on test data
    tsk.test <- makeClassifTask(id = "weights", data = task.data.test, target = "y")
    pred.test <- predict(m.fit, tsk.test)
    
    return(list(pred.train = res$pred$data, 
                var.imp    = var.imp, 
                weightmat  = weights.mat,
                pred.test  = pred.test$data 
                ))    
  }
  else{  # fixed weights in (0,1)
    nu <- as.numeric(method)
    weights.mat <- sapply(levels(train.set$group),
                            function(g) ifelse(train.set$group == g, 1, nu))
    
    # Observations in the training set belonging to subgroup g obtain weight 1 in
    # column g and otherwise weight nu.
    return(list(weightmat = weights.mat))
  }
}
