#' @title Fit Cox model
#'
#' @description
#' Function to fit (weighted) l1-penalized Cox model
#' 
#' @param data [list]
#' list with training data
#' "covar": design matrix
#' "surv": two-column matrix with columns named 'time' and 'status'
#' @param weights [vector]
#' vector with individual weights for each training sample
#' @param pen.index [vector]
#' penalty value for each covariate
#' 
#' @return [list]
#' list of length 1 with Cox model fit (object of class WrappedModel)

CoxModelBuilder <- function(data, weights, pen.index){  
  require(mlr)

  if ( sum(weights == 1L) == nrow(data$surv) ) weights <- NULL
  
  task.data           <- data.frame(data$surv, data$covar)
  colnames(task.data) <- make.names(colnames(task.data))
  tsk                 <- makeSurvTask("coxmodel", task.data, target = c("time", "status"))

  lrn   <- makeLearner("surv.cvglmnet", alpha = 1L, s = "lambda.min", 
                     penalty.factor = pen.index, standardize = FALSE)  

  m.fit <- tryCatch( train(lrn, tsk, weights=weights),  
                    error = function(e) NULL )
  
  return(list(modellfit = m.fit))  
}
