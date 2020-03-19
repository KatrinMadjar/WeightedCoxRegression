#' @title Stratified subsampling
#' 
#' @description
#' Stratified subsampling for real data application, stratified according to event indicator 
#' and subgroup indicator.

#' @param data [list] 
#' List with real data
#' @param ratio [numeric(1)]
#' Ratio for training set
#' 
#' @return [list]
#' List of length 2 with training and test data

StratifiedSubsampling <- function(data, ratio = .632, ...) {   
  study <- levels(data$group)  
  train <- list()
  test  <- list()
  for(i in 1:length(study)) {
      index1 <- which(data$group == study[i]) 
      status <- data$surv[index1, 2L]  
      complete <- complete.cases(data$surv[index1, ], data$clin[index1, ], data$genes[index1, ], data$group[index1])
      index2   <- index1[complete]
      train[[i]] <- unlist(lapply(split(index2, status[complete]),
                         function(id) sample(id, floor(length(id) * ratio))), use.names = FALSE )
      test[[i]]  <- setdiff(index2, train[[i]])
  }
  
  Traindata = data
  Traindata[c("surv","clin","genes")] = lapply(Traindata[c("surv","clin","genes")], function(x) x[unlist(train),] )
  Traindata$group = Traindata$group[unlist(train)]
  
  Testdata = data
  Testdata[c("surv","clin","genes")] = lapply(Testdata[c("surv","clin","genes")], function(x) x[unlist(test),] )
  Testdata$group = Testdata$group[unlist(test)]
  
  return(list("Traindata" = Traindata, "Testdata" = Testdata))
}
