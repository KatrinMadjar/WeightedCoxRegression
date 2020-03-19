#' Simulate survival and gene expression data
#' Contains two functions for data simulation

# Weibull parameters for simulation of survival data
load("data/Weibull_param.RData")


#' Generate one data set with multivariate normally-distributed gene expression data
#' and Weibull-distributed survival data
#' 
#' @param n [integer(1)]
#' Number of observations
#' @param p [integer(1)]
#' Number of genes in total (prognostic genes + noise)
#' @param effect [vector]
#' Vector with true effects of first prognostic genes (remaining genes have no effect)
#' @param mw [vector]
#' Vector of length p with mean values of each gene
#' @param surv.e [list]
#' List of length 2 with 'scale' and 'shape' parameter for Weibull distribution of
#' event times
#' @param surv.c [list]
#' List of length 2 with 'scale' and 'shape' parameter for Weibull distribution of 
#' censoring times
#' 
#' @return [data.frame]
#' Data.frame of dimension nx(p+2) with columns 'time', 'status', and genes

simdat <- function(n, p, effect, mw, surv.e, surv.c ){  
  require(MASS)
  library(survival)
  
  k  <- length(effect)                       # number of prognostic genes
  bX <- c(effect, rep(0, p - k))            # vector of true effects for all genes
  
  sigma <- diag(1, p)             # assume uncorrelated genes
  X     <- mvrnorm(n, mw, sigma)  # generate gene expression data
  
  # simulate event times from Weibull distribution
  dt   <- (-log(runif(n)) * (1/surv.e$scale) * exp(-X %*% bX))^(1/surv.e$shape)
  
  # simulate censoring times from Weibull distribution
  cens <- rweibull(n, shape = surv.c$shape, scale = ((surv.c$scale)^(-1/surv.c$shape)))
  
  # observed time and status for each observation
  status <- ifelse(dt <= cens, 1, 0)
  time   <- pmin(dt, cens)
  
  data   <- data.frame( time, status, X)
  return(data)
}


#' Generate data for two differently distributed groups, each with 2 subgroups
#'
#' @param N [vector]
#' Vector of length 4 with number of observations per subgroup
#' @param P [integer(1)]
#' Number of genes per subgroup
#' @param effect1 [vector]
#' Vector with true effects of first prognostic genes (remaining genes have no effect) in group 1
#' @param effect2 [vector]
#' Vector with true effects of first prognostic genes (remaining genes have no effect) in group 2
#' @param mw1 [vector]
#' Vector of length P with mean values of each gene in group 1
#' @param mw2 [vector]
#' Vector of length P with mean values of each gene in group 2
#' @param surv.e1 [list]
#' List of length 2 with "scale" and "shape" parameter for Weibull distribution of
#' event times in group 1
#' @param surv.c1 [list]
#' List of length 2 with "scale" and "shape" parameter for Weibull distribution of 
#' censoring times in group 1
#' @param surv.e2 [list]
#' List of length 2 with "scale" and "shape" parameter for Weibull distribution of
#' event times in group 2
#' @param surv.c2 [list]
#' List of length 2 with "scale" and "shape" parameter for Weibull distribution of 
#' censoring times in group 2
#' 
#' @return [list]
#' List of length 4 with:
#' "surv": two-column matrix with columns named "time" and "status"
#' "genes": p-column matrix with gene expression data
#' "group": factor vector with subgroup levels
#' "cohorts": two-column matrix with columns named "cohort" und "sub_cohort"

sim.cohorts <- function(N, P, effect1, effect2, mw1, mw2, surv.e1, surv.c1, surv.e2, surv.c2)
{
  library(BBmisc)
  ngsc <- 2  # number of sub_cohorts per group 
  ng   <- 2  # number of groups
  
  cohorts <- expand.grid(sub_cohort = LETTERS[1:ngsc], cohort = 1:ng )
  cohorts.expanded <- cohorts[rep(seq_row(cohorts), N), 2:1]
  
  # generate data of first group
  group1 <- simdat(n = sum(N[1:ngsc]), p = P, 
                   effect = effect1, mw = mw1, 
                   surv.e = surv.e1, surv.c = surv.c1 )
  
  # generate data of second group
  group2 <- simdat(n = sum(N[(ngsc+1):(ngsc*ng)]), p = P, 
                   effect = effect2, mw = mw2, 
                   surv.e = surv.e2, surv.c = surv.c2 )
  
  data <- list('surv'  = cbind( time = c(group1$time, group2$time), 
                               status = c(group1$status, group2$status)), 
               'genes' = rbind(group1[,-(1:2)], group2[,-(1:2)]), 
               'group' = as.factor(paste0(cohorts.expanded$cohort, cohorts.expanded$sub_cohort)),
               'cohorts' = cohorts.expanded)
  return(data)
}

