# Run the simulation with all settings described in the paper

library("batchtools")

# Batchtools registry to save all information regarding computational jobs and all results
reg <- makeExperimentRegistry(file.dir = "Simulation" )  

# Parameter combinations for data simulation
VP <- expand.grid(n = c(seq(20,100,10),200,500,1000),
                  p = c(12,100,200),  
                 epsilon = c(seq(0,0.5,0.1),1) )

# Set true effects of first prognostic genes in both groups:
beta1 = c(1,1,0,0, -0.5, 0.5, 0.75,0.25, -1,-1,-0.75,-0.25)
beta2 = c(0,0,1,1,  0.5,-0.5, 0.25,0.75, -1,-1,-0.75,-0.25)

# "Problem" function in batchtools that creates problem instance (here: generation of training
# and test data) used in the "algorithmic" function below. 
prob.fun = function( n, p, epsilon, ...){
  source("DataSimulation.R")  
  
  mu1 = 4+2*epsilon # mean value for medium effect of absolute size 0.5 or 0.75
  mu2 = 4+4*epsilon # mean value for strong effect of absolute size 1
  
  Traindata = sim.cohorts(N = rep(n,4), P = p,  
                           effect1 = beta1, effect2 = beta2,
                           # Set mean values of all genes according to true effects:
                           mw1 = c( c(mu2,mu2,4,4, mu1,mu1, mu1,4,mu2,mu2,mu1,4), rep(4,p-length(beta1)) ), 
                           mw2 = c( c(4,4,mu2,mu2, mu1,mu1, 4,mu1,mu2,mu2,mu1,4), rep(4,p-length(beta1)) ),   
                           # Weibull parameters for simulation of event and censoring times in both groups
                           # (we use the same parameters for event and censoring times to obtain about 
                           # 50% censoring rates)
                           surv.e1 = Surv.e[[1]], surv.c1 = Surv.e[[1]], 
                           surv.e2 = Surv.e[[2]], surv.c2 = Surv.e[[2]])
  
  # Generate test data of the same size and with the same parameters as the training data
  Testdata = sim.cohorts(N = rep(n,4), P = p,  
                          effect1 = beta1,
                          effect2 = beta2,
                          mw1 = c( c(mu2,mu2,4,4, mu1,mu1, mu1,4,mu2,mu2,mu1,4), rep(4,p-length(beta1)) ), 
                          mw2 = c( c(4,4,mu2,mu2, mu1,mu1, 4,mu1,mu2,mu2,mu1,4), rep(4,p-length(beta1)) ),   
                          surv.e1 = Surv.e[[1]], surv.c1 = Surv.e[[1]], 
                          surv.e2 = Surv.e[[2]], surv.c2 = Surv.e[[2]])

  return(list("Traindata" = Traindata, "Testdata" = Testdata))
}

prob.designs = vector("list", 1)
prob.designs[[1]] = VP
names(prob.designs)[1] = "Simdata"

addProblem(name  = "Simdata",
           data  = NULL,
           fun   = prob.fun,
           seed  = 20181128
)

source("SurvivalWrapper.R") # Algorithm function 
addAlgorithm(name  = "survival", fun = SurvivalWrapper)

# Parameters for algorithm function: different Cox model types
surv.design = CJ(
  method.weights = c("all", "sub", "lasso", "ridge", "rF", seq(0.1,0.9,0.1))
)

algo.designs = list(survival = surv.design)

# Add experiments to batchtools registry (combinations of problems and algorithms) to define jobs
addExperiments(prob.designs  = prob.designs,
               algo.designs  = algo.designs,
               repls = 100 # Number of subsampling (training and test) sets
)

submitJobs(ids = findNotSubmitted())  # Run (all) jobs or run single job with job ID 1 (example):
submitJobs(1)