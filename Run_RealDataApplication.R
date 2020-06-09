library("batchtools")

# Load preprocessed data of lung cancer cohorts, including survival data (list element "surv"),
# clinical covariates ("clin"), genetic covariates ("genes") and cohort membership of each 
# patient ("group"). The last three list elements refer to indices of genetic covariates (columns
# in "genes") used as gene filters.
load("data/Realdata_LC.RData")  

# Batchtools registry to save all information regarding computational jobs and all results
reg = makeExperimentRegistry(file.dir = "Lungcancer", 
                             seed     = 24420 # for Biometrical Journal
#                            seed     = 251157152 # for arXiv Submission
                             )  

# "Problem" function in batchtools that creates problem instance (here: stratified subsampling 
# and generation of training and test data) used in the "algorithmic" function below. 
source("StratifiedSubsampling.R")

prob.designs = vector("list", 1)
prob.designs[[1]] = data.frame(ratio = .632) 
names(prob.designs)[1] = "subsampling"

addProblem(name = "subsampling", 
           data = Data, 
           fun = StratifiedSubsampling, 
           seed  = 244 # for Biometrical Journal
#          seed  = 16032018 # for arXiv Submission
           )

source("SurvivalWrapper.R") # Algorithm function 
addAlgorithm(name  = "survival", fun = SurvivalWrapper)

# Parameters for algorithm function (different Cox model types and gene filters, different sets
# of covariates are also possible, but here only "genes" are used)
surv.design = CJ(
  method.weights = c("all", "sub", "lasso", "ridge", "rF", seq(0.1,0.9,0.1)),
  gene.filter    = c("all", "Top1000Varianz", "prognostic"),
  covariates     = "genes" # can be one of c("genes", "clin", "clin+genes")
)

algo.designs = list(survival = surv.design)

# Add experiments to batchtools registry (combinations of problems and algorithms) to define jobs
addExperiments(prob.designs  = prob.designs,
               algo.designs  = algo.designs,
               repls = 100L  # Number of subsampling (training and test) sets
)

submitJobs(ids = findNotSubmitted())  # Run (all) jobs or run single job with job ID 1 (example):
submitJobs(1)
