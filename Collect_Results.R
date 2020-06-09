# Collect, combine and save results  

library(batchtools)

# Collect simulation results or data application?
sim = TRUE # FALSE 

####
if(sim) {
  reg = loadRegistry(file.dir = "Simulation")
  
  prob.design  = expand.grid(n = c(seq(20,100,10),200,500,1000),
                             p = c(12,100,200),  
                             epsilon = c(seq(0,0.5,0.1),1) )
  surv.design  = data.frame(model = c("all", "sub", "lasso", "ridge", "rF", seq(0.1,0.9,0.1)))
  args         = CJ(k = 1:nrow(prob.design), l = 1:nrow(surv.design)) 
  pars         = cbind.data.frame(prob.design[args$k,], model = surv.design[args$l,], stringsAsFactors = F)
  
  auc.train = acc.train = auc.test = acc.test = matrix(ncol=100, nrow=nrow(pars))

  nf = "Results_Sim" # new folder name to save results
  
####
}else{
  surv.design = CJ(
    method.weights = c("all", "sub", "lasso", "ridge", "rF", seq(0.1,0.9,0.1)),
    gene.filter    = c("all", "Top1000Varianz", "prognostic")
  )
  pars = surv.design
  
  nf = "Results_App" # new folder name to save results
}

####

weights = cox.betas = cox.CI = vector("list", length = nrow(pars))

for(j in 1:nrow(pars)){
  
  if(sim){
    ids = findExperiments(algo.pars = (method.weights == pars[j,"model"]), reg = reg)
    ids = findExperiments(ids = ids, prob.pars = (p == pars[j,"p"] & n == pars[j,"n"] & epsilon == pars[j,"epsilon"]),
                          reg = reg)
  }else{
    ids = findExperiments(algo.pars = (method.weights == pars$method.weights[j] & 
                                         gene.filter == pars$gene.filter[j] ),
                          reg = reg)
  }
  
  # 100 subsampling IDs for a specific parameter configuration
  xx = reduceResultsList(ids = ids, reg = reg)
  
  ####
  # Weights estimation:
  if(mw %in% c("rF", "lasso", "ridge") ){
    
    # For the simulations, we can calculate AUC and accuracy based on the test and cross-validated training data
    if(sim){
      # CV performance on training sets:
      pred = lapply(xx, function(y) y$weights.fit$pred.train)
      
      truth = lapply(pred, function(y){tr = y$truth; tr = as.factor(substr(tr, 1, 1)); tr } ) 
      pprob = lapply(pred, function(y){pr = y[,grep("prob.", colnames(y))]; pr = rowSums(pr[,3:4]); pr } ) 
      resp = lapply(pred, function(y){re = y$response; re = as.factor(substr(re, 1, 1)); re } )  
      
      auc.train[j,] = sapply(1:length(truth), function(i) measureAUC(probabilities = pprob[[i]], truth = truth[[i]], positive = "2") )
      acc.train[j,] = sapply(1:length(truth), function(i) measureACC(truth = truth[[i]], response = factor(resp[[i]],levels=c("1","2")) ) ) 
      
      # Performance on test sets
      pred = lapply(xx, function(y) y$weights.fit$pred.test)
      
      truth = lapply(pred, function(y){tr = y$truth; tr = as.factor(substr(tr, 1, 1)); tr } )  
      pprob = lapply(pred, function(y){pr = y[,grep("prob.", colnames(y))]; pr = rowSums(pr[,3:4]); pr } )  
      resp = lapply(pred, function(y){re = y$response; re = as.factor(substr(re, 1, 1)); re } )  
      
      auc.test[j,] = sapply(1:length(truth), function(i) measureAUC(probabilities = pprob[[i]], truth = truth[[i]], positive = "2") )
      acc.test[j,] = sapply(1:length(truth), function(i) measureACC(truth = truth[[i]], response = factor(resp[[i]],levels=c("1","2")) ) ) 
    }
    
    # Estimated weights
    wei = lapply(1:length(xx), function(i){ 
      mat = xx[[i]]$weights.fit$weightmat; 
      tr = xx[[i]]$weights.fit$pred.train$truth[order(xx[[i]]$weights.fit$pred.train$id,decreasing=F)]; 
      cbind.data.frame(mat, truth = tr, iter = i) 
      })
    wei = do.call("rbind",wei)
    weights[[j]] = wei
    
  }else{  # no weights estimation
    if(sim){
      auc.train[j,] = acc.train[j,] = auc.test[j,] = acc.test[j,] = rep(NA, 100)
    }
  }
  
  ####
  # Cox model results:
  # Betas 
  bet = lapply(xx, function(y) do.call("rbind", y$cox.Beta) )
  bet2 = lapply(1:length(bet), function(i){yy = data.frame(iter = rep(i,nrow(bet[[i]])), subgroup = rownames(bet[[i]]), bet[[i]]); yy})
  bet2 = do.call("rbind", bet2)
  cox.betas[[j]] = bet2
  
  # C index
  ci = lapply(xx, function(y) do.call("rbind", y$cox.CI) )
  ci2 = lapply(1:length(ci), function(i){yy = data.frame(iter = rep(i,nrow(ci[[i]])), subgroup = rownames(ci[[i]]), ci[[i]]); yy})
  ci2 = do.call("rbind", ci2)
  cox.CI[[j]] = ci2
  
  print(j)
}

dir.create(nf)

if(sim){
  save(pars, auc.train, acc.train, auc.test, acc.test, file = paste0(nf, "/AUC_ACC.RData"))
}
save(pars, weights, file = paste0(nf, "/Weights.RData"))
save(pars, cox.betas, file = paste0(nf, "/Betas_CoxM.RData"))
save(pars, cox.CI, file = paste0(nf, "/Cindex.RData"))
