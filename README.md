# Weighted Cox Regression

R source code and data associated with the publication *Madjar K and Rahnenführer J (2020): Weighted Cox regression for the prediction of heterogeneous patient subgroups. arXiv:2003.08965*.

Our method is focused on the situation of predefined, possibly heterogenous subgroups of patients with available survival endpoint and high-dimensional molecular measurements such as gene expression
data, with the aim of obtaining a separate risk prediction model for each subgroup.
For this purpose, we propose a l1-penalized Cox regression model with a weighted version of the partial likelihood that includes patients of all subgroups but assigns them individual weights.
The weights are estimated from the data such that patients who are likely to belong to the subgroup of interest receive a higher weight.
We compare our approach to fixed weights and unweighted Cox models in simulations and application to real lung cancer cohorts.


The main file to run the simulation study is **Run_Simulation.R** and the main file for the real data application is **Run_RealDataApplication.R**.


## Overview of R files:

#### Run_Simulation.R

Setup and running of the simulation, including computation and evaluation of our proposed weighted Cox model, as well as the two standard unweighted Cox models and the weighted Cox model with fixed weights.
All simulations are run using the R package batchtools for parallelization.

#### Run_RealDataApplication.R

Setup and running of the real data application, including computation and evaluation of our proposed weighted Cox model, as well as the two standard unweighted Cox models and the weighted Cox model with fixed weights.
All models can be run with different gene filters and covariate sets (here only genetic covariates are included).
All settings are run using the R package batchtools for parallelization.

#### DataSimulation.R

Helper functions for generation of simulated (training and test) data. Only needed in simulation study, not real data application.

#### StratifiedSubsampling.R

Function for stratified random subsampling and generation of training and test data ("problem" function in batchtools, see file **Run_RealDataApplication.R**). Only needed in real data application, not in simulation study.

#### SurvivalWrapper.R

Wrapper function for running of a specific type of Cox model with, in the case of the real data application, a specific covariate set and gene filter.
The cox model is fitted based on the training data and evaluated based on the test data.
This function includes the following three functions.

#### DataPreparation.R

Function for preparation of the training and test data for model fitting and evaluation: standardization of numerical covariates, application of the gene filter to the genetic covariates (if defined), and definition of a penalty value for each covariate (only interesting in the case of a combination of clinical and genetic covariates, where only genes are penalized and clinical variables are included as mandatory (unpenalized) into the model).

#### EstimateWeights.R

Function to determine individual weights for each patient in the training set and in each subgroup model (either estimated weights as proposed by us with different classification methods or fixed weights as proposed by Weyer and Binder, 2015).

#### CoxModelBuilder.R

Function to fit (weighted) l1-penalized Cox model.

## Data:

#### Weibull_param.RData

RData object with parameters of the Weibull distribution used for data simulation (see file **DataSimulation.R**)

#### Realdata_LC.RData

RData object with preprocessed data of lung cancer cohorts, including information on survival endpoint, clinical covariates, genetic covariates, cohort membership of each patient, and gene filters (see file **Run_RealDataApplication.R**)
