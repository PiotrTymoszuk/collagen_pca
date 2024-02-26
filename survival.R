# Analysis of prognostic relevance of the collagen-related gene expression
#
# 1) 'Univariable' analyses of prognostic relevance of single genes for OS 
# and RFS. The variables are stratified by their respective optimal cutoffs 
# corresponding to the largest difference in survival as assessed by 
# Mentel-Henszel test. Significant effects are defined as significant genes 
# shared by at least four cohorts

# tools -------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)
  
  library(clustTools)

  library(survival)
  library(survminer)
  library(glmnet)
  library(coxExtensions)
  library(kmOptimizer)

  library(survivalsvm)
  library(randomForestSRC)
  
  library(furrr)
  library(soucer)
  library(ggtext)

  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  set_rownames <- trafo::set_rownames
  
  c('./tools/globals.R', 
    './tools/functions.R', 
    './tools/svm_tools.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# analysis data and globals --------
  
  insert_msg('Analysis data and globals')
  
  source_all('./survival scripts/globals.R', 
             message = TRUE, crash = TRUE)
  
# Univariable analysis ---------
  
  insert_msg('Univariable analysis')
  
  ## resorting to cached versions, because the analyses take a while
  
  access_cache(cache_path = './cache/rfs_cut.RData', 
               script_path = './survival scripts/univariable_rfs.R', 
               message = 'Loading cached univariable survival analysis results')

# Multi-parameter survival modeling ---------
  
  insert_msg('Multi-parameter modeling')
  
  list(cache_path = c('./cache/elnet_surv.RData', 
                      './cache/ridge_surv.RData', 
                      './cache/lasso_surv.RData', 
                      './cache/svm_surv.RData', 
                      './cache/svm_imp.RData', 
                      './cache/rf_surv.RData', 
                      './cache/gbm_surv.RData', 
                      './cache/surv_multi.RData'), 
       script_path = c('./survival scripts/elastic_net.R', 
                       './survival scripts/ridge.R', 
                       './survival scripts/lasso.R', 
                       './survival scripts/svm_survival.R', 
                       './survival scripts/svm_importance.R', 
                       './survival scripts/rf.R', 
                       './survival scripts/gbm.R', 
                       './survival scripts/clinic.R'), 
       message = c('Loading chached results of Elastic Net Cox modeling', 
                   'Loading chached results of Ridge Cox modeling', 
                   'Loading chached results of LASSO Cox modeling', 
                   'Loading cached results of SVM modeling', 
                   'Loading cached importance testing for the SVM survival score', 
                   'Loading cached results of RF modeling', 
                   'Loading cached results of GBM modeling', 
                   'Loading cahced results for modeling with clinical variables')) %>% 
    pwalk(access_cache)
  
  ## summary of the stats and plots for the multi-parameter modeling 
  ## testing for the clinic confounder in the pooled GEO cohort
  
  c('./survival scripts/summary.R', 
    './survival scripts/plots.R', 
    './survival scripts/clinic_plots.R', 
    './survival scripts/cohort.R') %>% 
    source_all(message = TRUE, crash = TRUE)

# END -------
  
  insert_tail()