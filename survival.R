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
  
  library(furrr)
  library(soucer)

  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  set_rownames <- trafo::set_rownames
  
  c('./tools/globals.R', 
    './tools/functions.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# Univariable analysis ---------
  
  insert_msg('Univariable analysis')
  
  ## resorting to cached versions, because the analyses take a while
  
  list(cache_path = c('./cache/os_cut.RData', 
                      './cache/rfs_cut.RData'), 
       script_path = c('./survival scripts/univariable_os.R', 
                       './survival scripts/univariable_rfs.R'), 
       message = paste('Loading cached univariable survival', 
                       'analysis results for', 
                       c('OS', 'RFS'))) %>% 
    pwalk(access_cache)
  
# Elastic Net Cox modeling ---------
  
  insert_msg('Elastic Net Cox modeling')
  
  access_cache(cache_path = './cache/coll_score.RData', 
               script_path = './survival scripts/collagen_score.R', 
               message = 'Loading chached results of Elastic Net Cox modeling')
  
  c('./survival scripts/collagen_score_plots.R', 
    './survival scripts/collagen_score_multi.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END -------
  
  insert_tail()