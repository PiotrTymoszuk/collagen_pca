# Association of the collagen-related gene expression with features of the tumor, 
# such as tissue type and aggressiveness measured by Gleason scoring 
# stratified by the ISUP rules (5 - 6, 7, 8+).
#
# 1) Analysis of the tumor and normal tissue stiffness by MRI. The diffusion
# capacities as hallmarks of the collagen content are compared be paired T test 
# with Cohen's d in paired measurements of the adjacent prostate and the tumor.
#
# 2) Differential expression of the collagen genes in the tumor and normal lung. 
# This is checked in paired specimens of the TCGA and GSE70768 data sets by 
# paired T test and Chen's d statistic.
#
# 3) Univariable comparison of gene expression between the ISUP strata of 
# Gleason scores with one-way ANOVA with eta-square effect size metric.
#
# 4) Ordinal Elastic Net regression of collagen gene expression in the Gleason 
# score strata - just a trial, since no meaningful multi-gene signature of the 
# Gleason score strata (ISUP criteria) could be established.
#
# 5) Correlation analysis (Spearman) of PSA at diagnosis and collagen gene 
# expression


# tools -------

  library(tidyverse)
  library(trafo)
  library(rlang)
  library(stringi)

  library(exda)
  library(microViz)
  library(caret)
  library(caretExtra)
  
  library(doParallel)
  library(furrr)
  library(soucer)

  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  set_rownames <- trafo::set_rownames
  extract <- clustTools::extract
  
  insert_head()
  
  c('./tools/globals.R', 
    './tools/functions.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# Analysis scripts --------
  
  insert_msg('Analysis scripts')
  
  c('./gleason scripts/mri.R', 
    './gleason scripts/normal_tumor.R', 
    './gleason scripts/univariable.R', 
    './gleason scripts/psa.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
  ## establishing a multi-gene signature
  
  access_cache(cache_path = './cache/gs_sign.RData', 
               script_path = './gleason scripts/signature.R', 
               message = 'Loading Elastic Net modeling results')
  
# END ------
  
  insert_tail()
