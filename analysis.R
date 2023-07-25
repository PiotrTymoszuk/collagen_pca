# Sources analysis scripts: 
# explorative data analysis and survival modeling

# tools ------

  library(plyr)
  library(tidyverse)
  library(trafo)
  library(rlang)
  library(stringi)

  library(ggtext)
  library(ggrepel)
  library(microViz)

  library(rsbml)
  library(exda)
  library(glmnet)
  library(caret)
  library(caretExtra)
  
  library(clustTools)
  library(mclust)
  library(sva)

  library(survival)
  library(survminer)
  library(kmOptimizer)
  library(coxExtensions)

  library(soucer)
  library(furrr)

  insert_head()
  
  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  map <- purrr::map
  cv <- clustTools::cv
  extract <- clustTools::extract
  
  
  source_all('./tools/project_tools.R', 
             crash = TRUE, message = TRUE)
  
# scripts ------
  
  insert_msg('soucing analysis scripts')
  
  ## characteristic of the cohorts and collagen-associated genes
  
  c('./analysis scripts/cohorts.R', ## descriptive characteristic of the study cohorts
    './analysis scripts/normal_tumor.R', ## expression of collagen genes in the tumor and normal prostate
    './analysis scripts/correlation.R') %>%  ## pairwise correlation of collagen pathway gene expression
    source_all(message = TRUE, crash = TRUE) %>% 
    print
  
  ## development of the survival-predicting Collagen Score
  ## association of the collagen gene expression with survival, clinics
  ## and infiltration
  
  c('./analysis scripts/collagen_score.R', ## development and properties of the Collagen Score
    './analysis scripts/collagen_score_plots.R', ## graphs for features of the Collagen Score
    './analysis scripts/clinic.R', ## Collagen genes and Collagen Score as a function of clinical features
    './analysis scripts/survival_uni.R', ## Uni-variable Cox modeling of survival
    './analysis scripts/survival_cut.R', ## survival in dichotomous collagen gene and Collagen Score strata
    './analysis scripts/survival_multi.R', ## multi-variable survival modeling, Collagen Score
    './analysis scripts/collagen_score_infiltration.R') %>% ## correlation of the Collagen Score and non-malignant cell infiltration)
    source_all(message = TRUE, crash = TRUE) %>% 
    print
  
  ## development and characteristic of the Collagen Clusters
  
  c('./analysis scripts/clust_devel.R', ## development of the collagen clusters
    './analysis scripts/collagen_clustering.R', ## semi-supervised clustering
    './analysis scripts/collagen_cluster_survival.R', ## diffferences in survival between the Collagen Clusters
    './analysis scripts/collagen_cluster_clinic.R', ## differences in clinics and Collagen score between the clusters
    './analysis scripts/collagen_cluster_infiltration.R') %>% ## infiltration in the clusters
    source_all(message = TRUE, crash = TRUE) %>% 
    print

# END -----
  
  insert_tail()
