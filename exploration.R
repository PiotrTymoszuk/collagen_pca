# Exploratory analysis of expression of collagen-related genes shared by all 
# transcriptomic datasets
#
# 1) Characteristic of the cohorts and samples
#
# 2) Normality of distribution of the collagen-related gene variables 
# (Shapiro-Wilk test) and their information content 
# (element frequency, variance and Gini coefficient)
#
# 3) Co-expression analysis of the genes in the cancer tissue. This analysis 
# step involves calculation of Euclidean distances between the genes and 
# multi-dimensional scaling as well as by a PCA
#
# 4) Assessment of clustering tendency of the collagen expression data set
# by Hopkins statistic and visualizations (heat maps and  UMAP)

# tools ------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)
  library(readxl)

  library(exda)
  library(microViz)
  library(clustTools)

  library(survival)
  library(survminer)

  library(ggrepel)
  library(ggtext)
  
  library(furrr)
  library(soucer)

  insert_head()
  
  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  set_rownames <- trafo::set_rownames
  
  c('./tools/globals.R', 
    './tools/functions.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# Analysis scripts ---------
  
  insert_msg('Analysis scripts')
  
  ## general distribution stats
  ## coexpression and clustering tendency analysis
  
  c('./exploration scripts/cohorts.R', 
    './exploration scripts/distribution.R', 
    './exploration scripts/coexpression.R', 
    './exploration scripts/clust_tendency.R') %>% 
    source_all(message = TRUE, crash = TRUE) %>% 
    print
  
  ## characteristic of the pooled GEO cohort
  
  c('./exploration scripts/pooled_geo.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END ------
  
  insert_tail()