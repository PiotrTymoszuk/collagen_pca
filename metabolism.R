# Analysis of biochemical reactions based on differential gene expression 
# between the collagen clusters.

# tools -----

  library(plyr)
  library(tidyverse)
  library(stringi)
  library(rlang)
  library(trafo)
  library(readxl)

  library(BiGGR)
  library(biggrExtra)

  library(exda)
  library(rstatix)
  library(effectsize)

  library(microViz)
  library(gseaTools)

  library(furrr)
  library(soucer)

  insert_head()
  
  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce

  source_all('./tools/project_tools.R', 
             message = TRUE, crash = TRUE)
  
# scripts ------
  
  insert_msg('Soucing the metabolis analysis scripts')
  
  ## identification of significantly regulated reactions
  ## working with the cache, because this analysis takes a while
  
  if(file.exists('./cache/meta.RData')) {
  
    insert_msg('Loading cached metabolic pathway modeling results')
    
    load('./cache/meta.RData')
  
  } else {
    
    c('./metabolism scripts/metab_testing.R') %>% 
      source_all(message = TRUE, crash = TRUE)
    
  }
  
  ## testing for enrichment of the metabolic subsystems
  
  if(file.exists('./cache/meta_sub.RData')) {
    
    insert_msg('Loading cached enrichment results')
    
    load('./cache/meta_sub.RData')
    
  } else {
    
    c('/metabolism scripts/metab_enrichment.R') %>% 
      source_all(message = TRUE, crash = TRUE)
    
  }
  
  ## analyses with the metabolic modeling results
  
  c('./metabolism scripts/metab_plots.R', ## plotting the results of metabolic modeling
    './metabolism scripts/metab_hg.R', ## plots of selected metabolic pathways
    './metabolism scripts/metab_genes.R') %>% ## expression of genes linked to metabolic reactions
    source_all(message = TRUE, crash = TRUE) %>% 
    print

# END -----
  
  insert_tail()