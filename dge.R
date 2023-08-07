# Differential gene expression scripts 

# tools ------
  
  library(plyr)
  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)

  library(microViz)
  library(exda)
  library(nomclust)
  library(clustTools)

  library(AnnotationDbi)
  library(org.Hs.eg.db)

  library(limma)
  library(SPIA)
  library(pathview)
  library(meta)
  library(GOSemSim)
  library(effectsize)

  library(soucer)
  library(furrr)

  insert_head()
  
  select <- dplyr::select
  reduce <- purrr::reduce
  
  source_all('./tools/project_tools.R', 
             message = TRUE, crash = TRUE)

# scripts ------

  insert_msg('Sourcing the scripts')
  
  ## differentially regulated genes, working with cached results
  
  if(file.exists('./cache/dge.RData')) {
    
    insert_msg('Loading cached differential gene expression data')
    
    load('./cache/dge.RData')
    
  } else {
    
    c('./dge scripts/dge_testing.R') %>% 
      source_all(message = TRUE, crash = TRUE)
    
  }
  
  ## differences in signaling pathways: working with cache
  ## because the process is extremely time consuming
  
  if(file.exists('./cache/dge_spia.RData')) {
    
    insert_msg('Loading cached signaling results')
    
    load('./cache/dge_spia.RData')
    
  } else {
    
    c('./dge scripts/dge_pathways.R') %>% 
      source_all(message = TRUE, crash = TRUE)
    
  }
  
  ## detailed analysis with the differentially regulated genes, GO terms
  ## and pathways

  c('./dge scripts/dge_plots.R', ## differential gene expression graphics
    './dge scripts/dge_go.R', ## GO enrichment
    './dge scripts/dge_pathview.R', ## visualization of the common regulated signaling pathways 
    './dge scripts/dge_reactome.R') %>% ## GSVA
    source_all(message = TRUE, crash = TRUE) %>% 
    print
  
# END -----
  
  insert_tail()