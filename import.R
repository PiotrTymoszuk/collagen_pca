# This script imports prostate cancer data from GEO (GSE16560, GSE40272, GSE70768 and GSE70769)
# and the TCGA PCA next generation RNA seq data
# Expression values log2(x) (Microarray) or log2(x + 1) (TCGA) transformed
#
# For modeling and clustering, the collagen gene dataset is normalized 
# with ComBat to suppress the cohort-specific effects

# toolbox ----

  library(plyr)
  library(tidyverse)
  library(readxl)
  library(purrr)
  library(rlang)
  library(trafo)
  
  library(immunedeconv)
  library(sva)
  library(exda)
  library(clustTools)
  
  library(furrr)
  library(doParallel)
  library(soucer)

  c('./tools/import_engine.R', 
    './tools/tcga_tools.R', 
    './tools/project_tools.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
  insert_head()
  
# data containers -----
  
  geo_data <- list() ## raw geo data
  tcga_data <- list() ## raw tcga data

  study_data <- list() ## cleared data in the uniform formt
  
# executing data import scripts, from scratch if not done before, saving the raw data ----
  
  insert_msg('Reading the raw expression and clinical data')
  
  if(('geo_data.RData' %in% list.files('./data')) & 
     ('tcga_data.RData' %in% list.files('./data'))) {
    
    load('./data/geo_data.RData')
    load('./data/tcga_data.RData')

  } else {
    
    c('./import scripts/geo_import.R', 
      './import scripts/tcga_import.R') %>% 
      walk(function(x) tryCatch(source(x), 
                                finally = go_proj_directory()))
    
    save(geo_data, file = './data/geo_data.RData')
    save(tcga_data, file = './data/tcga_data.RData')
    
  }
  
# data clearing: creating uniform data tibbles with whole genome expression and clinical data ----- 
  
  insert_msg('Data clearing')
  
  if('study_data.RData' %in% list.files('./data')) {
    
    load('./data/study_data.RData')

  } else {
    
    c('./clearing scripts/data_clearing_geo.R', 
      './clearing scripts/data_clearing_tcga.R') %>% 
      source_all(message = TRUE, crash = TRUE)
    
    save(study_data, file = './data/study_data.RData')
    
  }
  
# Fetching the infiltration data -------
  
  insert_msg('Infiltration data')
  
  if(file.exists('./data/infiltration.RData')) {
    
    load('./data/infiltration.RData')
    
  } else  {
    
    source_all('./clearing scripts/infiltration.R', 
               message = TRUE, crash = TRUE)
    
  }
  
  infil$xcell_types <- infil$xcell[[1]]$cell_type
  infil$mcp_counter_types <- infil$mcp_counter[[1]]$cell_type
  
  infil[c('xcell', 'mcp_counter')] <- 
    infil[c('xcell', 'mcp_counter')] %>% 
    map(~map(.x, column_to_rownames, 'cell_type') %>% 
          map(t) %>% 
          map(as.data.frame) %>% 
          map(rownames_to_column, 'patient_id') %>% 
          map(as_tibble))
  
# Fetching the Reactome signatures ------
  
  insert_msg('Fetching the Reactome signatures')
  
  if(file.exists('./data/reactome.RData')) {
    
    insert_msg('Loading the cached signature dataset')
    
    load('./data/reactome.RData')
    
  } else {
    
    source_all('./clearing scripts/reactome.R', 
               message = TRUE, crash = TRUE)
    
  }

# globals setup -----
  
  insert_msg('Globals setup')
  
  c('./tools/globals.R') %>% 
    source_all(message = TRUE, crash = TRUE)

# ComBat normalization of the collagen dataset ------
  
  insert_msg('Combat normalization')
  
  c('./clearing scripts/combat.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END -----
  
  insert_tail()
