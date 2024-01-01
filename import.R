# This script imports prostate cancer data from GEO 
# (GSE16560, GSE54460, GSE70768 and GSE70769, GSE116918)
# and the TCGA and DKFZ RNA seq data obteined from cBioportal
# Expression values log2(x) (Microarray) or log2(x + 1) (RNAseq) transformed

# toolbox ----

  library(tidyverse)
  library(readxl)
  library(rlang)
  library(trafo)
  library(stringi)
  
  library(immunedeconv)
  library(sva)
  library(exda)
  library(microViz)
  library(gseaTools)
  
  library(GEOquery)
  library(biggrExtra)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
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
  
# executing data import scripts, from scratch if not done before, saving the raw data ----
  
  insert_msg('Reading the expression and clinical data')

  for(i in names(globals$study_labels)) {
    
    access_cache(cache_path = paste0('./data/', i, '.RData'), 
                 script_path = paste0('./import scripts/', 
                                      globals$study_labels[i], '.R'),
                 message = paste('Loading cleared data for:', i))
    
  }
  
  ## control for presence of all requested collagen-related transcripts
  
  globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(names) %>% 
    map(~globals$genes_interest$gene_symbol[!globals$genes_interest$gene_symbol %in% .x])
  
# Infiltration --------
  
  insert_msg('Infiltration')
  
  list(cache_path = c('./data/xcell.RData', 
                      './data/mcp.RData'), 
       script_path = c('./import scripts/xCell.R', 
                       './import scripts/MCPcounter.R'), 
       message = paste('Loading cached infiltration data for:', 
                       c('xCell', 'MCP counter'))) %>% 
    pwalk(access_cache)
  
# Reactome and recon signatures -------
  
  insert_msg('Signature scores')
  
  list(cache_path = c('./data/reactome.RData',
                      './data/recon.RData'), 
       script_path = c('./import scripts/reactome.R', 
                       './import scripts/recon.R'), 
       message = paste('Loading cached ssGSEA scores for', 
                       c('Reactome pathways', 
                         'Recon subsystems'))) %>% 
    pwalk(access_cache)
  
# Treatment of genes used for clustering with COMBAT -------
  
  insert_msg('COMBAT for the clustering genes')
  
  ## done only on request: 
  ## abstaining from the batch correction since it distorts 
  ## the gene co-expression patterns (nearest neighbors)
  ## as investigated by MDS and PCA during exploratory data analysis
  
  #source_all('./import scripts/combat.R', 
   #          message = TRUE, crash = TRUE)

# END -----
  
  insert_tail()
