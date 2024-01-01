# Characteristic of the collagen clusters

# tools --------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)

  library(clustTools)
  library(microViz)
  library(exda)

  library(limma)
  library(decoupleR)
  library(SPIA)
  library(biggrExtra)

  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GOSemSim)

  library(survival)
  library(survminer)

  library(pathview)

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
  
# Analysis globals -------
  
  insert_msg('Analysis globals')
  
  source_all('./characteristic scripts/globals.R', 
             message = TRUE, crash = TRUE)
  
# Clinical characteristic and survival -------
  
  insert_msg('Clinical characteristic and survival')
  
  c('./characteristic scripts/clinic.R', 
    './characteristic scripts/os.R', 
    './characteristic scripts/rfs.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# Immunity -------
  
  insert_msg('Immunity and microenvironment')
  
  c('./characteristic scripts/xcell.R', 
    './characteristic scripts/mcp_counter.R') %>% 
    source_all(message = TRUE, crash = TRUE) %>% 
    print
  
# Differential gene expression, biology and signaling --------
  
  insert_msg('Biology and differential gene expression')
  
  ## gene set variation analyses
  
  c('./characteristic scripts/reactome.R', 
    './characteristic scripts/recon.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
  ## analyses based on differential gene expression
  
  list(cache_path = c('./cache/ana_dge.RData', 
                      './cache/ana_spia.RData', 
                      './cache/ana_go.RData', 
                      './cache/ana_meta.RData'), 
       script_path = c('./characteristic scripts/dge.R',
                       './characteristic scripts/spia.R', 
                       './characteristic scripts/go.R', 
                       './characteristic scripts/metabolism.R'), 
       message = c('Cached differential gene expression analysis', 
                   'Cached signaling pathway analysis with SPIA', 
                   'Cached GO enrichment results', 
                   'Cached Monte Carlo modeling of metabolic reactions')) %>% 
    pwalk(access_cache)
  
  ## decoupleR analyses
  
  c('./characteristic scripts/progeny.R', 
    './characteristic scripts/collectri.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
  ## visualization of the differential gene expression, signaling 
  ## and metabolic modeling
  
  c('./characteristic scripts/dge_plots.R', 
    './characteristic scripts/pathview.R', 
    './characteristic scripts/metabolism_plots.R', 
    './characteristic scripts/metabolic_genes.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# Mutations ------
  
  insert_msg('Differences in mutation counts and rates')
  
  c('./characteristic scripts/tmb.R', 
    './characteristic scripts/mutations.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END ------
  
  insert_tail()