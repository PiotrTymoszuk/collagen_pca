# Analysis of mutation counts and frequencies in the collagen clusters
# done for the TCGA cohort, Mutect2 algorithm
# the dataset is obtained from XENA 
# (https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PRAD.mutect2_snv.tsv.gz)

# tools ------- 

  library(plyr)
  library(tidyverse)
  library(trafo)
  library(rlang)

  library(xena)
  
  library(exda)
  library(rstatix)
  library(clustTools)
  
  library(soucer)
  library(furrr)

  library(ggrepel)
  library(ggtext)

  library(AnnotationDbi)
  library(org.Hs.eg.db)

  insert_head()

  c('./tools/globals.R', 
    './tools/project_tools.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
  select <- dplyr::select
  reduce <- purrr::reduce
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  mut_globals <- list()
  
  ## cluster assignment 
  
  mut_globals$assignment <- coll_clust$clust_obj$tcga$clust_assignment %>% 
    set_names(c('patient_id', 'clust_id')) %>% 
    mutate(clust_id = factor(clust_id, 
                             rev(levels(coll_clust$clust_obj$tcga$clust_assignment$clust_id))))
  
  ## reading the mutation data/Mutect
  
  mut_globals$mutect_tbl <- 
    load_xena(path = './data/tcga mutations/TCGA-PRAD.mutect2_snv.tsv') %>% 
    filter(!effect %in% c('synonymous_variant', 
                          'intergenic_variant', 
                          'intron_variant')) %>% 
    mutate(patient_id  = stri_extract(sample_id, 
                                      regex = '^\\w{4}-\\w{2}-\\w{4}'))
  
  ## reading the raw copy number variant (CNV) computed
  ## for the TCGA cohort by GISTIC, appending with the official gene symbols
  ## omitting entries without any valid official symbol
  
  mut_globals$gistic_tbl <- 
    read_tsv('./data/tcga mutations/TCGA-PRAD.gistic.tsv') %>% 
    mutate(ensembl_id = `Gene Symbol`, 
           ensembl_id = stri_split_fixed(ensembl_id, 
                                         pattern = '.', 
                                         simplify = TRUE)[, 1], 
           gene_symbol = mapIds(org.Hs.eg.db, 
                                keys = ensembl_id, 
                                column = 'SYMBOL', 
                                keytype = 'ENSEMBL', 
                                multiVals = 'first')) %>% 
    filter(!is.na(gene_symbol), 
           !`Gene Symbol` %in% c('ENSG00000261832.3', 
                                 'ENSG00000249624.7', 
                                 'ENSG00000280987.1', 
                                 'ENSG00000261130.4', 
                                 'ENSG00000272916.4', 
                                 'ENSG00000254093.7', 
                                 'ENSG00000279195.1', 
                                 'ENSG00000258790.1', 
                                 'ENSG00000258653.3', 
                                 'ENSG00000275740.1', 
                                 'ENSG00000255330.7', 
                                 'ENSG00000273167.1', 
                                 'ENSG00000255508.6', 
                                 'ENSG00000270011.5')) %>% 
    relocate(gene_symbol)
  
  ## manual correction of unambiguous gene symbols
  
  mut_globals$symbol_conv <-
    c('ENSG00000228696.7' = 'ARL17B', 
      'ENSG00000169627.7' = 'BOLA2B', 
      'ENSG00000226023.5' = 'CT47A6', 
      'ENSG00000176797.3' = 'DEF103A', 
      'ENSG00000186599.7' = 'DEF105B', 
      'ENSG00000198129.2' = 'DEF107B', 
      'ENSG00000177257.2' = 'DEF4B', 
      'ENSG00000256966.5' = 'H3BUX3', 
      'ENSG00000278662.3' = 'GOLGA6L10', 
      'ENSG00000237541.3' = 'HLA-DQA2', 
      'ENSG00000178934.4' = 'LGALS7B', 
      'ENSG00000221995.5' = 'TIAF1', 
      'ENSG00000111215.10' = 'PRR', 
      'ENSG00000205572.8' = 'SERF1B', 
      'ENSG00000205571.11' = 'SMN2', 
      'ENSG00000274512.3' = 'TBC1D3L', 
      'ENSG00000269226.6' = 'TSMB15C', 
      'ENSG00000236125.3' = 'USP17L4')
  
  for(i in names(mut_globals$symbol_conv)) {
    
    mut_globals$gistic_tbl[mut_globals$gistic_tbl$`Gene Symbol` == i, 'gene_symbol'] <- 
      mut_globals$symbol_conv[i]
    
  }

  rm(i)
  
# Analysis scripts ------
  
  insert_msg('Analysis scripts')
  
  ## cached mutation indexes, frequencies and counts
  
  if(file.exists('./cache/mut_tables.RData')) {
    
    insert_msg('Loading cached mutation tables')
    
    load('./cache/mut_tables.RData')
    
  } else {
    
    c('./mutation scripts/tables.R') %>% 
      source_all(message = TRUE, crash = TRUE)
    
  }
  
  c('./mutation scripts/counts.R',
    './mutation scripts/mutations.R', 
    './mutation scripts/cnv.R', 
    './mutation scripts/collagen.R') %>% 
    source_all(message = TRUE, crash = TRUE) %>% 
    print
  
# END -----
  
  insert_tail()