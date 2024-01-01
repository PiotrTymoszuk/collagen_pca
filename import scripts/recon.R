# Recon subsystem signatures

  insert_head()

# container ------
  
  recon <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# signature database -------
  
  insert_msg('Signature database')

  recon$reactions <- extract_subsystems(Recon2D)
  
  recon$genes <- extract_genes(Recon2D)
  
  recon$db <- recon[c("genes", "reactions")] %>% 
    reduce(left_join, by = 'react_id') %>% 
    blast(subsystem) %>% 
    map(~.x$entrez_id) %>% 
    map(reduce, union)
  
  ## conversion to official gene symbols
  
  recon$db <- recon$db %>% 
    map(mapIds, 
        x = org.Hs.eg.db, 
        keytype = 'ENTREZID', 
        column = 'SYMBOL')
  
# Expression data frames ------
  
  insert_msg('Expression data frames')

  recon$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor')
  
  ## frames to be used later for annotation of the infiltration estmates
  
  recon$annotation <- recon$expression %>% 
    map(select, sample_id, patient_id)
  
  recon$expression <- recon$expression %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, -patient_id, -tissue_type)
  
# ssGSEA scores ------
  
  insert_msg('ssGSEA scores')
  
  recon$signatures <- recon$expression %>% 
    future_map(~calculate.default(recon$db, data = .x))
  
  recon$signatures <- 
    map2(recon$annotation, recon$signatures, cbind) %>% 
    map(as_tibble)
 
  recon$signatures <- recon$signatures %>% 
    map(~set_colnames(.x, make.names(names(.x))))
  
# signature lexicon -------
  
  insert_msg('Signature lexicon')
  
  recon$cmm_signatures <- recon$signatures %>% 
    map(select, -sample_id, -patient_id) %>% 
    map(names) %>% 
    reduce(intersect)
  
  recon$lexicon <- recon$db %>% 
    tibble(variable = make.names(names(.)), 
           label = names(.), 
           genes = .) %>% 
    filter(variable %in% recon$cmm_signatures)
  
# caching the data sets ------
  
  insert_msg('Caching the datasets')
  
  recon <- recon[c("signatures", "lexicon")]
  
  save(recon, file = './data/recon.RData')
  
# END -----
  
  plan('sequential')
  
  insert_tail()