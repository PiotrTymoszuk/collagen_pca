# Reactome signatures

  insert_head()

# container ------
  
  reactome <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# signature database -------
  
  insert_msg('Signature database')

  reactome$db <- load_dbsig('./data/signatures/msigdb.v7.5.1.symbols.gmt')
  
  reactome$db <- reactome$db %>% 
    filter(stri_detect(sign_name, regex = '^REACTOME'))
  
# Expression data frames ------
  
  insert_msg('Expression data frames')

  reactome$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor')
  
  ## frames to be used later for annotation of the infiltration estmates
  
  reactome$annotation <- reactome$expression %>% 
    map(select, sample_id, patient_id)
  
  reactome$expression <- reactome$expression %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, -patient_id, -tissue_type)
  
# ssGSEA scores ------
  
  insert_msg('ssGSEA scores')
  
  reactome$signatures <- reactome$expression %>% 
    future_map(~calculate.dbsig(reactome$db, data = .x))
  
  reactome$signatures <- 
    map2(reactome$annotation, reactome$signatures, cbind) %>% 
    map(as_tibble)
  
# signature lexicon -------
  
  insert_msg('Signature lexicon')
  
  reactome$cmm_signatures <- reactome$signatures %>% 
    map(select, -sample_id, -patient_id) %>% 
    map(names) %>% 
    reduce(intersect)
  
  reactome$lexicon <- reactome$db %>% 
    transmute(variable = sign_name, 
              label = stri_replace(variable, 
                                   regex = '^REACTOME_', 
                                   replacement = ''), 
              label = stri_replace_all(label, 
                                       fixed = '_', 
                                       replacement = ' '), 
              genes = genes) %>% 
    filter(variable %in% reactome$cmm_signatures)
  
# caching the data sets ------
  
  insert_msg('Caching the datasets')
  
  reactome <- reactome[c("signatures", "lexicon")]
  
  save(reactome, file = './data/reactome.RData')
  
# END -----
  
  plan('sequential')
  
  insert_tail()