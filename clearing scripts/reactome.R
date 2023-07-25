# Reactome signatures

  insert_head()
  
# tools -----
  
  library(gseaTools)
  
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
  
  reactome$expression <- study_data %>% 
    map(~.x$expression)
  
  for(i in names(reactome$expression)) {
    
    if(!'tissue' %in% names(reactome$expression[[i]])) next

    reactome$expression[[i]] <- reactome$expression[[i]] %>% 
      filter(tissue == 'tumor')
    
  }
  
  reactome$ids <- reactome$expression %>% 
    map(~.x$patient_id)
  
  reactome$expression <- 
    map2(reactome$expression, 
         map(study_data, ~.x$annotation$symbol), 
         ~.x[.y])
  
# ssGSEA scores ------
  
  insert_msg('ssGSEA scores')
  
  reactome$signatures <- reactome$expression %>% 
    future_map(~calculate.dbsig(reactome$db, data = .x))
  
  reactome$signatures <- 
    map2(reactome$signatures, 
         reactome$ids, 
         ~mutate(.x, patient_id = .y)) %>% 
    map(relocate, patient_id)
  
# signature lexicon -------
  
  insert_msg('Signature lexicon')
  
  reactome$cmm_signatures <- reactome$signatures %>% 
    map(select, -patient_id) %>% 
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
  
# caching the datasets ------
  
  insert_msg('Caching the datasets')
  
  reactome <- reactome[c("signatures", "lexicon")]
  
  save(reactome, file = './data/reactome.RData')
  
# END -----
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()