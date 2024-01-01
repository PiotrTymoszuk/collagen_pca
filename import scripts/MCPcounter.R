# Immune deconvolution with mcp

  insert_head()
  
# container -------
  
  mcp <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')

  plan('multisession')
    
# expression data -------
  
  insert_msg('Expression data')
  
  mcp$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor')
  
  ## frames to be used later for annotation of the infiltration estmates
  
  mcp$annotation <- mcp$expression %>% 
    map(select, sample_id, patient_id)
  
  ## 2^pow transformation and matrix conversion
  
  mcp$expression <- mcp$expression %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, -patient_id, -tissue_type) %>% 
    map(as.matrix)
  
  mcp$expression <- mcp$expression %>% 
    map(~2^.x) %>% 
    map(t)
  
# immune deconvolution -------
  
  insert_msg('Immune deconvolution')
  
  mcp$cell_fractions <- 
    list(gene_expression = mcp$expression,
         arrays = globals$study_arrays[names(mcp$expression)]) %>%
    future_pmap(deconvolute,
                method = 'mcp_counter')
  
# Formatting of the results --------
  
  insert_msg('Formatting the ouput')
  
  mcp$cell_fractions <- mcp$cell_fractions %>% 
    map(column_to_rownames, 'cell_type') %>% 
    map(t) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map2(mcp$annotation, ., 
         left_join, by = 'sample_id')
  
# Caching the results --------
  
  insert_msg('Caching the results')
  
  mcp <- mcp$cell_fractions
  
  save(mcp, file = './data/mcp.RData')
  
# END -----
  
  plan('sequential')
  
  insert_tail()