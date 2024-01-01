# Immune deconvolution with xCell

  insert_head()
  
# container -------
  
  xcell <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')

  plan('multisession')
    
# expression data -------
  
  insert_msg('Expression data')
  
  xcell$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor')
  
  ## frames to be used later for annotation of the infiltration estmates
  
  xcell$annotation <- xcell$expression %>% 
    map(select, sample_id, patient_id)
  
  ## 2^pow transformation and matrix conversion
  
  xcell$expression <- xcell$expression %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, -patient_id, -tissue_type) %>% 
    map(as.matrix)
  
  xcell$expression <- xcell$expression %>% 
    map(~2^.x) %>% 
    map(t)
  
# immune deconvolution -------
  
  insert_msg('Immune deconvolution')
  
  xcell$cell_fractions <- 
    list(gene_expression = xcell$expression,
         arrays = globals$study_arrays[names(xcell$expression)]) %>%
    future_pmap(deconvolute,
                method = 'xcell',
                expected_cell_types = c('B cell', 
                                        'T cell CD8+', 
                                        'T cell CD4+ (non-regulatory)', 
                                        'T cell regulatory (Tregs)', 
                                        'Myeloid dendritic cell', 
                                        'Macrophage', 
                                        'Monocyte', 
                                        'Neutrophil', 
                                        'NK cell', 
                                        'Cancer associated fibroblast', 
                                        'Endothelial cell'))
  
# Formatting of the results --------
  
  insert_msg('Formatting the ouput')
  
  xcell$cell_fractions <- xcell$cell_fractions %>% 
    map(column_to_rownames, 'cell_type') %>% 
    map(t) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map2(xcell$annotation, ., 
         left_join, by = 'sample_id')
  
# Caching the results --------
  
  insert_msg('Caching the results')
  
  xcell <- xcell$cell_fractions
  
  save(xcell, file = './data/xcell.RData')
  
# END -----
  
  plan('sequential')
  
  insert_tail()