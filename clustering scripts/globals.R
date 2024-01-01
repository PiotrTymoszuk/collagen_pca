# Clustering globals: identity expression and Z-scores of 
# the collagen-related genes 

  insert_head()
  
# container --------
  
  clust_globals <- list()
  
# collagen gene expression --------
  
  insert_msg('Collagen gene expression')
  
  ## genes
  
  clust_globals$variables <- globals$genes_interest$gene_symbol
  
  ## identity expression
  
  clust_globals$data_identity <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(select, 
        sample_id, 
        all_of(clust_globals$variables)) %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(column_to_rownames, 'sample_id')
  
  ## Z-scores
  
  clust_globals$data <- clust_globals$data_identity %>% 
    map(center_data, 'mean')
  
# clustering distances ------
  
  insert_msg('Distance metrics')
  
  clust_globals$dists <- c('euclidean', 
                       'squared_euclidean', 
                       'manhattan', 
                       'cosine')
  
# N numbers ----------
  
  insert_msg('N numbers')
  
  clust_globals$n_numbers <- clust_globals$data %>% 
    map(nrow)
  
  clust_globals$cohort_caps <- clust_globals$n_numbers %>% 
    map2_chr(names(.), ., 
             ~paste0(globals$study_labels[.x], '\nn = ', .y)) %>% 
    set_names(names(clust_globals$data))
  
# END ------
  
  insert_tail()