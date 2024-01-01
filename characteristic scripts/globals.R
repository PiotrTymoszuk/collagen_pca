# Analysis globals: cluster assignment, whole transcriptome expression and 
# n numbers of observations in the clusters. 

  insert_head()
  
# container -------
  
  ana_globals <- list()
  
# cluster assignment -------
  
  insert_msg('Cluster assignment')
  
  ana_globals$assignment <- clust_semi$assignment
  
# transcriptome -------
  
  insert_msg('Whole transcriptome expression')
  
  ## gene variables
  
  ana_globals$annotation <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$annotation) %>% 
    map(filter, 
        !duplicated(gene_symbol), 
        !duplicated(entrez_id))
  
  ana_globals$genes <- ana_globals$annotation %>% 
    map(~.x$gene_symbol)
  
  ## log2-expression
  
  ana_globals$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map2(ana_globals$genes, 
         ~.x[c('sample_id', .y)]) %>% 
    map2(ana_globals$assignment, ., 
         left_join, by = 'sample_id') %>% 
    map(filter, !is.na(clust_id))
  
# N numbers of observations in the clusters -------
  
  insert_msg('N numbers of observations in the clusters')
  
  ana_globals$n_numbers <- clust_semi$clust_obj %>% 
    map(ngroups)
  
# END --------
  
  insert_tail()