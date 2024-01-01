# Distribution tests for the collagen variables. 
# Done for Z-scores, since they will be used for MDS, PCA, clustering tendency
# assessment, clustering and survival modeling.

  insert_head()
  
# container --------
  
  distr <- list()
  
# parallel backend --------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# expression data frames -------
  
  insert_msg('Expression data')
  
  distr$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(globals$genes_interest$gene_symbol)) %>% 
    map(center_data, 'mean')
  
# Normality by Shapiro-Wilk test --------
  
  insert_msg('Normality by Shapiro-Wilk test')
  
  distr$normality <- distr$expression %>% 
    future_map(explore, 
               what = 'normality', 
               .options = furrr_options(seed = TRUE)) %>% 
    map(arrange, test_stat)
  
# Variance and information content ------
  
  insert_msg('Variance and Gini coefficients')
  
  distr$stats <- distr$expression %>% 
    map(minimum_shift) %>% 
    map(distr_stats) %>% 
    map(arrange, gini_coef)
    
# END ------
  
  distr <- distr[c("normality", "stats")]
  
  plan('sequential')
  
  insert_tail()