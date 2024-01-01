# Permutation importance of the clustering factors

  insert_head()
  
# container -------
  
  clust_imp <- list()
  
# importance -------
  
  insert_msg('Importance')
  
  clust_imp$test <- clust_semi$clust_obj$tcga %>% 
    impact(n_iter = 100, 
           seed = 12345, 
           .parallel = TRUE)
  
# Plotting --------
  
  insert_msg('Plotting')
  
  clust_imp$plot <- clust_imp$test %>% 
    plot.importance(plot_title = paste('Importance of clustering factors,', 
                                       globals$study_labels["tcga"]), 
                    plot_subtitle = 'Permutation importance, n = 100 iterations', 
                    cust_theme = globals$common_theme, 
                    point_size = 1, 
                    point_alpha = 0.25) + 
    theme(axis.text.y = element_text(face = 'italic'))
  
# Caching the results ------
  
  insert_msg('Caching the results')
  
  save(clust_imp, file = './cache/clust_imp.RData')
  
# END -----
  
  insert_tail()