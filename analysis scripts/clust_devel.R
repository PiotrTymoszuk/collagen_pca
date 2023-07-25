# Clustering of the participants: finding the best-performing clustering 
# algorithm in the training TCGA cohort
#
# Working with the Combat-normalized collagen gene dataset

  insert_head()
  
# container ----
  
  clust_dev <- list()
  
# Parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals ------
  
  insert_msg('Globals setup')
  
  ## analysis tables with mean-centered collagen genes
  
  clust_dev$variables <- globals$genes_interest$gene_symbol

  clust_dev$analysis_tbl <- combatch$adjusted_data %>% 
    map(~.x[c('patient_id', clust_dev$variables)]) %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(column_to_rownames, 'patient_id') %>% 
    map(center_data, type = 'mean')
  
  ## clustering distances
  
  clust_dev$dists <- c('euclidean', 
                       'sumofsquares', 
                       'manhattan', 
                       'cosine')
  
# clustering tendency of the data sets ------
  
  insert_msg('Clustering tendencies')
  
  clust_dev$clust_tend <- 
    list(data = clust_dev$analysis_tbl, 
         n = floor(0.5 * map_dbl(clust_dev$analysis_tbl, nrow))) %>% 
    future_pmap(get_clust_tendency, 
                .options = furrr_options(seed = TRUE))
  
# clustering algorithms -----
  
  insert_msg('Building clustering objects')
  
  ## Ward's hierarchical clustering

  clust_dev$algos[paste0('hcl_', clust_dev$dists)] <- 
    clust_dev$dists %>% 
    map(~hcluster(data = clust_dev$analysis_tbl$tcga, 
                  distance_method = .x,
                  k = 3, 
                  hc_method = 'ward.D2'))
  
  ## K-MEANS
  
  clust_dev$algos[paste0('kmeans_', clust_dev$dists)] <- 
    clust_dev$dists %>% 
    map(~kcluster(data = clust_dev$analysis_tbl$tcga, 
                  distance_method = .x,
                  k = 3, 
                  clust_fun = 'kmeans'))
  
  ## PAM
  
  clust_dev$algos[paste0('pam_', clust_dev$dists)] <- 
    clust_dev$dists %>% 
    map(~kcluster(data = clust_dev$analysis_tbl$tcga, 
                  distance_method = .x,
                  k = 3, 
                  clust_fun = 'pam'))

# Variances of the algorithms ------
  
  insert_msg('Clustering variances')
  
  clust_dev$variance <- clust_dev$algos %>% 
    map(clustTools::var) %>% 
    map_dbl(~.x$frac_var) %>% 
    compress(names_to = 'algorithm', 
             values_to = 'clust_variance')

# Cross-validation ------
  
  insert_msg('Cross-validation')

  clust_dev$cv <- clust_dev$algos %>% 
    future_map(cv.clust_analysis, 
               nfolds = 10, 
               kNN = 7, 
               simple_vote = FALSE, 
               resolve_ties = TRUE, 
               .parallel = FALSE, 
               .options = furrr_options(seed = TRUE))
  
  clust_dev$cv <- clust_dev$cv %>% 
    map(~.x$summary) %>% 
    compress(names_to = 'algorithm')

# common table with the results and plotting ------
  
  insert_msg('Commmon result table and plotting')
  
  clust_dev$results <- left_join(clust_dev$variance, 
                                 clust_dev$cv, 
                                 by = 'algorithm') %>% 
    mutate(accuracy = 1 - mean_error, 
           algo_lab = toupper(stri_extract(algorithm, regex = 'pam|kmeans|hcl')), 
           algo_lab = paste(algo_lab, 
                            stri_split_fixed(algorithm, 
                                             pattern = '_', 
                                             simplify = TRUE)[, 2], 
                            sep = ', '))
  
  clust_dev$plot <- clust_dev$results[c('algo_lab', 
                                        'clust_variance', 
                                        'accuracy')] %>% 
    pivot_longer(cols = c('clust_variance', 'accuracy'), 
                 names_to = 'statistic', 
                 values_to = 'value') %>% 
    ggplot(aes(x = value, 
               y = reorder(algo_lab, value), 
               fill = statistic)) + 
    geom_bar(stat = 'identity', 
             color = 'black', 
             position = position_dodge(0.9)) + 
    scale_fill_manual(values = c(clust_variance = 'darkolivegreen4', 
                                 accuracy = 'steelblue3'), 
                      labels = c(clust_variance = 'Clustering variance', 
                                 accuracy = 'Accuracy, 10-fold CV'), 
                      name = '') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Clustering algorithm performance', 
         subtitle = 'TCGA training cohort', 
         x = 'Statistic value')
      
# END ------
  
  plan('sequential')