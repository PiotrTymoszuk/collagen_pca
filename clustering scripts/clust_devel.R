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

# clustering algorithms -----
  
  insert_msg('Building clustering objects')
  
  ## Ward's hierarchical clustering

  clust_dev$algos[paste0('hcl_', clust_globals$dists)] <- 
    clust_globals$dists %>% 
    map(~hcluster(data = clust_globals$data$tcga, 
                  distance_method = .x,
                  k = 2, 
                  hc_method = 'ward.D2'))
  
  ## K-MEANS
  
  clust_dev$algos[paste0('kmeans_', clust_globals$dists)] <- 
    clust_globals$dists %>% 
    map(~kcluster(data = clust_globals$data$tcga, 
                  distance_method = .x,
                  k = 2, 
                  clust_fun = 'kmeans'))
  
  ## PAM
  
  clust_dev$algos[paste0('pam_', clust_globals$dists)] <- 
    clust_globals$dists %>% 
    map(~kcluster(data = clust_globals$data$tcga, 
                  distance_method = .x,
                  k = 2, 
                  clust_fun = 'pam'))

# Numeric stats in the training data set ------
  
  insert_msg('Training data set stats')
  
  clust_dev$train_stats <- clust_dev$algos %>% 
    map(summary) %>% 
    compress(names_to = 'algorithm') %>% 
    mutate(dataset = 'training')
  
# Cross-validation ------
  
  insert_msg('Cross-validation')

  clust_dev$cv_stats <- clust_dev$algos %>% 
    future_map(cv, 
               nfolds = 5, 
               kNN = 11, 
               simple_vote = FALSE, 
               resolve_ties = TRUE, 
               .parallel = FALSE, 
               .options = furrr_options(seed = TRUE))
  
  clust_dev$cv_stats <- clust_dev$cv_stats %>% 
    map(summary) %>% 
    compress(names_to = 'algorithm') %>% 
    mutate(dataset = 'cv') %>% 
    select(algorithm, dataset, ends_with('mean'))
  
  names(clust_dev$cv_stats) <- names(clust_dev$cv_stats) %>% 
    stri_replace(fixed = '_mean', replacement = '')

# common table with the results ------
  
  insert_msg('Commmon result table')
  
  clust_dev$stats <- full_rbind(clust_dev$train_stats, 
                                clust_dev$cv_stats) %>% 
    mutate(algo_lab = toupper(stri_extract(algorithm, regex = 'pam|kmeans|hcl')), 
           dist_lab = stri_extract(algorithm, 
                                   regex = paste(clust_globals$dists,
                                                 collapse = '|')), 
           dist_lab = stri_replace(dist_lab, fixed = '_', replacement = ' '), 
           algo_lab = paste(algo_lab, dist_lab, sep = ', '))
  
# Bar plot of the numeric statistics --------
  
  insert_msg('Plotting')
  
  clust_dev$plot_data <- clust_dev$stats %>% 
    mutate(frac_flaw_neighborhood = 1 - frac_np) %>% 
    pivot_longer(cols = all_of(c('sil_width', 
                                 'frac_misclassified', 
                                 'frac_var', 
                                 'frac_flaw_neighborhood', 
                                 'accuracy')), 
                 names_to = 'statistic', 
                 values_to = 'value') %>% 
    mutate(statistic = factor(statistic, 
                              c('sil_width', 
                                'frac_misclassified', 
                                'frac_var', 
                                'frac_flaw_neighborhood', 
                                'accuracy')), 
           dataset = factor(dataset, c('training', 'cv')))
  
  clust_dev$plot <- clust_dev$plot_data %>% 
    ggplot(aes(x = value, 
               y = reorder(algo_lab, value), 
               fill = dataset)) + 
    geom_bar(stat = 'identity', 
             color = 'black', 
             position = position_dodge(0.9)) + 
    scale_y_discrete(label = clust_labeller) + 
    scale_fill_manual(values = c(training = 'indianred3', 
                                 cv = 'steelblue'), 
                      labels = c(training = 'training', 
                                 cv = 'cross-validation'), 
                      name = 'Data set') + 
    facet_grid(. ~ statistic, 
               scales = 'free', 
               labeller = as_labeller(c(sil_width = 'mean\nsilhouette\nwidth', 
                                        frac_misclassified = 'fraction\nmisclassified', 
                                        frac_var = 'fraction\nexplained\nvariance', 
                                        frac_flaw_neighborhood = 'fraction\nmisassigned\nneighbors', 
                                        accuracy = 're-assignment\naccuracy'))) + 
    globals$common_theme + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_markdown()) + 
    labs(title = 'Clustering algorithm performance', 
         subtitle = 'TCGA training cohort', 
         x = 'Statistic value')
  
# Cachaing the results --------
  
  insert_msg('Caching the results')
  
  clust_dev <- clust_dev[c("data", "algos", "stats", "plot")]
  
  save(clust_dev, file = './cache/clust_dev.RData')
      
# END ------
  
  plan('sequential')
  
  insert_tail()