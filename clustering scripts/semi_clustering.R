# Semi-supervised clustering 
# (PAM/cosine algorithm, trained in the TCGA cohort)

  insert_head()
  
# container ------
  
  clust_semi <- list()
  
# Parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals: the trained clustering structure and analysis tables -----
  
  insert_msg('Globals')
  
  ## clustering variables
  
  clust_semi$variables <- clust_globals$variables
  
  ## optimal clustering structure in the training cohort

  clust_semi$clust_obj$tcga <- clust_dev$algos$pam_cosine %>% 
    dplyr::rename(c('1' = 'Collagen low', 
                    '2' = 'Collagen hi'))

# Semi supervised clustering -------
  
  insert_msg('Projection of the custering structures')
  
  set.seed(12345)
  
  ## iterative kNN prediction: choice of the best k value
  
  clust_semi$predictions <- 
    list(newdata = clust_globals$data[names(clust_globals$data) != 'tcga']) %>% 
    pmap(prediter, 
         x = clust_semi$clust_obj$tcga, 
         select_stat = 'variance', 
         max_k = 30, 
         .parallel = TRUE, 
         resolve_ties = TRUE, 
         simple_vote = FALSE)
  
  ## extraction of plots, best tunes and the final clustering structures
  
  clust_semi$tune_plots <- clust_semi$predictions %>% 
    map(plot) %>% 
    map2(., globals$study_labels[names(clust_semi$predictions)], 
         ~.x + labs(title = .y))
  
  clust_semi$best_tunes <- clust_semi$predictions %>% 
    map(extract, 'best_tune') %>% 
    compress(names_to = 'cohort')
  
  clust_semi$clust_obj <- 
    c(clust_semi$clust_obj["tcga"], 
      map(clust_semi$predictions, extract))
  
  clust_semi$clust_obj <- clust_semi$clust_obj[names(clust_globals$data)]

# Cluster assignment tables -------
  
  insert_msg('Cluster assignment frames')
  
  clust_semi$assignment <- clust_semi$clust_obj %>% 
    map(extract, 'assignment') %>% 
    map(set_names, c('sample_id', 'clust_id'))  
  
# Diagnostic plots of the training clustering structure -------
  
  insert_msg('Diagnostic plots, the training cohort')
  
  ## WSS curve and silhouette statistic
  
  clust_semi$diagnostic_plots[c('wss', 
                                'silhouette')] <- 
    plot(clust_semi$clust_obj$tcga, 
         cust_theme = globals$common_theme)

# Numeric stats of cluster performance  ------
  
  insert_msg('Clustering stats in the semi-supervised setting')
  
  clust_semi$stats <- clust_semi$clust_obj %>% 
    map(summary) %>% 
    compress(names_to = 'cohort') %>% 
    left_join(clust_semi$best_tunes, by = 'cohort') %>% 
    mutate(frac_flaw_neighborhood = 1 - frac_np, 
           dataset = ifelse(cohort == 'tcga', 
                            'training', 'test'))
  
# Bar plots with the numeric stats -------
  
  insert_msg('Bar plots of the numeric stats')
  
  clust_semi$stat_plots <- 
    list(x = c('sil_width', 'frac_misclassified', 'frac_var', 'frac_flaw_neighborhood'), 
         y = c('Cluster separation', 
               'Misclassification rate', 
               'Explained variance', 
               'Misassigned neighbors'),
         z = c('mean silhouette width', 
               'fraction negative silhouette widths', 
               'fraction explained variance', 
               '1 - neighborhood preservation rate')) %>% 
    pmap(function(x, y, z) clust_semi$stats %>% 
           ggplot(aes(x = .data[[x]], 
                      y = reorder(cohort, .data[[x]]), 
                      fill = dataset)) + 
           geom_bar(stat = 'identity', 
                    color = 'black') +
           geom_text(aes(label = signif(.data[[x]], 2)), 
                     color = 'white', 
                     hjust = 1.2, 
                     size = 2.75) + 
           scale_fill_manual(values = c(test = 'steelblue', 
                                        training = 'indianred3'), 
                             name = 'Data set') + 
           scale_y_discrete(labels = clust_globals$cohort_caps) + 
           globals$common_theme + 
           theme(axis.title.y = element_blank()) + 
           labs(title = y, 
                x = z)) %>% 
    set_names(c('sil_width', 'frac_misclassified', 'frac_var', 'frac_np'))
  
# Cluster distribution ------
  
  insert_msg('Clusetr distribution')
  
  ## n numbers
  
  clust_semi$n_numbers <- clust_semi$clust_obj %>% 
    map(ngroups) %>% 
    map(arrange, desc(clust_id)) %>% 
    map(mutate,
        percent = n/sum(n) * 100, 
        x_pos = cumsum(percent) - 0.5 * percent)
  
  ## n numbers to be displayed in the plot legends
  
  clust_semi$n_legends <- clust_semi$n_numbers %>% 
    map(mutate, 
        clust_id = stri_extract(clust_id, regex = 'low|hi')) %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = '\nn = ')) %>% 
    map2(., clust_semi$n_numbers, 
         ~set_names(.x, .y[[1]]))
  
  ## stack plots
  
  clust_semi$n_plot <- clust_semi$n_numbers %>% 
    compress(names_to = 'cohort') %>% 
    ggplot(aes(x = percent, 
               y = cohort, 
               fill = clust_id)) + 
    geom_bar(position = 'stack', 
             stat = 'identity', 
             color = 'black') + 
    geom_label(aes(label = signif(percent, 2), 
                   x = x_pos), 
               size = 2.5, 
               show.legend = FALSE) + 
    scale_y_discrete(labels = clust_globals$cohort_caps, 
                     limits = names(clust_semi$clust_obj)) + 
    scale_fill_manual(values = globals$cluster_colors, 
                      labels = c('Collagen hi' = 'hi', 
                                 'Collagen low' = 'low'), 
                      name = '') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Collagen cluster distribution', 
         x = '% of cohort')
  
# Distance diagnostic plots for the clustering structures ------
  
  insert_msg('Distance heat maps and MDS')
  
  ## distance heat maps
  
  clust_semi$dist_heat_maps <- clust_semi$clust_obj %>% 
    map(plot,
        type = 'heat_map', 
        cust_theme = globals$common_theme) %>% 
    map2(., globals$study_labels[names(clust_semi$clust_obj)], 
         ~.x + 
           labs(title = .y) + 
           theme(axis.text = element_blank(), 
                 axis.title = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks = element_blank()) +
           scale_fill_gradient2(low = 'firebrick',
                                mid = 'white', 
                                high = 'steelblue', 
                                midpoint = 1, 
                                limits = c(0, 2), 
                                oob = scales::squish,
                                name = 'cosine distance'))

  ## MDS of the distance matrix
  
  clust_semi$dist_mds_plots <- clust_semi$clust_obj %>% 
    map(plot, 
        type = 'components', 
        with = 'distance', 
        kdim = 2, 
        red_fun = 'mds', 
        cust_theme = globals$common_theme)
  
  clust_semi$dist_mds_plots <- 
    list(x = clust_semi$dist_mds, 
         y = globals$study_labels[names(clust_semi$dist_mds_plots)], 
         z = clust_semi$n_legends) %>% 
    pmap(function(x, y, z) x + 
           labs(title = y) + 
           scale_fill_manual(values = globals$cluster_colors,
                             labels = z, 
                             name = '') + 
           theme(plot.tag = element_blank()))
  
# UMAP plots ----------
  
  insert_msg('UMAP plots')
  
  ## the training UMAP layout: used for predictions
  
  clust_semi$train_umap <- clust_semi$clust_obj$tcga %>% 
    components(kdim = 2, 
               red_fun = 'umap', 
               with = 'data')
  
  ## plots and styling
  
  clust_semi$data_umap_plots <- clust_semi$clust_obj %>% 
    map(plot, 
        type = 'components', 
        kdim = 2, 
        with = 'data', 
        red_fun = 'umap',
        train_object = clust_semi$train_umap,
        cust_theme = globals$common_theme)
  
  clust_semi$data_umap_plots <- 
    list(x = clust_semi$data_umap_plots, 
         y = globals$study_labels[names(clust_semi$data_umap_plots)], 
         z = clust_semi$n_legends) %>% 
    pmap(function(x, y, z) x + 
           labs(title = y) + 
           scale_fill_manual(values = globals$cluster_colors,
                             labels = z, 
                             name = '') + 
           theme(plot.tag = element_blank()))
  
# Homologous cross-distances -------
  
  insert_msg('Homologous cross-distances')
  
  clust_semi$homo_cross_dist <- clust_semi$clust_obj %>% 
    map(cross_distance)
  
  ## mean distances between the clusters and their heat maps
  
  clust_semi$homo_cross_stats <- clust_semi$homo_cross_dist %>% 
    map(summary)
  
  clust_semi$homo_cross_plots <- clust_semi$homo_cross_dist %>% 
    map(plot, 
        type = 'mean', 
        cust_theme = globals$common_theme) %>% 
    map2(., globals$study_labels[names(clust_semi$homo_cross_dist)], 
         ~.x + 
           labs(title = .y, 
                x = .y, 
                y = .y))
  
# Heterologous cross-distances to the trainig clusters ------
  
  insert_msg('Heterologous distances')
  
  clust_semi$hetero_cross_dist <- clust_semi$clust_obj %>% 
    map(cross_distance, 
        y = clust_semi$clust_obj$tcga)
  
  ## mean distances between the clusters and their heat maps
  
  clust_semi$hetero_cross_stats <- clust_semi$hetero_cross_dist %>% 
    map(summary)
  
  clust_semi$hetero_cross_plots <- clust_semi$hetero_cross_dist %>% 
    map(plot, 
        type = 'mean', 
        cust_theme = globals$common_theme) %>% 
    map2(., globals$study_labels[names(clust_semi$homo_cross_dist)], 
         ~.x + 
           labs(title = paste(.y, globals$study_labels["tcga"], sep = ' vs '),
                x = globals$study_labels["tcga"], 
                y = .y))
  
# caching the results --------
  
  insert_msg('Caching the results')
  
  clust_semi$predictions <- NULL
  clust_semi$train_umap <- NULL
  clust_semi$hetero_cross_dist <- NULL
  clust_semi$homo_cross_dist <- NULL
  
  clust_semi <- compact(clust_semi)
  
  save(clust_semi, file = './cache/clust_semi.RData')

# END ------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()