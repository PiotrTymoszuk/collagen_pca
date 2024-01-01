# Co-regulation of collagen gene expression in the tumor tissue
# Euclidean distances and MDS, as well as by mean-centered PCA

  insert_head()
  
# container -----
  
  coex <- list()
  
# analysis data frames ------
  
  insert_msg('Analysis data frames')
  
  coex$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(globals$genes_interest$gene_symbol)) %>% 
    map(center_data, 'mean')
  
# Numbers of observations -------
  
  insert_msg('Numbers of observations')
  
  coex$n_numbers <- coex$expression %>% 
    map_dbl(nrow) %>% 
    paste('n =', .)
  
# MDS -----
  
  insert_msg('MDS')
  
  coex$mds_objects <- coex$expression %>% 
    map(t) %>% 
    map(reduce_data, 
        kdim = 2, 
        distance_method = 'euclidean', 
        red_fun = 'mds')
  
# Plots ---------
  
  insert_msg('Plots')
  
  coex$mds_plots <- 
    list(x = coex$mds_objects, 
         plot_subtitle =  coex$n_numbers, 
         point_color = globals$study_colors[names(coex$mds_objects)]) %>% 
    pmap(plot, 
         cust_theme = globals$common_theme)
  
  ## appending with titles and gene symbols
  
  coex$mds_plots <- 
    list(x = coex$mds_plots, 
         y = globals$study_colors[names(coex$mds_plots)]) %>% 
    pmap(function(x, y) x + 
           geom_text_repel(aes(label = observation), 
                           size = 2.5,
                           fontface = 'italic') + 
           theme(plot.tag = element_blank()) + 
           labs(title = y))
  
# PCA --------
  
  insert_msg('PCA')
  
  ## as inferred from scree plots, the first 4 dimensions
  ## capture 3/4 of the variance
  
  coex$pca_objects <- coex$expression %>% 
    map(reduce_data, 
        kdim = 6, 
        red_fun = 'pca')
  
# PCA plots --------
  
  insert_msg('PCA plots')
  
  ## scree plots of the dimensions' variances
  
  coex$pca_scree_plots <- coex$pca_objects %>% 
    map(plot, 
        type = 'scree', 
        cust_theme = globals$common_theme) %>% 
    map2(., globals$study_labels[names(coex$pca_objects)], 
         ~.x + 
           theme(plot.tag = element_blank()) + 
           labs(title = .y))
  
  ## loadings
  
  coex$pca_loadings_plots <- 
    list(x = coex$pca_objects, 
         plot_subtitle =  coex$n_numbers, 
         point_color = globals$study_colors[names(coex$mds_objects)]) %>% 
    pmap(plot, 
         type = 'loadings', 
         label_points = FALSE, 
         cust_theme = globals$common_theme)
  
  ## adding cohort-specific colors and labeling points 
  ## with gene symbols
  
  coex$pca_loadings_plots <- 
    list(x = coex$pca_loadings_plots, 
         y = globals$study_colors[names(coex$mds_objects)], 
         z = globals$study_labels[names(coex$mds_objects)]) %>% 
    pmap(function(x, y, z) x + 
           geom_text_repel(aes(label = variable), 
                           size = 2.5, 
                           fontface = 'italic', 
                           color = y) + 
           labs(title = z) + 
           theme(plot.tag = element_blank()))

  
  
# END ------
  
  coex$expression <- NULL
  coex$n_numbers <- NULL
  
  coex <- compact(coex)
  
  insert_tail()