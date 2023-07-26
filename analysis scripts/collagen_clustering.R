# Semi-supervised clustering 
# (PAM/Manhattan algorithm, trained in the TCGA cohort)

  insert_head()
  
# container ------
  
  coll_clust <- list()
  
# Parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals: the trained clustering structure and analysis tables -----
  
  insert_msg('Globals')
  
  ## clustering variables
  
  coll_clust$variables <- clust_dev$variables
  
  ## optimal clustering structure

  coll_clust$clust_obj$tcga <- clust_dev$algos$pam_manhattan
  
  ## renaming
  
  coll_clust$clust_obj$tcga$clust_assignment <- 
    coll_clust$clust_obj$tcga$clust_assignment  %>% 
    mutate(clust_id = car::recode(clust_id, 
                                  "'1' = 'Collagen int'; 
                                  '2' = 'Collagen hi'; 
                                  '3' = 'Collagen low'"), 
           clust_id = factor(clust_id, 
                             c('Collagen hi', 
                               'Collagen int', 
                               'Collagen low')))
  
  ## tables used for the clustering structure assignment
  
  coll_clust$clust_tbl <- clust_dev$analysis_tbl
  
  ## n numbers
  
  coll_clust$n_numbers <- coll_clust$clust_tbl %>% 
    map(nrow)
  
  coll_clust$cohort_caps <- coll_clust$n_numbers %>% 
    map2_chr(., names(.), 
             ~paste0(globals$study_labels[.y], '\nn = ', .x))

# characteristic of the training clustering structure -------
  
  insert_msg('Diagnostic plots, the training cohort')
  
  ## WSS curve and silhouette statistic
  
  coll_clust$diagnostic_plots[c('wss', 
                                'silhouette')] <- 
    plot(coll_clust$clust_obj$tcga, 
         cust_theme = globals$common_theme)
  
# Semi supervised clustering -------
  
  insert_msg('Projection of the custering structures')
  
  set.seed(12345)
  
  coll_clust$semi_clust <- 
    list(newdata = coll_clust$clust_tbl[c('GSE16560', 
                                          'GSE40272', 
                                          'GSE70768', 
                                          'GSE70769')]) %>% 
    pmap(adapt_knn, 
         object = coll_clust$clust_obj$tcga, 
         type = 'propagation', 
         simple_vote = FALSE, 
         resolve_ties = TRUE, 
         max_kNN = 27)
  
  coll_clust$clust_obj[c('GSE16560', 
                          'GSE40272', 
                          'GSE70768', 
                          'GSE70769')] <- 
    coll_clust$semi_clust[c('GSE16560', 
                            'GSE40272', 
                            'GSE70768', 
                            'GSE70769')] %>% 
    map(~.x$object)
  
  coll_clust$clust_obj <- coll_clust$clust_obj[names(coll_clust$clust_tbl)]

# Clustering variances ------
  
  insert_msg('Clustering variances in the semi-supervised setting')
  
  ## clustering variance
  
  coll_clust$clust_variance <- coll_clust$clust_obj %>% 
    map(clustTools::var) %>% 
    map_dbl(~.x$frac_var) %>% 
    compress(names_to = 'cohort', 
             values_to = 'clust_variance') %>% 
    mutate(cohort_lab = coll_clust$cohort_caps[cohort], 
           type = ifelse(cohort == 'tcga', 'training', 'test'), 
           kNN = ifelse(cohort == 'tcga', 
                        NA,
                        map_dbl(coll_clust$semi_clust, ~.x$kNN)))
  
  ## bar plot with the fractions of explained variances
  
  coll_clust$variance_plot <- coll_clust$clust_variance %>% 
    ggplot(aes(x = clust_variance, 
               y = reorder(cohort_lab, clust_variance), 
               fill = type)) + 
    geom_bar(stat = 'identity', 
             color = 'black') + 
    geom_text(aes(label = signif(clust_variance, 2)), 
              size = 2.75, 
              hjust = -0.4) + 
    scale_fill_manual(values = c(test = 'steelblue', 
                                 training = 'darkolivegreen4'), 
                      name = 'Cohort') + 
    scale_x_continuous(limits = c(0, 0.53), 
                       breaks = seq(0, 0.5, by = 0.1)) +
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Clustering variance', 
         subtitle = 'PAM/Manhattan semi-supervised clustering', 
         x = 'Clustering variance')
  
# Cluster distribution ------
  
  insert_msg('Clusetr distribution')
  
  ## n numbers
  
  coll_clust$n_numbers <- coll_clust$clust_obj %>% 
    map(ngroups) %>% 
    map(arrange, desc(clust_id)) %>% 
    map(mutate,
        percent = n/sum(n) * 100, 
        y_pos = cumsum(percent) - 0.5 * percent)
  
  ## n numbers to be displayed in the plot legends
  
  coll_clust$n_legends <- coll_clust$n_numbers %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
              paste, sep = ', n = ')) %>% 
    map2(., coll_clust$n_numbers, 
         ~set_names(.x, .y[[1]]))
  
  ## stack plots
  
  coll_clust$n_plot <- coll_clust$n_numbers %>% 
    compress(names_to = 'cohort') %>% 
    ggplot(aes(x = cohort, 
               y = percent, 
               fill = clust_id)) + 
    geom_bar(position = 'stack', 
             stat = 'identity', 
             color = 'black') + 
    geom_label(aes(label = signif(percent, 2), 
                   y = y_pos), 
               size = 2.5, 
               show.legend = FALSE) + 
    scale_x_discrete(labels = coll_clust$cohort_caps) + 
    scale_fill_manual(values = globals$cluster_colors, 
                      labels = c('Collagen hi' = 'hi', 
                                 'Collagen int' = 'int', 
                                 'Collagen low' = 'low'), 
                      name = '') + 
    globals$common_theme + 
    theme(axis.title.x = element_blank()) + 
    labs(title = 'Collagen cluster distribution', 
         y = '% of cohort')
  
# Diagnostic plots for the clustering structures ------
  
  insert_msg('Diagnostic plots for the clusetring objects')
  
  ## distance heat maps
  
  coll_clust$dist_heat_maps <- coll_clust$clust_obj %>% 
    map(plot,
        type = 'heat_map', 
        cust_theme = globals$common_theme) %>% 
    map2(., globals$study_labels[names(coll_clust$clust_obj)], 
         ~.x + 
           labs(title = .y) + 
           theme(axis.text = element_blank(), 
                 axis.title = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks = element_blank()) +
           scale_fill_gradient2(low = 'firebrick',
                                mid = 'white', 
                                high = 'steelblue', 
                                midpoint = 30, 
                                limits = c(0, 60), 
                                oob = scales::squish))
  
  ## UMAP layouts
  
  coll_clust$umap_layouts <- coll_clust$clust_obj %>% 
    map(plot, 
        type = 'components', 
        with = 'data', 
        kdim = 2, 
        red_fun = 'umap', 
        cust_theme = globals$common_theme)
  
  ## MDS of the distance matrix
  
  coll_clust$dist_mds <- coll_clust$clust_obj %>% 
    map(plot, 
        type = 'components', 
        with = 'distance', 
        kdim = 2, 
        red_fun = 'mds', 
        cust_theme = globals$common_theme)
  
  ## styling
  
  for(i in c('umap_layouts', 'dist_mds')) {
    
    coll_clust[[i]] <- 
      list(x = coll_clust[[i]], 
           y = globals$study_labels[names(coll_clust$clust_obj)], 
           z = coll_clust$n_legends) %>% 
      pmap(function(x, y, z) x + 
             labs(title = y) + 
             scale_fill_manual(values = globals$cluster_colors,
                               labels = z, 
                               name = '') + 
             theme(plot.tag = element_blank()))
    
  }

# Cluster assignment tables and clustering features -------
  
  insert_msg('Cluster assignment and clustering feature tables')
  
  coll_clust$assignment <- coll_clust$clust_obj %>% 
    map(extract, 'assignment') %>% 
    map(set_names, c('patient_id', 'clust_id'))
  
  ## normalized table: for representation in plot panels
  
  coll_clust$analysis_tbl <- coll_clust$clust_tbl %>% 
    map(rownames_to_column, 'patient_id') %>% 
    map2(coll_clust$assignment, ., 
         left_join, by = 'patient_id')
  
  ## identity table: for descriptive stats and single variable plots
  
  coll_clust$identity_tbl <- study_data %>% 
    map(~.x$expression) %>% 
    map(extract_tumor_samples) %>% 
    map(~.x[c('patient_id', coll_clust$variables)]) %>% 
    map2(coll_clust$assignment, ., 
         left_join, by = 'patient_id')
  
# Testing for the differences in clustering features between the clusters -----
  
  insert_msg('Testing for the differences in clustering factors')
  
  ## descriptive stats
  
  coll_clust$stats <- coll_clust$identity_tbl %>% 
    future_map(explore, 
               variables = coll_clust$variables, 
               split_factor = 'clust_id', 
               what = 'table', 
               pub_styled = TRUE, 
               .options = furrr_options(seed = TRUE)) %>% 
    map(reduce, left_join, by = 'variable') %>% 
    map(set_names, 
        c('variable', levels(coll_clust$identity_tbl[[1]]$clust_id)))
  
  ## one-way ANOVA, appending the table with gene labels - significant
  ## effects are labeled in bold
  
  coll_clust$test <- coll_clust$identity_tbl %>% 
    future_map(~compare_variables(.x, 
                                  variables = coll_clust$variables, 
                                  split_factor = 'clust_id', 
                                  what = 'eff_size', 
                                  types = 'etasq', 
                                  ci = FALSE, 
                                  pub_styled = FALSE, 
                                  adj_method = 'BH'), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(mutate, 
        eff_size = paste('\u03B7\u00B2 =', signif(estimate, 2)), 
        plot_cap = paste(eff_size, significance, sep = ', '), 
        plot_lab = ifelse(p_adjusted < 0.05, 
                          paste0('<b><em>', variable, '</em></b>'), 
                          paste0('<em>', variable, '</em>')))
  
  ## result table
  
  coll_clust$result_tbl <- 
    map2(coll_clust$stats, 
         map(coll_clust$test, 
             ~.x[c('variable', 'significance', 'eff_size')]), 
         left_join, by = 'variable') %>% 
    map(format_summ_tbl) %>% 
    map2(., 
         map(coll_clust$n_numbers, 
             mutate, clust_id = as.character(clust_id)), 
        ~full_rbind(tibble(variable = 'Samples, n', 
                           !!.y[[1]][1] := .y[[2]][1], 
                           !!.y[[1]][2] := .y[[2]][2], 
                           !!.y[[1]][3] := .y[[2]][3]), 
                    .x))
  
  ## violin plots for single variables
  
  for(i in names(coll_clust$test)) {
    
    coll_clust$violin_plots[[i]] <- 
      list(variable = coll_clust$test[[i]]$variable, 
           plot_title = coll_clust$test[[i]]$variable %>% 
             paste0('<em>', ., '</em>') %>% 
             paste(globals$study_labels[[i]], sep = ', ') %>% 
             paste0('<b>', ., '</b>'), 
           plot_subtitle = coll_clust$test[[i]]$plot_cap) %>% 
      future_pmap(plot_variable, 
                  coll_clust$identity_tbl[[i]], 
                  split_factor = 'clust_id', 
                  type = 'violin', 
                  cust_theme = globals$common_theme, 
                  y_lab = expression('log'[2] * ' expression'), 
                  x_n_labs = TRUE, 
                  .options = furrr_options(seed = TRUE)) %>% 
      map(~.x + 
            scale_fill_manual(values = globals$cluster_colors, 
                              name = '') + 
            theme(plot.title = element_markdown())) %>% 
      set_names(coll_clust$test[[i]]$variable)
    
  }
  
# Heat map limits determined by the gradient between the Collagen high and collagen low clusters -----
  
  insert_msg('Heat map limits')
  
  ## plotting order for the heat map is based on the effect size
  ## of the Collagen high vs Collagen low difference measured
  ## by Cohen's d
  
  coll_clust$two_test <- coll_clust$analysis_tbl %>% 
    map(filter, clust_id %in% c('Collagen hi', 'Collagen low')) %>% 
    map(mutate, clust_id = droplevels(clust_id)) %>% 
    future_map(~compare_variables(.x, 
                                  variables = coll_clust$variables, 
                                  split_factor = 'clust_id', 
                                  what = 'eff_size', 
                                  types = 'cohen_d', 
                                  ci = FALSE, 
                                  pub_styled = FALSE, 
                                  adj_method = 'BH'), 
               .options = furrr_options(seed = TRUE))
  
  ## Heat map limits
  
  coll_clust$hm_limits <- coll_clust$two_test %>% 
    map(arrange, estimate) %>% 
    map(~.$variable)

# Heat maps of the clustering features ------
  
  insert_msg('Heat maps')
  
  ## clustering objects with modified cluster names
  
  coll_clust$clust_hm <- coll_clust$clust_obj
  
  for(i in names(coll_clust$clust_hm)) {
    
    coll_clust$clust_hm[[i]]$clust_assignment <- 
      coll_clust$clust_hm[[i]]$clust_assignment %>% 
      mutate(clust_id = stri_replace(clust_id, 
                                     fixed = 'Collagen ', 
                                     replacement = ''), 
             clust_id = factor(clust_id, c('low', 'int', 'hi')))
    
  }
  
  ## heat maps
  
  coll_clust$hm_plots <- 
    list(x_object = coll_clust$clust_hm, 
         plot_title = paste('Clustering factor levels,', 
                            globals$study_labels[names(coll_clust$clust_obj)])) %>% 
    pmap(plot_clust_hm, 
         x_lab = 'Tumor sample', 
         fill_lab = 'Z-score', 
         cust_theme = globals$common_theme)
  
  coll_clust$hm_plots <- 
    list(x = coll_clust$hm_plots, 
         y = coll_clust$hm_limits, 
         z = coll_clust$test) %>% 
    pmap(function(x, y, z) x + 
           scale_y_discrete(limits = y, 
                            labels = function(lab) exchange(lab, z, value = 'plot_lab'))  +
           scale_fill_gradient2(low = 'steelblue', 
                                mid = 'black', 
                                high = 'firebrick', 
                                midpoint = 0, 
                                limits = c(-3, 3), 
                                oob = scales::squish, 
                                name = 'Z-score') + 
           theme(axis.text.y = element_markdown()))

# Ribbon plots ----
  
  insert_msg('Ribbon plots')
  
  ## significant effects highlighted in bold
  
  coll_clust$ribbon_plots <- 
    list(data = coll_clust$analysis_tbl, 
         plot_title = paste('Clustering factor levels,', 
                            globals$study_labels[names(coll_clust$clust_obj)]), 
         plot_subtitle = map(coll_clust$hm_plots, ~.x$labels$tag)) %>% 
    pmap(draw_stat_panel, 
         variables = coll_clust$variables, 
         split_factor = 'clust_id', 
         stat = 'mean', 
         err_stat = '2se', 
         form = 'line', 
         x_lab = 'mean Z score \u00B1 2\u00D7SEM', 
         cust_theme = globals$common_theme) %>% 
    map2(., coll_clust$test, 
         ~.x + 
           scale_y_discrete(labels = function(x) exchange(x, .y, value = 'plot_lab')) + 
           geom_vline(xintercept = 0, 
                      linetype = 'dashed') + 
           scale_fill_manual(values = globals$cluster_colors, 
                             name = '') + 
           scale_color_manual(values = globals$cluster_colors, 
                              name = '') + 
           theme(axis.title.y = element_blank(), 
                 axis.text.y = element_markdown()))
  
  ## adding the gene function facets
  ## arranging the genes by the difference between the high and low cluster
  
  for(i in names(coll_clust$ribbon_plots)) {
    
    coll_clust$ribbon_plots[[i]]$data <- 
      left_join(coll_clust$ribbon_plots[[i]]$data, 
                globals$genes_interest[c('gene_symbol', 'gene_group')] %>% 
                  set_names('variable', 'gene_group'), 
                by = 'variable') %>% 
      mutate(variable = factor(variable, coll_clust$hm_limits[[i]]), 
             gene_group = stri_replace(gene_group, 
                                       fixed = ' ', 
                                       replacement = '\n'))
    
  }
  
  coll_clust$ribbon_plots <- coll_clust$ribbon_plots %>% 
    map(~.x + 
          facet_grid(gene_group ~ ., 
                     scales = 'free', 
                     space = 'free'))

# Heat map of normalized means ------
  
  insert_msg('Heat map of normalized means')
  
  ## variables and the plotting order
  
  coll_clust$mean_hm$plot_order <- coll_clust$two_test %>% 
    map(select, variable, estimate) %>% 
    compress(names_to = 'cohort') %>% 
    summarise(estimate = mean(estimate), .by = variable) %>% 
    arrange(estimate)
  
  ## data: mean normalized expression per clusters in the cohorts
  
  coll_clust$mean_hm$data <- coll_clust$analysis_tbl %>% 
    map(select, -patient_id) %>% 
    map(blast, clust_id) %>% 
    map(map, select, -clust_id) %>% 
    map(map, colMeans) %>% 
    map(map, 
        compress, 
        names_to = 'variable', 
        values_to = 'mean_exp') %>% 
    map(compress, 
        names_to = 'clust_id') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(clust_id = as.character(clust_id), 
           clust_id = stri_extract(clust_id, regex = 'low|int|hi'), 
           clust_id = factor(clust_id, c('low', 'int', 'hi')))
  
  ## plotting data: appending with the plotting order and gene classification
  
  coll_clust$mean_hm$data <- 
    left_join(coll_clust$mean_hm$data, 
              coll_clust$mean_hm$plot_order, 
              by = 'variable') %>% 
    left_join(globals$genes_interest[c('gene_symbol', 'gene_group')] %>% 
                set_names(c('variable', 'gene_group')), 
              by = 'variable') %>% 
    mutate(gene_group = stri_replace(gene_group, 
                                     fixed = ' ', 
                                     replacement = '\n'))
    
  
  ## heat map plot
  
  coll_clust$mean_hm$plot <- coll_clust$mean_hm$data %>% 
    ggplot(aes(x = cohort, 
               y = reorder(variable, estimate), 
               fill = mean_exp)) + 
    geom_tile() + 
    facet_grid(gene_group ~ clust_id, 
               scales = 'free', 
               space = 'free') + 
    scale_fill_gradient2(low = 'steelblue', 
                         mid = 'black', 
                         high = 'firebrick', 
                         midpoint = 0, 
                         limits = c(-1.5, 1.5), 
                         name = 'Mean Z-score', 
                         oob = scales::squish) + 
    scale_x_discrete(labels = globals$study_labels) + 
    globals$common_theme + 
    theme(axis.title = element_blank(),
          axis.text.x = element_text(hjust = 1, 
                                     vjust = 0.5, 
                                     angle = 90), 
          axis.text.y = element_text(face = 'italic')) + 
    labs(title = 'Collagen clusters of prostate cancers', 
         subtitle = 'Mean collagen pathway gene expression levels')

# END ------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()