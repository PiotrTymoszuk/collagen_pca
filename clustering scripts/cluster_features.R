# Comparison of collagen-related gene expression (cluster-defining factors) 
# between the collagen clusters (two-tailed T test with Cohen's d effect size). 
# FDR adjustment is applied to each cohort's results. Factors with pFDR < 0.05 
# and d >= 0.2 (at least weak effect size) are considered significant. 
# Common significant factors are shared by at least five cohorts excluding 
# GSE165060.

  insert_head()
  
# container -------
  
  clust_ft <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# expression data: normalized and identity, n numbers --------
  
  insert_msg('Expression data')
  
  clust_ft$variables <- clust_globals$variables
  
  clust_ft$assignment <- clust_semi$assignment
  
  clust_ft$data_zscores <- clust_globals$data %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map2(clust_ft$assignment, ., 
         left_join, by = 'sample_id')
    
  clust_ft$data_identity <- clust_globals$data_identity %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map2(clust_ft$assignment, ., 
         left_join, by = 'sample_id')
  
  clust_ft$n_numbers <- clust_semi$n_numbers
  
# gene classification -------
  
  insert_msg('Gene classification')
  
  clust_ft$gene_class <- 
    globals$genes_interest[c('gene_symbol', 'gene_group')] %>% 
    set_names('variable', 'gene_group') %>% 
    mutate(gene_group = stri_replace(gene_group, 
                                     fixed = ' ', 
                                     replacement = '\n'), 
           gene_group = stri_replace(gene_group, 
                                     fixed = 'adhesion', 
                                     replacement = 'Adh.'), 
           gene_group = factor(gene_group, 
                               c('Pro', 
                                 'collagen\nmodification', 
                                 'ECM\ncomponent', 
                                 'ECM\nprocessing',
                                 'Adh.')))
  
# descriptive stats -------
  
  insert_msg('Descriptive stats')
  
  clust_ft$stats <- clust_ft$data_identity %>% 
    future_map(explore, 
               variables = clust_ft$variables, 
               split_factor = 'clust_id', 
               what = 'table', 
               pub_styled = TRUE, 
               .options = furrr_options(seed = TRUE)) %>% 
    map(reduce, left_join, by = 'variable') %>% 
    map(set_names, 
        c('variable', levels(clust_ft$data_identity[[1]]$clust_id)))
  
# testing -------
  
  insert_msg('Testing')
  
  clust_ft$test <- clust_ft$data_identity %>% 
    future_map(~compare_variables(.x, 
                                  variables = clust_ft$variables, 
                                  split_factor = 'clust_id', 
                                  what = 'eff_size', 
                                  types = 'cohen_d', 
                                  ci = FALSE, 
                                  pub_styled = FALSE, 
                                  adj_method = 'BH'), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(mutate, 
        eff_size = paste(estimate_name, signif(estimate, 2), sep = ' = '), 
        plot_cap = paste(eff_size, significance, sep = ', '), 
        plot_lab = ifelse(p_adjusted < 0.05, 
                          html_bold(html_italic(variable)), 
                          html_italic(variable)))
  
# significant factors -------
  
  insert_msg('Significant factors')
  
  clust_ft$significant <- clust_ft$test %>% 
    map(filter, 
        p_adjusted < 0.05, 
        abs(estimate) >= 0.2) %>% 
    map(~.x$variable)
  
  clust_ft$common_significant <- 
    clust_ft$significant[names(clust_ft$significant) != 'gse16560'] %>% 
    shared_features(m = 5)
  
# Result table for the manuscript -------
  
  insert_msg('Result table for the manuscript')
  
  clust_ft$result_tbl <- 
    map2(clust_ft$stats, 
         map(clust_ft$test, 
             ~.x[c('variable', 'significance', 'eff_size')]), 
         left_join, by = 'variable') %>% 
    map(format_summ_tbl) %>% 
    map2(., 
         map(clust_ft$n_numbers, 
             mutate, clust_id = as.character(clust_id)), 
         ~full_rbind(tibble(variable = 'Samples, n', 
                            !!.y[[1]][1] := .y[[2]][1], 
                            !!.y[[1]][2] := .y[[2]][2]), 
                     .x))
  
# Violin plots for single variables -------
  
  insert_msg('Violin plots for single clustering factors')
  
  for(i in names(clust_ft$test)) {
    
    clust_ft$violin_plots[[i]] <- 
      list(variable = clust_ft$test[[i]]$variable, 
           plot_title = clust_ft$test[[i]]$variable %>% 
             paste0('<em>', ., '</em>') %>% 
             paste(globals$study_labels[[i]], sep = ', ') %>% 
             paste0('<b>', ., '</b>'), 
           plot_subtitle = clust_ft$test[[i]]$plot_cap) %>% 
      future_pmap(plot_variable, 
                  clust_ft$data_identity[[i]], 
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
      set_names(clust_ft$test[[i]]$variable)
    
  }
  
# Heat map limits determined by the gradient between the Collagen high and collagen low clusters -----
  
  insert_msg('Heat map limits')
  
  ## plotting order for the heat map is based on the effect size
  ## of the Collagen high vs Collagen low difference measured
  ## by Cohen's d
  
  clust_ft$hm_limits <- clust_ft$test %>% 
    map(arrange, estimate) %>% 
    map(~.$variable)
  
  ## and for the heat maps of the mean Z scores
  
  clust_ft$mean_plot_order <- clust_ft$test %>% 
    map(select, variable, estimate) %>% 
    compress(names_to = 'cohort') %>% 
    summarise(estimate = mean(estimate), .by = variable) %>% 
    arrange(estimate)
  
# Heat maps of levels of the clustering factors -------
  
  insert_msg('Heat maps of the clustering factors')
  
  clust_ft$hm_plots <- 
    list(x_object = clust_semi$clust_obj, 
         plot_title = paste('Clustering factor levels,', 
                            globals$study_labels[names(clust_semi$clust_obj)])) %>% 
    pmap(plot_clust_hm, 
         x_lab = 'Tumor sample', 
         fill_lab = 'Z-score', 
         cust_theme = globals$common_theme)
  
  ## setting the levels and highlighting of the significant effects
  
  clust_ft$hm_plots <- 
    list(x = clust_ft$hm_plots, 
         y = clust_ft$hm_limits, 
         z = clust_ft$test) %>% 
    pmap(function(x, y, z) x + 
           labs(subtitle = x$labels$tag) + 
           scale_y_discrete(limits = y, 
                            labels = function(lab) exchange(lab, z, value = 'plot_lab'))  +
           scale_fill_gradient2(low = 'steelblue', 
                                mid = 'black', 
                                high = 'firebrick', 
                                midpoint = 0, 
                                limits = c(-3, 3), 
                                oob = scales::squish, 
                                name = 'Z-score') + 
           theme(axis.text.y = element_markdown(),
                 axis.text.x = element_blank(), 
                 axis.ticks.x = element_blank(), 
                 axis.line.x = element_blank()))
  
# Ribbon plots ----
  
  insert_msg('Ribbon plots')
  
  ## significant effects highlighted in bold
  
  clust_ft$ribbon_plots <- 
    list(data = clust_ft$data_zscores, 
         plot_title = paste('Clustering factor levels,', 
                            globals$study_labels[names(clust_ft$data_zscores)]), 
         plot_subtitle = map(clust_ft$hm_plots, ~.x$labels$subtitle)) %>% 
    pmap(draw_stat_panel, 
         variables = clust_ft$variables, 
         split_factor = 'clust_id', 
         stat = 'mean', 
         err_stat = '2se', 
         form = 'line', 
         x_lab = 'mean Z score \u00B1 2\u00D7SEM', 
         cust_theme = globals$common_theme) %>% 
    map2(., clust_ft$test, 
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
  
  for(i in names(clust_ft$ribbon_plots)) {
    
    clust_ft$ribbon_plots[[i]]$data <- 
      left_join(clust_ft$ribbon_plots[[i]]$data, 
                clust_ft$gene_class, 
                by = 'variable') %>% 
      mutate(variable = factor(variable, clust_ft$hm_limits[[i]]))
    
  }
  
  clust_ft$ribbon_plots <- clust_ft$ribbon_plots %>% 
    map(~.x + 
          facet_grid(gene_group ~ ., 
                     scales = 'free', 
                     space = 'free'))
  
# Cohort-wise mean Z scores of the clustering factors -------
  
  insert_msg('Average Z scores')
  
  clust_ft$mean_data <- clust_ft$data_zscores %>% 
    map(select, -sample_id) %>% 
    map(blast, clust_id, .skip = TRUE) %>% 
    map(map, colMeans) %>% 
    map(map, 
        compress, 
        names_to = 'variable', 
        values_to = 'mean_exp') %>% 
    map(compress, 
        names_to = 'clust_id') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(clust_id = as.character(clust_id), 
           clust_id = stri_extract(clust_id, regex = 'low|hi'), 
           clust_id = factor(clust_id, c('low', 'hi')))
  
  ## plotting data: appending with the plotting order and gene classification
  
  clust_ft$mean_data <- 
    left_join(clust_ft$mean_data, 
              clust_ft$mean_plot_order, 
              by = 'variable') %>% 
    left_join(clust_ft$gene_class, by = 'variable')
  
# Heat map of the mean Z scores -------
  
  insert_msg('Mean Z score heat map')
  
  clust_ft$mean_hm <- clust_ft$mean_data %>% 
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
                         limits = c(-1, 1), 
                         name = 'Mean Z-score', 
                         oob = scales::squish) + 
    scale_x_discrete(labels = globals$study_labels, 
                     limits = names(globals$study_labels)) + 
    globals$common_theme + 
    theme(axis.title = element_blank(),
          axis.text.x = element_text(hjust = 1, 
                                     vjust = 0.5, 
                                     angle = 90), 
          axis.text.y = element_text(face = 'italic')) + 
    labs(title = 'Collagen clusters of prostate cancers', 
         subtitle = 'Mean collagen-related gene expression levels')

# END -----
  
  clust_ft$gene_class <- NULL
  clust_ft$mean_plot_order <- NULL
  clust_ft$mean_data <- NULL
  clust_ft$assignment <- NULL
  clust_ft$data_zscores <- NULL
  clust_ft$data_identity <- NULL
  clust_ft$hm_limits <- NULL
  
  clust_ft <- compact(clust_ft)

  rm(i)
  
  plan('sequential')
  
  insert_tail()