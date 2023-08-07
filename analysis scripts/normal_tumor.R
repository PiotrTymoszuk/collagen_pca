# Comparing collagen gene expression between matched normal and tumor 
# samples with paired T test

  insert_head()
  
# container -----
  
  norm_tumor <- list()
  
# parallel backend -----
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals ------
  
  insert_msg('Globals')
  
  ## analysis tables with the matched samples
  ## arranging by the patient's ID to make the data 
  ## compatible with paired T test
  
  norm_tumor$analysis_tbl <- 
    study_data[c('GSE70768', 
                 'tcga')] %>% 
    map(~.x$expression) %>% 
    map(select, 
        patient_id, tissue, 
        all_of(globals$genes_interest$gene_symbol)) %>% 
    map(~filter(.x, 
                patient_id %in% .x$patient_id[duplicated(.x$patient_id)])) %>% 
    map(mutate, 
        tissue = factor(tissue, c('tumor', 'benign'))) %>% 
    map(arrange, patient_id)
  
  ## n numbers of the normal - tumor pairs
  
  norm_tumor$n_numbers <- norm_tumor$analysis_tbl %>% 
    map(filter, !duplicated(patient_id)) %>% 
    map(nrow)
  
  norm_tumor$n_tags <- norm_tumor$n_numbers %>% 
    map(~paste('n =', .x))
  
# descriptive stats ------
  
  insert_msg('Descriptive stats')
  
  norm_tumor$stats <- norm_tumor$analysis_tbl %>% 
    future_map(~explore(.x, 
                        split_factor = 'tissue', 
                        variables = globals$genes_interest$gene_symbol, 
                        what = 'table', 
                        pub_styled = TRUE), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(reduce, left_join, by = 'variable') %>% 
    map(set_names, c('variable', levels(norm_tumor$analysis_tbl[[1]]$tissue)))
  
# serial testing ------
  
  insert_msg('Serial testing')
  
  norm_tumor$test_results <- norm_tumor$analysis_tbl %>% 
    future_map(~compare_variables(.x, 
                                  variables = globals$genes_interest$gene_symbol, 
                                  split_factor = 'tissue', 
                                  what = 'eff_size', 
                                  types = 'paired_cohen_d', 
                                  exact = FALSE, 
                                  ci = FALSE, 
                                  pub_styled = FALSE, 
                                  adj_method = 'BH'), 
               .options = furrr_options(seed = TRUE)) %>% 
    map2(., names(.), 
         ~mutate(.x,  
                 significant = ifelse(p_adjusted < 0.05, 'yes', 'no'), 
                 regulation = ifelse(significant == 'no', 
                                     'ns', 
                                     ifelse(estimate > 0, 
                                            'upregulated', 
                                            'downregulated')), 
                 regulation = factor(regulation, 
                                     c('upregulated', 'downregulated', 'ns')), 
                 plot_lab = ifelse(significant == 'yes', variable, NA), 
                 est_lab = paste(estimate_name, signif(abs(estimate), 2), sep = ' = '), 
                 plot_cap = paste('n =', n/2), 
                 plot_cap = paste(plot_cap, est_lab, significance, sep = ', '), 
                 cohort = globals$study_labels[.y], 
                 plot_title = paste0('<b><em>', variable, 
                                     '</em>, ', cohort, '</b>')))

# significantly regulated transcripts ------
  
  insert_msg('Significantly regulated transcripts')
  
  ## significantly up- and downregulated transcripts
  
  norm_tumor$significant <- norm_tumor$test_results %>% 
    map(filter, significant == 'yes') %>% 
    map(blast, regulation) %>% 
    map(map, ~.x$variable) %>% 
    transpose
  
  ## transcripts regulated in both cohorts
  
  norm_tumor$common <- norm_tumor$significant %>% 
    map(reduce, intersect)

# volcano plots -------
  
  insert_msg('Volcano plots')
  
  norm_tumor$volcano_plots <- 
    list(x = norm_tumor$test_results, 
         y = globals$study_labels[names(norm_tumor$test_results)], 
         z = norm_tumor$n_tags) %>% 
    pmap(function(x, y, z) ggplot(x, 
                                  aes(x = estimate, 
                                      y = -log10(p_adjusted), 
                                      color = regulation)) + 
           geom_vline(xintercept = 0, 
                      linetype = 'dashed') + 
           geom_hline(yintercept = -log10(0.05), 
                      linetype = 'dashed') + 
           geom_point(size = 2, 
                      shape = 16) + 
           geom_text_repel(aes(label = plot_lab), 
                           size = 2.5, 
                           fontface = 'italic', 
                           show.legend = FALSE) + 
           scale_color_manual(values = c(upregulated = 'firebrick', 
                                         downregulated = 'steelblue', 
                                         ns = 'gray60'), 
                              name = 'Regulation, tumor vs benign') + 
           globals$common_theme + 
           theme(plot.tag = element_blank()) + 
           labs(title = y, 
                subtitle = z, 
                x = expression("Effect size, cohen's d"), 
                y = expression('-log'[10] * 'pFDR')))
  
# Venn plots for the common up- and downregulated genes -----
  
  insert_msg('Venn plots')
  
  ## displaying common regulated genes (at least two cohorts!)
  
  norm_tumor$venn_plots <- norm_tumor$significant %>% 
    map(~set_names(.x, globals$study_labels[names(.x)])) %>% 
    list(plotting_lst = ., 
         plot_title = c('Upregulated transcripts', 
                        'Downregulated transcripts'), 
         text = norm_tumor$common) %>% 
    pmap(plot_venn, 
         show_text = FALSE, 
         colors = globals$study_colors[names(norm_tumor$analysis_tbl)], 
         plot_subtitle = 'Tumor vs benign', 
         fct_per_line = 1, 
         rel_widths = c(0.8, 0.2), 
         fontface = 'italic')
  
# single plots ------
  
  insert_msg('Plots for single variables')
  
  for(i in names(norm_tumor$analysis_tbl)) {
    
    norm_tumor$plots[[i]] <- 
      list(variable = norm_tumor$test_results[[i]]$variable, 
           plot_subtitle = norm_tumor$test_results[[i]]$plot_cap, 
           plot_title = norm_tumor$test_results[[i]]$plot_title) %>% 
      pmap(plot_variable, 
           norm_tumor$analysis_tbl[[i]], 
           split_factor = 'tissue', 
           type = 'paired', 
           y_lab = expression('log'[2] * ' expression'), 
           cust_theme = globals$common_theme) %>% 
      map(~.x + 
            scale_fill_manual(values = c(tumor = 'coral3', 
                                         benign = 'darkolivegreen4'), 
                              labels = c(tumor = 'Tumor', 
                                         benign = 'Benign'), 
                              name = '') +
            scale_x_discrete(labels = c(tumor = 'Tumor', 
                                        benign = 'Benign')) + 
            theme(plot.tag = element_blank(), 
                  plot.title = element_markdown())) %>% 
      set_names(globals$genes_interest$gene_symbol)
    
  }

# Ribbon plots, common genes ------
  
  insert_msg('Ribbon plots, common genes')
  
  ## variables and plotting order determined by the effect size metric
  
  norm_tumor$ribbon_plots$variables <- norm_tumor$common %>% 
    reduce(union)
  
  norm_tumor$ribbon_plots$plot_order <- norm_tumor$test_results %>% 
    map(filter, variable %in% norm_tumor$ribbon_plots$variables) %>% 
    map(select, variable, estimate) %>% 
    compress(names_to = 'cohort') %>% 
    summarise(estimate = mean(estimate), .by = variable) %>% 
    arrange(estimate)
  
  ## axis labels: all are significant, no extra labeling!
  
  norm_tumor$ribbon_plots$ax_labs <- norm_tumor$test_results %>% 
    map(filter, variable %in% norm_tumor$ribbon_plots$variables) %>% 
    map(mutate, 
        ax_lab = paste0('<em>', variable, '</em>'), 
        significance = ifelse(stri_detect(significance, fixed = 'ns'), 
                              'ns', significance), 
        ax_lab = paste0(ax_lab, '<br>', est_lab, ', ', significance), 
        #ax_lab = ifelse(significant == 'yes', 
         #               paste0('<b>', ax_lab, '</b>'), 
          #              ax_lab)
        ) %>% 
    map(~set_names(.x$ax_lab, .x$variable))
  
  ## plotting data/Z-scores and plots
  
  norm_tumor$ribbon_plots$data <- norm_tumor$analysis_tbl %>% 
    map(select, tissue, all_of(norm_tumor$ribbon_plots$variables)) %>% 
    map(map_dfc, function(x) if(is.numeric(x)) scale(x)[, 1] else x)
  
  norm_tumor$ribbon_plots$plots <- 
    list(data = norm_tumor$ribbon_plots$data, 
         plot_title = globals$study_labels[names(norm_tumor$ribbon_plots$data)], 
         plot_subtitle = paste('n =', 
                               map_dbl(norm_tumor$ribbon_plots$data, ~nrow(.x)/2))) %>% 
    pmap(draw_stat_panel, 
         variables = norm_tumor$ribbon_plots$variables, 
         split_factor = 'tissue', 
         stat = 'mean', 
         err_stat = '2se', 
         form = 'line', 
         cust_theme = globals$common_theme, 
         x_lab = 'Mean Z-score \u00B1 2 \u00D7 SEM') %>% 
    map2(., norm_tumor$ribbon_plots$ax_labs,
         ~.x + 
           scale_y_discrete(limits = norm_tumor$ribbon_plots$plot_order$variable, 
                            labels = .y) + 
           scale_fill_manual(values = c(tumor = 'coral3', 
                                        benign = 'darkolivegreen4'), 
                             labels = c(tumor = 'Tumor', 
                                        benign = 'Benign'), 
                             name = '') + 
           scale_color_manual(values = c(tumor = 'coral3', 
                                         benign = 'darkolivegreen4'), 
                              labels = c(tumor = 'Tumor', 
                                         benign = 'Benign'), 
                              name = '') + 
           theme(axis.title.y = element_blank(), 
                 axis.text.y = element_markdown()))
  
# Heat map of the means, common genes -------
  
  insert_msg('Heat map of the mean expression levels')
  
  ## plotting data, variables and plot order
  
  norm_tumor$hm_plot$variables <- norm_tumor$ribbon_plots$variables
  
  norm_tumor$hm_plot$plot_order <- norm_tumor$ribbon_plots$plot_order
  
  norm_tumor$hm_plot$data <- norm_tumor$ribbon_plots$data %>% 
    map(blast, tissue) %>% 
    map(map, select, -tissue) %>% 
    map(map, colMeans) %>% 
    map(map, 
        compress, 
        names_to = 'variable', 
        values_to = 'mean_z_score') %>% 
    map(compress, 
        names_to = 'tissue') %>% 
    compress(names_to = 'cohort')
  
  ## heat map
  
  norm_tumor$hm_plot$plot <- norm_tumor$hm_plot$data %>% 
    ggplot(aes(x = variable, 
               y = cohort, 
               fill = mean_z_score)) + 
    facet_grid(tissue ~ .) + 
    geom_tile() + 
    scale_x_discrete(limits = norm_tumor$hm_plot$plot_order$variable,
                     position = 'top') + 
    scale_fill_gradient2(low = 'steelblue', 
                         mid = 'black', 
                         high = 'firebrick', 
                         midpoint = 0, 
                         name = 'Mean Z-score') + 
    scale_y_discrete(labels = globals$study_labels) + 
    globals$common_theme + 
    theme(axis.title = element_blank(), 
          axis.line = element_blank(),
          axis.text.x = element_text(face = 'italic', 
                                     hjust = 0, 
                                     vjust = 0, 
                                     angle = 45)) + 
    labs(title = 'Differentially expressed genes, tumor vs benign', 
         subtitle = map2_chr(globals$study_labels[names(norm_tumor$n_tags)], 
                             norm_tumor$n_tags, 
                             paste, sep = ': ') %>% 
           paste(collapse = ', '))
  
# END ------
  
  plan('sequential')
  
  rm(i)
  
  insert_tail()