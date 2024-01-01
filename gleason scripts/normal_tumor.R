# Comparing collagen gene expression between matched normal and tumor 
# samples with paired T test and Cohen's d. 
#
# Differentially regulated transcripts are identified by the pFDR < 0.05 and 
# Cohen's d >= 0.02 cutoffs. 

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
  
  norm_tumor$genes <- globals$genes_interest$gene_symbol
  
  norm_tumor$analysis_tbl <- 
    list(gse70768 = gse70768, 
         tcga = tcga) %>% 
    map(~.x$expression) %>% 
    map(select, 
        patient_id, tissue_type, 
        all_of(norm_tumor$genes)) %>% 
    map(~filter(.x, patient_id %in% .x$patient_id[duplicated(.x$patient_id)])) %>% 
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
                        split_factor = 'tissue_type', 
                        variables = norm_tumor$genes, 
                        what = 'table', 
                        pub_styled = TRUE), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(format_desc)
  
# serial testing ------
  
  insert_msg('Serial testing')
  
  norm_tumor$test <- norm_tumor$analysis_tbl %>% 
    future_map(~compare_variables(.x, 
                                  variables = norm_tumor$genes, 
                                  split_factor = 'tissue_type', 
                                  what = 'eff_size', 
                                  types = 'paired_cohen_d', 
                                  exact = FALSE, 
                                  ci = FALSE, 
                                  pub_styled = FALSE, 
                                  adj_method = 'BH'), 
               .options = furrr_options(seed = TRUE))
  
# Formatting: appending the results with labels displayed later in the plots -----
  
  insert_msg('Formatting the testing results')
  
  norm_tumor$test <- norm_tumor$test %>% 
    map2(., names(.), 
         ~mutate(.x,  
                 significant = ifelse(p_adjusted < 0.05, 'yes', 'no'), 
                 regulation = ifelse(significant == 'no', 
                                     'ns', 
                                     ifelse(estimate >= 0.2, 'upregulated', 
                                            ifelse(estimate <= -0.2, 
                                                   'downregulated', 'ns'))), 
                 regulation = factor(regulation, 
                                     c('upregulated', 'downregulated', 'ns')), 
                 plot_lab = ifelse(significant == 'yes', variable, NA), 
                 est_lab = paste(estimate_name, 
                                 signif(abs(estimate), 2), 
                                 sep = ' = '), 
                 plot_cap = paste('n =', n/2), 
                 plot_cap = paste(plot_cap, 
                                  est_lab, 
                                  significance, sep = ', '), 
                 cohort = globals$study_labels[.y], 
                 plot_title = html_italic(variable), 
                 plot_title = paste(plot_title, cohort, sep = ', '), 
                 plot_title = html_bold(plot_title)))

# significantly regulated transcripts ------
  
  insert_msg('Significantly regulated transcripts')
  
  ## significantly up- and downregulated transcripts
  
  norm_tumor$significant <- norm_tumor$test %>% 
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>% 
    map(blast, regulation) %>% 
    map(map, ~.x$variable) %>% 
    transpose
  
  ## transcripts regulated in both cohorts
  
  norm_tumor$common <- norm_tumor$significant %>% 
    map(reduce, intersect)

# Volcano plots -------
  
  insert_msg('Volcano plots')

  norm_tumor$volcano_plots <- 
    list(x = norm_tumor$test, 
         y = globals$study_labels[names(norm_tumor$test)], 
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
  
# Single plots ------
  
  insert_msg('Plots for single variables')
  
  for(i in names(norm_tumor$analysis_tbl)) {
    
    norm_tumor$plots[[i]] <- 
      list(variable = norm_tumor$test[[i]]$variable, 
           plot_subtitle = norm_tumor$test[[i]]$plot_cap, 
           plot_title = norm_tumor$test[[i]]$plot_title) %>% 
      pmap(plot_variable, 
           norm_tumor$analysis_tbl[[i]], 
           split_factor = 'tissue_type', 
           type = 'paired', 
           y_lab = expression('log'[2] * ' expression'), 
           cust_theme = globals$common_theme) %>% 
      map(~.x + 
            scale_fill_manual(values = c(tumor = 'coral3', 
                                         normal = 'darkolivegreen4'), 
                              labels = c(tumor = 'Tumor', 
                                         normal = 'Benign'), 
                              name = '') +
            scale_x_discrete(labels = c(tumor = 'Tumor', 
                                        normal = 'Benign')) + 
            theme(plot.tag = element_blank(), 
                  plot.title = element_markdown())) %>% 
      set_names(norm_tumor$genes)
    
  }

# Plotting meta for ribbon plots and heat maps --------
  
  insert_msg('Plotting meta for ribbon and heat map plots')
  
  ## variables and plotting order determined by the effect size metric
  
 # norm_tumor$plot_variables <- norm_tumor$common %>% 
  #  reduce(union)
  
  norm_tumor$plot_order <- norm_tumor$test %>% 
    map(filter, variable %in% norm_tumor$genes) %>% 
    map(select, variable, estimate) %>% 
    compress(names_to = 'cohort') %>% 
    summarise(estimate = mean(estimate), .by = variable) %>% 
    arrange(estimate) %>% 
    .$variable
  
  ## axis labs
  
 # norm_tumor$ax_labs <- norm_tumor$test %>% 
  #  map(mutate, 
   #     ax_lab = html_italic(variable), 
    #    significance = ifelse(stri_detect(significance, fixed = 'ns'), 
     #                         'ns', significance), 
      #  ax_lab = paste0(ax_lab, '<br>', est_lab, ', ', significance)) %>% 
    #map(~set_names(.x$ax_lab, .x$variable))
  
  norm_tumor$short_ax_labs <- norm_tumor$test %>% 
    map(mutate, 
        ax_lab = html_italic(variable), 
        ax_lab = ifelse(regulation %in% c('upregulated', 'downregulated'), 
                        html_bold(ax_lab), ax_lab)) %>% 
    map(~set_names(.x$ax_lab, .x$variable))
  
  ## plotting data: differentially regulated genes only, Z-scores
  
  norm_tumor$plot_data <- norm_tumor$analysis_tbl %>% 
    map(select, 
        tissue_type, 
        all_of(norm_tumor$genes)) %>% 
    map(map_dfc, 
        function(x) if(is.numeric(x)) scale(x)[, 1] else x)
  
# Ribbon plots, all genes ------
  
  insert_msg('Ribbon plots, all genes')
  
  ## all genes are presented as requested by the study team, 
  ## to spare place, we're highlighting significant effects in bold

  norm_tumor$ribbon_plots <- 
    list(data = norm_tumor$plot_data, 
         plot_title = globals$study_labels[names(norm_tumor$plot_data)], 
         plot_subtitle = norm_tumor$n_tags) %>% 
    pmap(draw_stat_panel, 
         variables = norm_tumor$genes, 
         split_factor = 'tissue_type', 
         stat = 'mean', 
         err_stat = '2se', 
         form = 'line', 
         cust_theme = globals$common_theme, 
         x_lab = 'Mean Z-score \u00B1 2 \u00D7 SEM')
  
  ## scale and axis adjustment
  
  norm_tumor$ribbon_plots <-  norm_tumor$ribbon_plots %>% 
    map2(., norm_tumor$short_ax_labs,
         ~.x + 
           scale_y_discrete(limits = norm_tumor$plot_order, 
                            labels = .y) + 
           scale_fill_manual(values = c(tumor = 'coral3', 
                                        normal = 'darkolivegreen4'), 
                             labels = c(tumor = 'Tumor', 
                                        normal = 'Benign'), 
                             name = '') + 
           scale_color_manual(values = c(tumor = 'coral3', 
                                         normal = 'darkolivegreen4'), 
                              labels = c(tumor = 'Tumor', 
                                         normal = 'Benign'), 
                              name = '') + 
           theme(axis.title.y = element_blank(), 
                 axis.text.y = element_markdown()))
  
# Heat map of the means, common genes -------
  
  insert_msg('Heat map of the mean expression levels')
  
  norm_tumor$hm_genes <- reduce(norm_tumor$common, union)

  norm_tumor$hm_data <- norm_tumor$plot_data %>% 
    map(select, tissue_type, all_of(norm_tumor$hm_genes)) %>% 
    map(blast, tissue_type) %>% 
    map(map, select, -tissue_type) %>% 
    map(map, colMeans) %>% 
    map(map, 
        compress, 
        names_to = 'variable', 
        values_to = 'mean_z_score') %>% 
    map(compress, 
        names_to = 'tissue_type') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(variable = factor(variable, norm_tumor$plot_order), 
           variable = droplevels(variable))
  
  ## heat map
  
  norm_tumor$hm_plot <- norm_tumor$hm_data %>% 
    ggplot(aes(y = variable, 
               x = cohort, 
               fill = mean_z_score)) + 
    facet_grid(. ~ tissue_type, 
               scales = 'free', 
               space = 'free', 
               labeller = as_labeller(c(normal = 'benign', tumor = 'tumor'))) + 
    geom_tile() + 
    scale_fill_gradient2(low = 'steelblue', 
                         mid = 'black', 
                         high = 'firebrick', 
                         midpoint = 0, 
                         name = 'Mean Z-score') + 
    scale_x_discrete(labels = globals$study_labels) + 
    globals$common_theme + 
    theme(axis.title = element_blank(), 
          axis.line = element_blank(),
          axis.text.y = element_text(face = 'italic')) + 
    labs(title = 'Differentially expressed genes, tumor vs benign', 
         subtitle = map2_chr(globals$study_labels[names(norm_tumor$n_tags)], 
                             norm_tumor$n_tags, 
                             paste, sep = ': ') %>% 
           paste(collapse = ', '))
  
# Result table for the manusucript -------
  
  insert_msg('Result table for the manuscript')
  
  ## appnding with the N numbers, there are only complete cases
  
  norm_tumor$result_tbl <- 
    map2(norm_tumor$stats, 
         map(norm_tumor$test, ~.x[c('variable', 'significance', 'est_lab')]), 
         left_join, by = 'variable') %>% 
    map(format_summ_tbl) %>% 
    map2(., norm_tumor$n_numbers, 
         ~full_rbind(tibble(variable = 'Samples, n', 
                            tumor = .y, 
                            normal = .y), 
                     .x))
  
  norm_tumor$result_tbl <- norm_tumor$result_tbl %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    select(cohort, variable, normal, tumor, significance, est_lab) %>% 
    set_names(c('Cohort', 'Variable', 'Benign', 'Tumor', 
                'Significance', 'Effect size'))
  
# END ------
  
  norm_tumor <- norm_tumor[c("n_numbers", 
                             "stats", 
                             "test", 
                             "significant", 
                             "common", 
                             "volcano_plots", 
                             "venn_plots", 
                             "plots", 
                             "ribbon_plots", 
                             "hm_plot", 
                             "result_tbl")]

  plan('sequential')
  
  rm(i)
  
  insert_tail()