# Comparing collagen gene expression between matched normal and tumor 
# samples with paired T test

  insert_head()
  
# container -----
  
  norm_tumor <- list()
  
# globals ------
  
  insert_msg('Globals')
  
  ## analysis tables with the matched samples
  ## arranging by the patient's ID to make the data 
  ## compatible with paired T test
  
  norm_tumor$analysis_tbl <- 
    study_data[c('GSE40272', 
                 'GSE70768', 
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
  
# serial testing ------
  
  insert_msg('Serial testing')
  
  norm_tumor$test_results <- norm_tumor$analysis_tbl %>% 
    map(test_two_groups, 
        split_fct = 'tissue', 
        variables = globals$genes_interest$gene_symbol,
        type = 't', 
        adj_method = 'BH', 
        paired = TRUE) %>% 
    map(mutate, 
        estimate = expect_x, 
        regulation = ifelse(significant == 'no', 
                            'ns', 
                            ifelse(estimate > 0, 
                                   'upregulated', 
                                   'downregulated')), 
        regulation = factor(regulation, 
                            c('upregulated', 'downregulated', 'ns')), 
        plot_lab = ifelse(significant == 'yes', response, NA), 
        plot_cap = ifelse(significant == 'yes', 
                          paste('p =', signif(p_adjusted, 2)), 
                          paste0('ns (p = ', signif(p_adjusted, 2), ')')))
  
# significantly regulated transcripts ------
  
  insert_msg('Significantly regulated transcripts')
  
  ## significantly up- and downregulated transcripts
  
  norm_tumor$significant <- norm_tumor$test_results %>% 
    map(filter, significant == 'yes') %>% 
    map(blast, regulation) %>% 
    map(map, ~.x$response) %>% 
    transpose
  
  ## transcripts regulated in all cohorts
  
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
                           fontface = 'italic') + 
           scale_color_manual(values = c(upregulated = 'firebrick', 
                                         downregulated = 'steelblue', 
                                         ns = 'gray60'), 
                              name = 'Regulation, tumor vs benign') + 
           globals$common_theme + 
           theme(plot.tag = element_blank()) + 
           labs(title = y, 
                subtitle = z, 
                x = expression('log'[2] * ' fold regulation, tumor vs benign'), 
                y = expression('-log'[10] * 'pFDR')))
  
# Venn plots for the common up- and downregulated genes -----
  
  insert_msg('Venn plots')
  
  norm_tumor$venn_plots <- norm_tumor$significant %>% 
    map(~set_names(.x, globals$study_labels[names(.x)])) %>% 
    map2(., c('Upregulated transcripts', 
              'Downregulated transcripts'), 
         ~plot_venn(plotting_lst = .x, 
                    colors = globals$study_colors[names(norm_tumor$analysis_tbl)], 
                    plot_title = .y, 
                    plot_subtitle = 'Tumor vs benign', 
                    fct_per_line = 1, 
                    rel_widths = c(0.8, 0.2), 
                    fontface = 'italic'))

# single plots ------
  
  insert_msg('Plots for single variables')
  
  for(i in names(norm_tumor$analysis_tbl)) {
    
    norm_tumor$plots[[i]] <- 
      list(variable = globals$genes_interest$gene_symbol, 
           plot_subtitle = paste(norm_tumor$n_tags, 
                                 norm_tumor$test_results[[i]]$plot_cap, 
                                 sep = ', '), 
           plot_title = paste0('<b><em>', 
                               globals$genes_interest$gene_symbol, 
                               '</em>, ', globals$study_labels[[i]], '</b>')) %>% 
      pmap(plot_variable, 
           norm_tumor$analysis_tbl[[i]], 
           split_factor = 'tissue', 
           type = 'paired', 
           y_lab = expression('log'[2] * 'expression'), 
           cust_theme = globals$common_theme) %>% 
      map(~.x + 
            scale_fill_manual(values = c(tumor = 'coral3', 
                                         benign = 'darkolivegreen3'), 
                              labels = c(tumor = 'Tumor', 
                                         benign = 'Benign'), 
                              name = '') +
            scale_x_discrete(labels = c(tumor = 'Tumor', 
                                        benign = 'Benign')) + 
            theme(plot.tag = element_blank(), 
                  plot.title = element_markdown())) %>% 
      set_names(globals$genes_interest$gene_symbol)
    
  }

# END ------
  
  rm(i)
  
  insert_tail()