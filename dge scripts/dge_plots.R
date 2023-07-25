# Visualizations for the differentially regulated genes

  insert_head()
  
# container ------
  
  dge_plots <- list()
  
# globals ------
  
  insert_msg('Globals')
  
  ## n numbers in the clusters
  
  dge_plots$n_numbers <- dge$analysis_tbl %>% 
    map(count, clust_id)
  
  dge_plots$n_caps <- dge_plots$n_numbers %>% 
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ')) %>% 
    map_chr(paste, collapse = ', ')
  
  ## plot title prefixes
  
  dge_plots$prefixes[c("dge_collagen_int", "dge_collagen_hi")] <- 
    c('Collagen int vs low', 
      'Collagen high vs low')
  
# N number plots -------  
  
  insert_msg('Numbers of differentialy regulated genes')
  
  ## stats
  
  dge_plots$dge_numbers$regulated <- 
    dge[c("dge_collagen_int", "dge_collagen_hi")] %>% 
    map(map, count, regulation) %>% 
    map(compress, names_to = 'cohort')
  
  dge_plots$dge_numbers$total <- dge$annotation %>% 
    map_dbl(nrow) %>% 
    compress(names_to = 'cohort',
             values_to = 'n_total')
  
  dge_plots$dge_numbers$regulated <- dge_plots$dge_numbers$regulated %>% 
    map(left_join, 
        dge_plots$dge_numbers$total, 
        by = 'cohort') %>% 
    map(mutate, 
        cohort_lab = paste(globals$study_labels[cohort], 
                           n_total, sep = '\nn = '), 
        percent = n/n_total * 100, 
        cohort = factor(cohort, names(dge$analysis_tbl)))
  
  ## plots
  
  dge_plots$dge_numbers$plots <- 
    list(x = dge_plots$dge_numbers$regulated, 
         y = dge_plots$prefixes) %>% 
    pmap(function(x, y) x %>% 
           ggplot(aes(x = percent, 
                      y = cohort, 
                      fill = regulation)) + 
           geom_bar(stat = 'identity', 
                    color = 'black', 
                    position = position_dodge(0.9)) + 
           scale_y_discrete(labels = function(lab) exchange(lab, 
                                                            x, 
                                                            key = 'cohort', 
                                                            value = 'cohort_lab')) + 
           scale_fill_manual(values = c(upregulated = 'firebrick', 
                                        downregulated = 'steelblue'), 
                             name = '') + 
           globals$common_theme +
           theme(axis.title.y = element_blank()) + 
           labs(title = y, 
                x = '% of analyzed transcripts'))
  
# volcano plots -------
  
  insert_msg('Volcano plots')
  
  for(i in c('dge_collagen_int', 'dge_collagen_hi')) {
    
    dge_plots$volcano_plots[[i]] <- 
      list(data = dge[[i]] %>% 
             map(select, -counfounder), 
           plot_title = paste(dge_plots$prefixes[[i]], 
                              globals$study_labels[names(dge$test_results)], 
                              sep = ', ')) %>% 
      pmap(plot_volcano, 
           regulation_variable = 'estimate', 
           p_variable = 'p_adjusted', 
           signif_level = 0.05, 
           regulation_level = log2(1.25), 
           x_lab = expression('log'[2] * ' fold regulation'), 
           y_lab = expression('-log'[10] * 'pFDR'), 
           top_significant = 20, 
           label_variable = 'gene_symbol', 
           label_type = 'text', 
           txt_size = 2.5, 
           txt_face = 'italic', 
           fill_title = 'Regulation\nvs Collagen low', 
           cust_theme = globals$common_theme) %>% 
      map(~.x + 
            geom_vline(xintercept = 0, linetype = 'dashed') + 
            labs(subtitle = .x$labels$tag %>% 
                   stri_replace(fixed = '\n', replacement = '')) + 
            theme(plot.tag = element_blank()))
      
    
  }

# Top 20 regulated genes ------
  
  insert_msg('Top 20 regulated genes')
  
  for(i in c('dge_collagen_int', 'dge_collagen_hi')) {
    
    dge_plots$top_plots[[i]] <- 
      list(data = dge[[i]], 
           plot_title = paste(dge_plots$prefixes, 
                              globals$study_labels[names(dge$test_results)], 
                              sep = ', '), 
           plot_subtitle = paste('Collagen clusters:', 
                                 dge_plots$n_caps)) %>% 
      pmap(plot_top, 
           regulation_variable = 'estimate', 
           label_variable = 'gene_symbol', 
           p_variable = 'p_adjusted', 
           regulation_level = 0, 
           lower_ci_variable = 'lower_ci', 
           upper_ci_variable = 'upper_ci', 
           top_regulated = 20, 
           fill_title = 'Regulation\nvs Collagen low', 
           x_lab = expression('log'[2] * ' fold regulation'), 
           cust_theme = globals$common_theme) %>% 
      map(~.x + theme(axis.text.y = element_text(face = 'italic')))
    
  }

# Upset plots -------
  
  insert_msg('Upset plots')
  
  dge_plots$upset_plots <- 
    list(plotting_lst = dge[c('entrez_collagen_int', 
                              'entrez_collagen_hi')], 
         plot_title = c('Genes regulated in Collagen int vs low tumors', 
                        'Genes regulated in Collagen high vs low tumors')) %>% 
    pmap(plot_upset, 
         plot_subtitle = '', 
         feat_txt = ' ', 
         y_lab = '# regulated transcripts', 
         label_common = TRUE, 
         rel_widths = c(0.99, 0.01))
  
# END -----
  
  rm(i)
  
  insert_tail()