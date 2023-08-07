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

# N number plots -------  
  
  insert_msg('Numbers of differentialy regulated genes')
  
  ## stats
  
  dge_plots$dge_numbers$regulated <- 
    dge$signif_results %>% 
    map(count, regulation) %>% 
    compress(names_to = 'cohort')

  dge_plots$dge_numbers$total <- dge$annotation %>% 
    map_dbl(nrow) %>% 
    compress(names_to = 'cohort',
             values_to = 'n_total')
  
  dge_plots$dge_numbers$regulated <- 
    left_join(dge_plots$dge_numbers$regulated, 
              dge_plots$dge_numbers$total, 
              by = 'cohort') %>% 
    mutate(cohort_lab = paste(globals$study_labels[cohort], 
                              n_total, sep = '\nn = '), 
           percent = n/n_total * 100, 
           cohort = factor(cohort, names(dge$analysis_tbl)))

  ## plot
  
  dge_plots$dge_numbers$plot <- dge_plots$dge_numbers$regulated %>% 
    ggplot(aes(x = percent, 
               y = cohort, 
               fill = regulation)) + 
    geom_bar(stat = 'identity', 
             color = 'black', 
             position = position_dodge(0.9)) + 
    scale_y_discrete(labels = function(lab) exchange(lab, 
                                                     dge_plots$dge_numbers$regulated, 
                                                     key = 'cohort', 
                                                     value = 'cohort_lab')) + 
    scale_fill_manual(values = c(upregulated = 'firebrick', 
                                 downregulated = 'steelblue'), 
                      name = '') + 
    globals$common_theme +
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Differential gene regulation, collagen hi vs low', 
         x = '% of analyzed transcripts')
  
# volcano plots -------
  
  insert_msg('Volcano plots')
  
  dge_plots$volcano_plots <- 
    list(data = dge$signif_results, 
         plot_title = globals$study_labels[names(dge$test_results)]) %>% 
    pmap(plot_volcano, 
         regulation_variable = 'estimate', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = 0, 
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

# Top 20 regulated genes ------
  
  insert_msg('Top 20 regulated genes')
  
  dge_plots$top_plots <- 
    list(data = dge$signif_results, 
         plot_title = globals$study_labels[names(dge$test_results)], 
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

# END -----

  insert_tail()