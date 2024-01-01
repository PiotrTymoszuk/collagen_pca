# Plots for results of the differential expression analysis.
# Genes used for the cluster definition are removed from the plots!!!

  insert_head()
  
# container --------
  
  ana_dgeplots <- list()

# globals ------
  
  insert_msg('Globals')
  
  ## n numbers in the clusters
  
  ana_dgeplots$n_numbers <- ana_globals$n_numbers
  
  ana_dgeplots$n_caps <- ana_dgeplots$n_numbers %>% 
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ')) %>% 
    map_chr(paste, collapse = ', ')
  
# Percentages of regulated transcriptome --------
  
  insert_msg('Percentages of regulated transcriptome')
  
  ## percentages of regulated transcriptome
  
  ana_dgeplots$dge_stats <- ana_dge$significant %>% 
    transpose %>% 
    map(map_dbl, length) %>% 
    map(compress, 
        names_to = 'regulation', 
        values = 'n') %>% 
    map2(., ana_globals$genes, 
         ~mutate(.x, 
                 n_total = length(.y),
                 percent = n/n_total * 100))
  
  ana_dge$ax_labs <- ana_dgeplots$dge_stats %>% 
    map2(names(.), ., 
         ~paste(globals$study_labels[.x], .y$n_total[1], 
                sep = '\nn = '))
  
  ana_dgeplots$dge_stats <- ana_dgeplots$dge_stats %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = factor(cohort, names(ana_dgeplots$dge_stats)))
  
  ## stack plot with the percentages
  
  ana_dgeplots$gene_number_plot <- ana_dgeplots$dge_stats %>% 
    mutate(percent = ifelse(regulation == 'downregulated', 
                            -percent, percent)) %>% 
    ggplot(aes(x = percent, 
               y = cohort, 
               fill = regulation)) + 
    geom_bar(stat = 'identity', 
             color = 'black') + 
    geom_text(aes(label = signif(abs(percent), 2),
                  x = percent * 0.8), 
              size = 2.5, 
              color = 'white') + 
    scale_x_continuous(labels = function(x) abs(x)) + 
    scale_y_discrete(labels = ana_dge$ax_labs) + 
    scale_fill_manual(values = c(upregulated = 'firebrick', 
                                 downregulated = 'steelblue'), 
                      name = '') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Regulated transcriptome, collagen hi vs low', 
         x = '% of analyzed genes')
  
# Volcano plots -------
  
  insert_msg('Volcano plots')
  
  ana_dgeplots$volcano_plots <- 
    map2(ana_dge$test, 
         map(transpose(ana_dge$significant), reduce, union), 
         ~filter(.x, response %in% .y)) %>% 
    map(filter, !response %in% globals$genes_interest$gene_symbol) %>% 
    list(data = ., 
         plot_title = globals$study_labels[names(ana_dge$test)]) %>% 
    pmap(plot_volcano, 
         regulation_variable = 'estimate', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = 0, 
         x_lab = expression('log'[2] * ' fold regulation'), 
         y_lab = expression('-log'[10] * 'pFDR'), 
         top_regulated = 20, 
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
  
  ana_dgeplots$top_plots <- 
    map2(ana_dge$test, 
         map(transpose(ana_dge$significant), reduce, union), 
         ~filter(.x, response %in% .y)) %>% 
    map(filter, !response %in% globals$genes_interest$gene_symbol) %>% 
    list(data = ., 
         plot_title = globals$study_labels[names(ana_dge$test)], 
         plot_subtitle = paste('Collagen clusters:', 
                               ana_dgeplots$n_caps)) %>% 
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