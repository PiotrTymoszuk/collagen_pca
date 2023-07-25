# Co-regulation of collagen gene expression in the tumor tissue
# Pearson's correlation

  insert_head()
  
# container -----
  
  corr <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals ------
  
  insert_msg('Globals')
  
  ## analysis table with the tumor samples only
  
  corr$analysis_tbl <- study_data %>% 
    map(~.$expression) %>% 
    map(extract_tumor_samples) %>% 
    map(~.x[globals$genes_interest$gene_symbol])
  
  ## variables in the plotting order
  
  corr$variables <- globals$genes_interest$gene_symbol %>% 
    sort
  
  ## gene pairs
  
  corr$pairs <- combn(corr$variables, 
                      m = 2, 
                      simplify = FALSE)
  
  ## n numbers
  
  corr$n_numbers <- corr$analysis_tbl %>% 
    map(nrow)
  
  corr$n_tags <- corr$n_numbers %>% 
    map(~paste('n =', .x))

# Serial correlation -----
  
  insert_msg('Correlation')
  
  corr$test_results <- corr$analysis_tbl %>% 
    future_map(function(data) corr$pairs %>% 
                 map_dfr(~correlate_variables(data, 
                                              variables = .x, 
                                              what = 'correlation', 
                                              type = 'pearson', 
                                              ci = TRUE, 
                                              pub_styled = FALSE)), 
               .options = furrr_options(seed = TRUE))
  
  ## multiple testing correction
  ## defining the significant correlations
  
  corr$test_results <- corr$test_results %>% 
    map(mutate, 
        p_adjusted = p.adjust(p_adjusted, 'BH'), 
        significant = ifelse(p_adjusted < 0.05, 'yes', 'no'), 
        regulation = ifelse(significant == 'no', 
                            'ns', 
                            ifelse(estimate > 0, 
                                   'positive', 'negative')), 
        regulation = factor(regulation, c('positive', 'negative', 'ns')), 
        pair = paste(variable1, variable2, sep = ':'))
  
  ## plot captions to be used in scatter plots
  ## and plot titles
  
  corr$test_results <- corr$test_results %>% 
    map2(., names(.), 
         ~mutate(.x, 
                 cohort = globals$study_labels[.y], 
                 plot_cap = paste0('\u03C1 = ', signif(estimate, 2), 
                                   ' [', signif(lower_ci, 2), ' - ', 
                                   signif(upper_ci, 2), ']'), 
                 plot_cap = paste(plot_cap, significance, sep = ', '), 
                 plot_cap = paste(plot_cap, n, sep = ', n = '), 
                 plot_title = paste0('<b><em>', variable1, 
                                     '</em> and <em>', variable2, '</em>, ', 
                                     cohort, '</b>')))
  
# Identification of significant correlations -------
  
  insert_msg('Significant correlations')
  
  corr$significant <- corr$test_results %>% 
    map(filter, significant == 'yes') %>% 
    map(dlply, 'regulation', function(x) x$pair) %>% 
    transpose
  
  ## correlations detected in all cohorts
  
  corr$common <- corr$significant %>% 
    map(reduce, intersect)
  
# Summary bubble-correlograms -------
  
  insert_msg('Correlograms')
  
  corr$bubble_plots <- 
    list(x = corr$test_results %>% 
           map(filter, significant == 'yes'), 
         y = globals$study_labels[names(corr$analysis_tbl)], 
         z = corr$n_tags) %>% 
    pmap(function(x, y, z) x %>% 
           ggplot(aes(x = variable1, 
                      y = variable2, 
                      fill = estimate, 
                      size = abs(estimate))) + 
           geom_point(shape = 21) + 
           scale_x_discrete(limits = corr$variables) + 
           scale_y_discrete(limits = corr$variables) + 
           scale_fill_gradient2(low = 'steelblue', 
                                mid = 'white', 
                                high = 'firebrick', 
                                midpoint = 0, 
                                limits = c(-1, 1), 
                                name = 'r') + 
           scale_radius(limits = c(0, 1), 
                        breaks = seq(0, 1, by = 0.2), 
                        name = 'abs(r)', 
                        range = c(1, 4)) + 
           globals$common_theme + 
           theme(axis.title = element_blank(), 
                 axis.text.x = element_text(hjust = 1, 
                                            angle = 90, 
                                            face = 'italic'), 
                 axis.text.y = element_text(face = 'italic')) + 
           labs(title = y, 
                subtitle = paste('Pearson correlation, significant pairs,', z)))
         
# Upset plots with the numbers of significant correlations -----
  
  insert_msg('Upset plots')
  
  corr$upset_plots <- 
    list(plotting_lst = corr$significant %>% 
           map(~set_names(.x, globals$study_labels[names(.x)])), 
         plot_title = c('Positive correlations', 
                        'Negative correlations'), 
         label_common = c(TRUE, FALSE), 
         rel_widths = list(c(2, 1), 
                           c(0.97, 0.03))) %>% 
    pmap(plot_upset, 
         plot_subtitle = 'Pairwise Pearson correlation', 
         fct_per_line = 2, 
         y_lab = '# gene pairs')
  
# Traditional scatter plots for each gene pair -------
  
  insert_msg('Gene pair scatter plots')

  for(i in names(corr$analysis_tbl)) {
    
    corr$plots[[i]] <- 
      list(variables = map2(corr$test_result[[i]]$variable1, 
                            corr$test_result[[i]]$variable2, c), 
           plot_title = corr$test_results[[i]]$plot_title, 
           plot_subtitle = corr$test_results[[i]]$plot_cap, 
           x_lab = paste0('<em>', corr$test_results[[i]]$variable1, 
                          '</em>, log<sub>2</sub> expression'), 
           y_lab = paste0('<em>', corr$test_results[[i]]$variable2, 
                          '</em>, log<sub>2</sub> expression')) %>% 
      future_pmap(plot_correlation, 
                  data = corr$analysis_tbl[[i]], 
                  type = 'correlation', 
                  point_color = globals$study_colors[[i]], 
                  cust_theme = globals$common_theme, 
                  .options = furrr_options(seed = TRUE, 
                                           packages = c('tidyverse', 
                                                        'exda'))) %>% 
      map(~.x + 
            theme(plot.tag = element_blank(), 
                  plot.title = element_markdown(), 
                  axis.title.x = element_markdown(), 
                  axis.title.y = element_markdown())) %>% 
      set_names(corr$test_results[[i]]$pair)
    
  }

# END -----
  
  plan('sequential')
  
  insert_tail()