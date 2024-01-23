# Plots for the GBM survival models employing collagen expression 
# and clinical markers

  insert_head()
  
# container -------
  
  surv_mplots <- list()
  
# Tuning process -------
  
  insert_msg('Tuning process: combined clinical/collagen model')
  
  surv_mplots$tuning <- surv_multi$tuning$summary %>% 
    mutate(best = ifelse(cv_deviance == min(cv_deviance), 
                         'yes', 'no')) %>% 
    ggplot(aes(x = shrinkage, 
               y = cv_deviance, 
               color = factor(n.trees), 
               fill = best)) + 
    geom_line(aes(group = factor(n.trees))) + 
    geom_point(shape = 21, 
               size = 1, 
               color = 'black') + 
    scale_fill_manual(values = c(no = 'steelblue', 
                                 yes = 'coral3'), 
                      name = 'best tune') + 
    facet_grid(interaction.depth ~ n.minobsinnode) + 
    globals$common_theme + 
    labs(title = 'Tuning the collagen/clinical GBM model', 
         subtitle = paste0('Best tune: ', 
                           'N trees = ', surv_multi$tuning$best_tune$n.trees[1], 
                           ', shrinkage = ', surv_multi$tuning$best_tune$shrinkage[1], 
                           ', interaction depth = ', surv_multi$tuning$best_tune$interaction.depth[1], 
                           ', minimal node size = ', surv_multi$tuning$best_tune$n.minobsinnode[1]), 
         x = 'N trees', 
         y = 'Deviance, CV')
  
# Performance stats for the algorithms: C-index and IBS ------
  
  insert_msg('Performance stats for the algorithms: C-index and IBS')
  
  surv_mplots$stat_plots <- 
    list(stats = surv_multi$stats %>% 
           blast(cohort), 
         plot_title = surv_globals$study_labels[levels(surv_multi$stats$cohort)])  %>% 
    pmap(plot_surv_stats, 
         palette = set_names(surv_globals$type_colors, 
                             unname(surv_globals$type_labels)), 
         labels = surv_globals$type_labels, 
         label_variable = 'model', 
         txt_size = 2.75, 
         color_variable = 'model', 
         plot_subtitle = 'GBM algoritm')
  
# Brier scores for the unique time points -------
  
  insert_msg('Brier score for the unique time points')
  
  for(i in names(surv_multi$brier_scores)) {
    
    surv_mplots$bs_plots[[i]] <- surv_multi$brier_scores[[i]] %>% 
      map(plot, cust_theme = globals$common_theme)
    
    surv_mplots$bs_plots[[i]] <- 
      list(x = surv_mplots$bs_plots[[i]], 
           y = names(surv_mplots$bs_plots[[i]])) %>% 
      pmap(function(x, y) x + 
             labs(title = paste(surv_globals$type_labels[[y]], 
                                surv_globals$study_labels[[i]],
                                sep = ', '), 
                  x = 'Min/max scaled relapse-free survival') + 
             scale_color_manual(values = c(reference = 'gray60', 
                                           training = 'coral3'), 
                                labels = c(reference = 'random prediction', 
                                           training = 'model'), 
                                name = ''))
    
  }
  
# GBM score tertiles ------
  
  insert_msg('GBM score tertiles')
  
  ## tertile data and surv fits: used lated for KM plots
  
  surv_mplots$tertile_data <- surv_multi$score_tbl
  
  for(i in names(surv_mplots$tertile_data)) {
    
    surv_mplots$tertile_data[[i]][c('clinic_score', 'full_score', 'collagen_score')] <- 
      surv_mplots$tertile_data[[i]][c('clinic_score', 'full_score', 'collagen_score')] %>% 
      map_dfc(cut_tertiles)
    
    surv_mplots$tertile_fits[[i]] <- 
      list(clinic = Surv(rfs_months, relapse) ~ clinic_score,
           collagen = Surv(rfs_months, relapse) ~ collagen_score, 
           full = Surv(rfs_months, relapse) ~ full_score) %>% 
      map(survminer::surv_fit, data = surv_mplots$tertile_data[[i]])
    
  }
  
  ## N numbers and numbers of events in the tertiles
  
  surv_mplots$tertile_n <- surv_mplots$tertile_data %>% 
    map(function(data) c(clinic = 'clinic_score', 
                         collagen = 'collagen_score', 
                         full = 'full_score') %>% 
          map(~count(data, .data[[.x]]))) %>% 
    map(map, set_names, c('score_cuts', 'n_total'))
  
  surv_mplots$tertile_events <- surv_mplots$tertile_data %>% 
    map(filter, relapse == 1) %>% 
    map(function(data) c(clinic = 'clinic_score', 
                         collagen = 'collagen_score', 
                         full = 'full_score') %>% 
          map(~count(data, .data[[.x]]))) %>% 
    map(map, set_names, c('score_cuts', 'n_events'))
  
  surv_mplots$tertile_n <- 
    map2(surv_mplots$tertile_n, 
         surv_mplots$tertile_events, 
         function(x, y) map2(x, y, left_join, by = 'score_cuts'))
  
  ## testing results
  
  surv_mplots$tertile_test <- surv_multi$tertile_test %>% 
    blast(cohort)

# Survival in the score tertiles -------
  
  insert_msg('Kaplan-Meier plots for survival in the score tertiles')

  for(i in names(surv_mplots$tertile_fits)) {
    
    surv_mplots$km_plots[[i]] <- 
      list(fit = surv_mplots$tertile_fits[[i]], 
           n_numbers = surv_mplots$tertile_n[[i]], 
           p_value = surv_mplots$tertile_test[[i]]$significance, 
           plot_title = paste(surv_globals$type_labels[names(surv_mplots$tertile_fits[[i]])], 
                              surv_globals$study_labels[[i]], 
                              sep = ', ')) %>% 
      pmap(plot_tertile_km, 
           x_lab = surv_globals$algo_xlabs["gbm"])
    
  }
  
# Variable importance: the clinical and full models --------
  
  insert_msg('Variable importance plots')
  
  surv_mplots$importance_plots <- 
    list(data = surv_multi$importance, 
         plot_title = paste(surv_globals$type_labels[names(surv_multi$importance)])) %>% 
    pmap(plot_surv_importance, 
         plot_subtitle = 'GBM, pooled GEO cohort', 
         imp_stat = 'rel_influence', 
         form = 'bar', 
         palette = c(positive = 'steelblue', 
                     negative = 'steelblue', 
                     ns = 'gray70'), 
         x_lab = expression(Delta * ' SSE'), 
         labeller = c(positive = 'top improvement', 
                      negative = 'top worsening')) %>% 
    map(~.x + 
          scale_y_discrete(labels = function(x) ifelse(x %in% surv_multi$gene_variables, 
                                                       html_italic(x), 
                                                       exchange(x, 
                                                                globals$clinical_lexicon))) + 
          theme(axis.text.y = element_markdown(face = 'plain'), 
                strip.background = element_blank(), 
                strip.text = element_blank(), 
                legend.position = 'none'))
  
# Plots of square errors: error structure --------
  
  insert_msg('Plots of square errors')
  
  surv_mplots$square_plots <- surv_multi$calibration %>% 
    map(map, 
        plot, 
        type = 'squares')
  
  ## plot titles and styling
  ## transposition in a more handy format
  
  for(i in names(surv_mplots$square_plots)) {
    
    surv_mplots$square_plots[[i]] <- 
      list(x = surv_mplots$square_plots[[i]], 
           y = paste(surv_globals$type_labels[names(surv_mplots$square_plots[[i]])], 
                     surv_globals$study_labels[[i]], 
                     sep = ', ')) %>% 
      pmap(function(x, y) x %>% 
             map(~.x + 
                   labs(title = y) + 
                   globals$common_theme + 
                   theme(panel.grid.major.x = element_blank())))
    
  }
  
  surv_mplots$square_plots <- surv_mplots$square_plots %>% 
    map(transpose) %>%
    transpose
  
  surv_mplots$square_plots$time <- 
    surv_mplots$square_plots$time %>% 
    map(map, ~.x + labs(x = surv_globals$algo_xlabs["gbm"]))
  
  surv_mplots$square_plots$observation <- 
    surv_mplots$square_plots$observation %>% 
    map(map, 
        ~.x + 
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank()) + 
          geom_hline(yintercept = 0.25, linetype = 'dashed'))
  
# END ------
  
  rm(i)
  
  insert_tail()