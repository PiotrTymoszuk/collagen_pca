# Plots for the results of multi-parameter survival modeling

  insert_head()
  
# container ------
  
  surv_plots <- list()
  
# Tuning process -------
  
  insert_msg('Plots of the tuning process')
  
  ## for the Cox-like models
  
  surv_plots$tuning[c('ridge', 'elnet', 'lasso')] <- 
    list(data = map(list(ridge_surv, elnet_surv, lasso_surv), 
                    ~.x$lambda_tbl), 
         best_tune = map(list(ridge_surv, elnet_surv, lasso_surv), 
                         ~.x$opt_lambda), 
         plot_title = surv_globals$algo_labels[c("ridge", "elnet", "lasso")]) %>% 
    pmap(plot_reg_cox_tuning)

  ## for SVM: 
  
  surv_plots$tuning$svm <- svm_surv$tuning$summary %>% 
    mutate(best = ifelse(c_index == max(c_index), 
                         'yes', 'no')) %>% 
    ggplot(aes(x = gamma.mu, 
               y = c_index, 
               fill = best)) + 
    geom_point(shape = 21, 
               size = 2) + 
    scale_x_continuous(trans = 'log', 
                       labels = function(x) signif(x, 2)) + 
    scale_fill_manual(values = c(no = 'steelblue', 
                                 yes = 'coral3'), 
                      name = 'best tune') + 
    globals$common_theme + 
    labs(title = surv_globals$algo_labels["svm"], 
         subtitle = paste('Best tune: \u03B3 =', 
                          signif(svm_surv$tuning$best_tune$gamma.mu, 3)), 
         x = expression(gamma),
         y = 'C-index, CV')
  
  ## for RF: 
  
  surv_plots$tuning$rf <- rf_surv$tuning$summary %>% 
    mutate(best = ifelse(c_index == max(c_index), 
                         'yes', 'no')) %>% 
    ggplot(aes(x = mtry, 
               y = c_index, 
               color = splitrule, 
               fill = best)) + 
    geom_line(aes(group = splitrule)) + 
    geom_point(shape = 21, 
               size = 1, 
               color = 'black') + 
    scale_fill_manual(values = c(no = 'steelblue', 
                                yes = 'coral3'), 
                     name = 'best tune') + 
    facet_grid(nsplit ~ nodesize) + 
    globals$common_theme + 
    labs(title = surv_globals$algo_labels["rf"], 
         subtitle = paste0('Best tune: ', 
                           'mtry = ', signif(rf_surv$tuning$best_tune$mtry[1], 2), 
                           ', splitrule: ', rf_surv$tuning$best_tune$splitrule[1], 
                           ', n splits = ', rf_surv$tuning$best_tune$nsplit[1], 
                           ', node size = ', rf_surv$tuning$best_tune$nodesize[1]), 
         x = 'mtry', 
         y = 'C-index, CV')
  
  ## for the GBM model
  
  surv_plots$tuning$gbm <- gbm_surv$tuning$summary %>% 
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
    labs(title = surv_globals$algo_labels["gbm"], 
         subtitle = paste0('Best tune: ', 
                           'N trees = ', gbm_surv$tuning$best_tune$n.trees[1], 
                           ', shrinkage = ', gbm_surv$tuning$best_tune$shrinkage[1], 
                           ', interaction depth = ', gbm_surv$tuning$best_tune$interaction.depth[1], 
                           ', minimal node size = ', gbm_surv$tuning$best_tune$n.minobsinnode[1]), 
         x = 'N trees', 
         y = 'Deviance, CV')
  
# Performance stats for the algorithms: C-index and IBS ------
  
  insert_msg('Performance stats for the algorithms: C-index and IBS')
  
  ## comparison of the cohorts
  
  surv_plots$stat_plots <- 
    list(stats = surv_summary$stats, 
         plot_subtitle = surv_globals$algo_labels[names(surv_summary$stats)] %>% 
           paste('algorithm')) %>% 
    pmap(plot_surv_stats)
  
  ## comparison of the algorithms within single cohorts
  
  surv_plots$algorithm_stat_plots <- surv_summary$stats %>% 
    compress(names_to = 'algorithm') %>% 
    mutate(cohort = factor(cohort, surv_summary$stats[[1]]$cohort)) %>% 
    blast(cohort) %>% 
    list(stats = ., 
         plot_title = surv_globals$study_labels[names(.)]) %>% 
    pmap(plot_surv_stats, 
         palette = set_names(surv_globals$algo_colors, 
                             unname(surv_globals$algo_labels)), 
         labels = surv_globals$algo_labels, 
         label_variable = 'algorithm', 
         color_variable = 'algorithm', 
         plot_subtitle = NULL) 
  
# Brier scores for the unique time points -------
  
  insert_msg('Brier score for the unique time points')
  
  for(i in names(surv_summary$brier_scores)) {
    
    surv_plots$bs_plots[[i]] <- surv_summary$brier_scores[[i]] %>% 
      map(plot, cust_theme = globals$common_theme)
    
    surv_plots$bs_plots[[i]] <- 
      list(x = surv_plots$bs_plots[[i]], 
           y = surv_globals$study_labels[names(surv_plots$bs_plots[[i]])]) %>% 
      pmap(function(x, y) x + 
             labs(title = y, 
                  subtitle = paste(surv_globals$algo_labels[[i]], 'algorithm'), 
                  x = surv_globals$algo_xlabs[[i]]) + 
             scale_color_manual(values = c(reference = 'gray60', 
                                           training = 'coral3'), 
                                labels = c(reference = 'random prediction', 
                                           training = 'model'), 
                                name = ''))
    
  }

# Survival in the score tertiles -------
  
  insert_msg('Kaplan-Meier plots for survival in the score tertiles')
  
  for(i in names(surv_summary$tertile_fits)) {
    
    surv_plots$km_plots[[i]] <- 
      list(fit = surv_summary$tertile_fits[[i]], 
           n_numbers = surv_summary$tertile_n[[i]], 
           p_value = surv_summary$tertile_test[[i]]$significance, 
           plot_title = surv_globals$study_labels[names(surv_summary$tertile_fits[[i]])]) %>% 
      pmap(plot_tertile_km, 
           x_lab = surv_globals$algo_xlabs[[i]])
    
  }
  
# Variable importance -------
  
  insert_msg('Variable importance')
  
  surv_plots$importance_plots <- 
    list(data = surv_summary$importance, 
         imp_stat = c(rep('coef', 3), 
                      rep('delta_c_index', 2), 
                      'rel_influence'), 
         labeller = list(c(positive = 'top unfavorable', 
                           negative = 'top favorable'), 
                         c(positive = 'top unfavorable', 
                           negative = 'top favorable'), 
                         c(positive = 'top unfavorable', 
                           negative = 'top favorable'), 
                         c(positive = 'top improvement', 
                           negative = 'top worsening'), 
                         c(positive = 'top improvement', 
                           negative = 'top worsening'), 
                         c(positive = 'top improvement', 
                           negative = 'top worsening')), 
         plot_subtitle = paste(surv_globals$algo_labels[names(surv_summary$importance)], 
                               'algorithm'), 
         x_lab = list(expression('log HR'[Ridge]), 
                      expression('log HR'[ElasticNet]), 
                      expression('log HR'[LASSO]), 
                      expression(Delta * ' C-index'), 
                      expression(Delta * ' C-index'), 
                      expression(Delta * ' SSE')), 
         n_top = c(rep(20, 5), 100), 
         palette = list(c(positive = 'firebrick', 
                          negative = 'steelblue', 
                          ns = 'gray70'), 
                        c(positive = 'firebrick', 
                          negative = 'steelblue', 
                          ns = 'gray70'), 
                        c(positive = 'firebrick', 
                          negative = 'steelblue', 
                          ns = 'gray70'), 
                        c(positive = 'firebrick', 
                          negative = 'steelblue', 
                          ns = 'gray70'), 
                        c(positive = 'firebrick', 
                          negative = 'steelblue', 
                          ns = 'gray70'), 
                        c(positive = 'steelblue', 
                          negative = 'steelblue', 
                          ns = 'gray70'))) %>% 
    pmap(plot_surv_importance, 
         form = 'bar') %>% 
    map(~.x + 
          geom_vline(xintercept = 0,
                     linetype = 'dashed') + 
          theme(legend.position = 'none'))
  
# Error structure: plots of square errors of prediction -------
  
  insert_msg('Plots of square errors of prediction')
  
  ## the square errors are available
  
  surv_plots$square_plots <- list(ridge = ridge_surv, 
                                  elnet = elnet_surv, 
                                  lasso = lasso_surv, 
                                  gbm = gbm_surv) %>% 
    map(~.x$calibration) %>% 
    map(map, 
        plot, 
        type = 'squares')
  
  ## plot titles and styling
  ## transposition in a more handy format
  
  for(i in names(surv_plots$square_plots)) {
    
    surv_plots$square_plots[[i]] <- 
      list(x = surv_plots$square_plots[[i]], 
           y = paste(surv_globals$algo_labels[[i]], 
                     surv_globals$study_labels[names(surv_plots$square_plots[[i]])], 
                     sep = ', ')) %>% 
      pmap(function(x, y) x %>% 
             map(~.x + 
                   labs(title = y) + 
                   globals$common_theme + 
                   theme(panel.grid.major.x = element_blank())))
    
  }
  
  surv_plots$square_plots <- surv_plots$square_plots %>% 
    map(transpose) %>%
    transpose
  
  for(i in names(surv_plots$square_plots$time)) {
    
    surv_plots$square_plots$time[[i]] <- 
      surv_plots$square_plots$time[[i]] %>% 
      map(~.x + labs(x = surv_globals$algo_xlabs[[i]]))
    
  }
  
  surv_plots$square_plots$observation <- 
    surv_plots$square_plots$observation %>% 
    map(map, 
        ~.x + 
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank()) + 
          geom_hline(yintercept = 0.25, linetype = 'dashed'))

# END -------
  
  rm(i)
  
  insert_tail()