# Visualization of the collagen score prediction properties and calibration

  insert_head()
  
# container ------
  
  cs_plots <- list()
  
# globals ------
  
  insert_msg('Globals')
  
  ## n numbers: total and events
  
  cs_plots$n_numbers <- 
    map2(coll_score$analysis_tbl, 
         c('death', rep('relapse', 4)), 
         ~count(.x, .data[[.y]])) %>% 
    map(~c(total = sum(.x$n), 
           events = .x$n[2]))

  cs_plots$n_tags <- cs_plots$n_numbers %>% 
    map(~paste0('\ntotal: n = ', .x['total'], 
                ', events: n = ', .x['events']))
  
  ## stats table
  
  cs_plots$stats_tbl <- coll_score$stats %>% 
    mutate(type = ifelse(cohort == 'tcga', 
                         'training', 'test'), 
           type = factor(type, c('training', 'test')), 
           response = ifelse(cohort == 'GSE16560', 'OS', 'RFS'), 
           cohort_lab = globals$study_labels[cohort], 
           cohort_lab = paste(cohort_lab, response, sep = ', '), 
           cohort_lab = paste0(cohort_lab, '\nn = ', n_complete), 
           c_lab = paste0(signif(c_index, 2), 
                          ' [', signif(lower_ci, 2), 
                          ' - ', signif(upper_ci, 2), ']'), 
           r_lab = signif(raw_rsq, 2), 
           ibs_lab = signif(ibs_model, 2))
  
  ## global calibration table
  
  cs_plots$cal_tbl <- coll_score$global_cal %>% 
    mutate(plot_cap = paste0('\u03C7\u00B2(DN) = ', 
                             signif(x2_dn, 2), 
                             ', df = ', df), 
           plot_cap = ifelse(p_value < 0.05, 
                             paste0(plot_cap, ', p =', signif(p_value, 2)), 
                             paste0(plot_cap, ', ns(p = ', 
                                    signif(p_value, 2), ')')))
  
  ## tertile test results
  
  cs_plots$tertile_test <- coll_score$tertile_test
  
  ## model type colors and KM plot labels
  
  cs_plots$model_colors <-
    c(test = 'steelblue', 
      training = 'darkolivegreen4')
  
  cs_plots$KM_x_titles <- 
    c('Overall survival, months', 
      rep('Relapse-free survival, months', 4))
  
# Training model coefficient plots ------
  
  insert_msg('Training model coeffcients')
  
  cs_plots$coef_plot <- coll_score$coefs %>% 
    ggplot(aes(x = exp_coef, 
               y = reorder(variable, exp_coef), 
               color = factor(sign(coef)), 
               size = abs(coef))) + 
    geom_vline(xintercept = 1, 
               linetype = 'dashed') + 
    geom_point(shape = 16) + 
    geom_text(aes(label = signif(exp_coef, 3)), 
              size = 2.75, 
              vjust = -1.5, 
              show.legend = FALSE) + 
    scale_color_manual(values = c('-1' = 'steelblue', 
                                  '1' = 'firebrick'), 
                       labels = c('-1' = 'favorable', 
                                  '1' = 'unfavorable'), 
                       name = 'Risk factor') + 
    guides(size = 'none') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(face = 'italic')) + 
    labs(title = 'Elastic Net coefficients, TCGA training cohort', 
         subtitle = paste('RFS, Elastic Net Cox model, \u03BB =', 
                          signif(coll_score$opt_lambda$lambda[1], 2)), 
         x = expression('HR'[ElasticNet]), 
         tag = cs_plots$n_tags$tcga)
  
# Plotting of C-indexes in the training and test cohorts --------
  
  insert_msg('Plotting the C-indexes')
  
  cs_plots$c_index_plot <- cs_plots$stats_tbl %>% 
    ggplot(aes(x = c_index, 
               y = reorder(cohort_lab, c_index), 
               color = type)) + 
    geom_vline(xintercept = 0.5, 
               linetype = 'dashed') + 
    geom_errorbarh(aes(xmin = lower_ci, 
                       xmax = upper_ci), 
                   height = 0) + 
    geom_point(shape = 16, 
               size = 2) + 
    geom_text(aes(label = c_lab), 
              size = 2.75, 
              vjust = -1.3, 
              hjust = 0.3) + 
    scale_color_manual(values = cs_plots$model_colors, 
                       name = 'Cohort') + 
    scale_x_continuous(limits = c(0.42, 1), 
                     breaks = seq(0.5, 1, by = 0.1)) +
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Collagen Score, predictive performance',
         subtitle = 'Elastic Net Cox model', 
         x = 'C-index \u00B1 95%CI')
  
# Plotting the R-squares ------
  
  insert_msg('Ploting of the R-squares')
  
  cs_plots$r_sq_plot <- cs_plots$stats_tbl %>% 
    ggplot(aes(x = raw_rsq, 
               y = reorder(cohort_lab, raw_rsq), 
               fill = type)) + 
    geom_bar(stat = 'identity', 
             color = 'black') + 
    geom_text(aes(label = r_lab), 
              size = 2.75, 
              hjust = -0.5) + 
    scale_fill_manual(values = cs_plots$model_colors, 
                      name = 'Cohort') + 
    scale_x_continuous(limits = c(0, 0.55), 
                       breaks = seq(0, 0.5, by = 0.1)) + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Collagen Score, explanatory performance',
         subtitle = 'Elastic Net Cox model', 
         x = 'R\u00B2')

# Plotting C-indexes versus IBS ------
  
  insert_msg('C-index vs IBS')
  
  cs_plots$c_ibs_plot <- cs_plots$stats_tbl %>% 
    ggplot(aes(x = c_index, 
               y = 1 - ibs_model, 
               color = type)) + 
    geom_hline(yintercept = 1 - 0.5^2, 
               linetype = 'dashed') + 
    geom_vline(xintercept = 0.5, 
               linetype = 'dashed') + 
    geom_point(shape = 16, 
               size = 2) + 
    geom_text_repel(aes(label = paste(globals$study_labels[cohort], 
                                      response, sep = ', ')), 
                    size = 2.75) + 
    scale_color_manual(values = cs_plots$model_colors, 
                       name = 'Cohort') +
    globals$common_theme + 
    labs(title = 'Collagen Score, predictive performance', 
         subtitle = 'Elastic Net Cox model', 
         x = 'C-index', 
         y = '1 - IBS')
  
# Kaplan-Meier Fit plots -------
  
  insert_msg('Kaplan-Meier fit plots')
  
  cs_plots$cox_fit_plots <- 
    list(x = coll_score$models, 
         title = globals$study_labels[names(coll_score$models)], 
         xlab = cs_plots$KM_x_titles) %>% 
    pmap(plot, 
         type = 'fit', 
         cust_theme = globals$common_theme, 
         conf.int = TRUE, 
         conf.int.fill = 'black', 
         conf.int.alpha = 0.07, 
         legend = 'right') %>% 
    map(~.x$plot + 
          labs(tag = .x$lablels$tag %>% 
                 paste0('\n', .)))
  
# Calibration plots -------
  
  insert_msg('Calibration plots')
  
  ## updating the captions with the numbers of complete observations
  ## and events
  
  cs_plots$cox_cal_plots <- 
    list(x = coll_score$calibration, 
         title = globals$study_labels[names(coll_score$calibration)], 
         xlab = cs_plots$KM_x_titles) %>% 
    pmap(plot, 
         palette = c('darkolivegreen4', 
                     'steelblue4', 
                     'firebrick4'), 
         KM_size = 0.4, 
         cox_size = 0.3, 
         cust_theme = globals$common_theme, 
         legend.title = 'Collagen Score') %>% 
    map(~.x + 
          labs(subtitle = paste(stri_replace(.x$labels$tag, 
                                             regex = '^\\n', 
                                             replacement = ''), 
                                stri_replace(.x$labels$subtitle, 
                                             regex = ',\\s{1}p.*$', 
                                             replacement = ''), 
                                sep = ', ')) + 
          theme(plot.tag = element_blank()))
  
  ## adding p values for differences in survival between
  ## the Collagen score tertiles
  
  cs_plots$cox_cal_plots <- 
    map2(cs_plots$cox_cal_plots, 
         paste0('Survival difference:\n', 
                cs_plots$tertile_test$significance), 
         ~.x + 
           annotate('text', 
                    label = .y, 
                    x = 0, 
                    y = 0.1, 
                    size = 2.75, 
                    hjust = 0, 
                    vjust = 0))
  
# Plots of Brier Scores as a function of unique timepoints ------
  
  insert_msg('Brier Score plots')
  
  ## base plots
  
  cs_plots$brier_plots <- coll_score$brier_scores %>% 
    map(plot, 
        cust_theme = globals$common_theme)
  
  ## appending with titles and IBS
  
  cs_plots$brier_plots <- 
    list(x = cs_plots$brier_plots, 
         y = globals$study_labels[names(cs_plots$brier_plots)], 
         z = paste('IBS =', cs_plots$stats_tbl$ibs_lab)) %>% 
    pmap(function(x, y, z) x + 
           labs(title = y, 
                subtitle = z) + 
           scale_color_manual(values = c(reference = 'gray40', 
                                         training = 'coral3'), 
                              labels = c(reference = 'reference', 
                                         training = 'model'), 
                              name = ''))
  
# END ------
  
  insert_tail()