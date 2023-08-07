# Uni-variable Cox modeling of relapse-free survival and overall survival 
# as a function of collagen gene expression and of the Collagen Score
# Z scores of the explanatory variables are used for modeling

  insert_head()
  
# container -----
  
  surv_cox <- list()
  
# globals -----
  
  insert_msg('Globals')
  
  ## variables and their labels
  
  surv_cox$variables <- c(globals$genes_interest$gene_symbol, 
                          'collagen_score')
  
  surv_cox$var_labs <- ifelse(surv_cox$variables == 'collagen_score', 
                              'Collagen Score', 
                              surv_cox$variables) %>% 
    set_names(surv_cox$variables) %>% 
    compress(names_to = 'variable', 
             values_to = 'label')
  
  surv_cox$var_labs <- surv_cox$var_labs %>% 
    mutate(label = ifelse(variable != 'collagen_score', 
                          paste0('<em>', label, '</em>'), 
                          label))
  
  ## analysis tables
  ## Z score of the genes and of the Collagen Score
  ## cases with the complete survival information
  
  surv_cox$analysis_tbl <- study_data %>% 
    map(~.x$expression) %>% 
    map(extract_tumor_samples) %>% 
    map(select, 
        patient_id, 
        any_of(c('death', 'vitality_fup')), 
        any_of(c('relapse', 'relapse_fup')), 
        all_of(globals$genes_interest$gene_symbol))
  
  surv_cox$analysis_tbl <- 
    map2(surv_cox$analysis_tbl, 
         map(coll_score$score_tbl, ~.x[c('patient_id', 'collagen_score')]), 
         right_join, by = 'patient_id')
  
  surv_cox$analysis_tbl <- 
    map2(surv_cox$analysis_tbl, 
         list(c('death', 'vitality_fup'), 
              c('relapse', 'relapse_fup'), 
              c('relapse', 'relapse_fup'), 
              c('relapse', 'relapse_fup'), 
              c('relapse', 'relapse_fup')), 
         ~filter(.x, 
                 !is.na(.data[[.y[1]]]), 
                 !is.na(.data[[.y[2]]])))
  
  for(i in names(surv_cox$analysis_tbl)) {
    
    surv_cox$analysis_tbl[[i]][surv_cox$variables] <- 
      surv_cox$analysis_tbl[[i]][surv_cox$variables] %>% 
      map_dfc(~scale(.x)[, 1])
    
  }

  ## n numbers: total and events
  
  surv_cox$n_numbers <- 
    map2(surv_cox$analysis_tbl, 
         c('death', rep('relapse', 4)), 
         ~count(.x, .data[[.y]])) %>% 
    map(~c(total = sum(.x$n), 
           events = .x$n[2]))

  surv_cox$n_tags <- surv_cox$n_numbers %>% 
    map(~paste0('total: n = ', .x['total'], 
                ', events: n = ', .x['events']))
  
  ## model formulas
  
  surv_cox$formulas <- 
    c(GSE16560 = 'Surv(vitality_fup, death)', 
      GSE70768 = 'Surv(relapse_fup, relapse)', 
      GSE70769 = 'Surv(relapse_fup, relapse)',
      GSE116918 = 'Surv(relapse_fup, relapse)', 
      tcga = 'Surv(relapse_fup, relapse)') %>% 
    map(function(response) surv_cox$variables %>% 
          map(~paste(response, .x, sep = ' ~ ')) %>% 
          map(as.formula) %>% 
          set_names(surv_cox$variables))

# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# Building the models ------
  
  insert_msg('Building the models')
  
  ## GSE16560: overall survival
  ## the remaining scores: relapse-free survival
  
  surv_cox$models <- 
    list(form_lst = surv_cox$formulas, 
         data_lst = surv_cox$analysis_tbl) %>% 
    pmap(function(form_lst, data_lst) form_lst %>% 
           map(~call('coxph', 
                     formula = .x, 
                     data = data_lst, x = TRUE)) %>% 
           map(eval) %>% 
           map(~as_coxex(.x, data = data_lst)))
  
# Checking the model assumptions, fit stats and inference -------
  
  insert_msg('Model assumptions, fit stats and inference')
  
  for(i in c('assumptions', 'fit', 'inference')) {
    
    surv_cox[[i]] <- surv_cox$models %>% 
      map(future_map, 
          summary.coxex, 
          type = i, 
          .options = furrr_options(seed = TRUE)) %>% 
      map(compress, names_to = 'variable')

  }
  
  surv_cox$stats <- surv_cox$fit
  
  surv_cox$fit <- NULL
  surv_cox <- compact(surv_cox)

  surv_cox$inference <- surv_cox$inference %>% 
    map(re_adjust) %>% 
    map(mutate, 
        significant = ifelse(p_adjusted < 0.05, 'yes', 'no'), 
        correlation = ifelse(significant == 'no', 
                             'ns', 
                             ifelse(estimate > 0, 'unfavorable', 'favorable')), 
        correlation = factor(correlation, c('unfavorable', 'favorable', 'ns')), 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
# Identification of significant favorable and unfavorable markers ------
  
  insert_msg('Significant markers')
  
  surv_cox$significant <- surv_cox$inference %>% 
    map(filter, significant == 'yes') %>% 
    map(blast, correlation) %>% 
    map(map, ~.x$variable)

# Plotting the C-indexes -------
  
  insert_msg('Plotting of the C-indexes')
  
  surv_cox$c_index_plots <- 
    list(x = surv_cox$stats, 
         y = globals$study_labels[names(surv_cox$stats)], 
         v = globals$study_colors[names(surv_cox$stats)], 
         resp = c('OS', rep('RFS', 4)), 
         z = surv_cox$n_tags) %>% 
    pmap(function(x, y, v, resp, z) x %>% 
           ggplot(aes(x = c_index, 
                      y = reorder(variable, c_index))) +
           geom_vline(xintercept = 0.5, 
                      linetype = 'dashed') + 
           geom_errorbarh(aes(xmin = lower_ci, 
                              xmax = upper_ci), 
                          height = 0, 
                          color = v) + 
           geom_point(shape = 16, 
                      size = 2, 
                      color = v) + 
           scale_y_discrete(labels = function(x) exchange(x, surv_cox$var_labs)) + 
           globals$common_theme + 
           theme(axis.title.y = element_blank(), 
                 axis.text.y = element_markdown()) + 
           labs(title = paste('Predictive performance,', y), 
                subtitle = paste0(resp, ', univariable Cox modeling, ', z), 
                x = paste('C-index \u00B1 95%CI')))
  
# Plotting the normalized hazard ratios ------  
  
  insert_msg('Forest plots of the hazard ratios')
  
  surv_cox$hr_plots <- 
    list(x = surv_cox$inference, 
         y = globals$study_labels[names(surv_cox$inference)], 
         resp = c('OS', rep('RFS', 4)), 
         z = surv_cox$n_tags) %>% 
    pmap(function(x, y, resp, z) ggplot(x, 
                                        aes(x = estimate, 
                                            y = reorder(variable, estimate), 
                                            color = correlation)) +
           geom_vline(xintercept = 1, 
                      linetype = 'dashed') + 
           geom_errorbarh(aes(xmin = lower_ci, 
                              xmax = upper_ci), 
                          height = 0) + 
           geom_point(shape = 16, 
                      size = 2) + 
           scale_y_discrete(labels = function(x) exchange(x, surv_cox$var_labs)) + 
           scale_color_manual(values = c(unfavorable = 'firebrick', 
                                         favorable = 'steelblue', 
                                         ns = 'gray60'), 
                              name = 'Risk marker') + 
           globals$common_theme + 
           theme(axis.title.y = element_blank(), 
                 axis.text.y = element_markdown()) + 
           labs(title = paste('Inference,', y), 
                subtitle = paste0(resp, ', univariable Cox modeling, ', z), 
                x = paste('HR \u00B1 95%CI, normalized')))

# END -----
  
  plan('sequential')
  
  rm(i)
  
  insert_tail()