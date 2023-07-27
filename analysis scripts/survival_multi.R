# Multi-variate survival modeling: 
# Collagen score corrected for Gleason score and age

  insert_head()
  
# container ------
  
  surv_multi <- list()
  
# globals -----
  
  insert_msg('Globals setup')
  
  ## analysis tables
  ## Gleason is coded as up to 5-6, 7 versus 8+
  
  surv_multi$analysis_tbl <- coll_score$score_tbl %>% 
    map(select, 
        patient_id, 
        collagen_score, 
        any_of(c('death', 'vitality_fup', 
                 'relapse', 'relapse_fup')))
  
  surv_multi$analysis_tbl <- study_data %>% 
    map(~.x$expression) %>% 
    map(extract_tumor_samples) %>% 
    map(select, 
        patient_id, 
        any_of(c('age', 'gleason'))) %>% 
    map(format_clinical) %>% 
    map(mutate, 
        gleason_factor = cut(gleason, 
                             c(-Inf, 6, 7, Inf), 
                             c('5 - 6', '7', '8+'))) %>% 
    map2(surv_multi$analysis_tbl, ., 
         left_join, 
         by = 'patient_id') %>% 
    map(~filter(.x, complete.cases(.x)))
  
  ## n numbers: total and events
  
  surv_multi$n_numbers <- 
    map2(surv_multi$analysis_tbl, 
         c('death', rep('relapse', 4)), 
         ~count(.x, .data[[.y]])) %>% 
    map(~c(total = sum(.x$n), 
           events = .x$n[2]))
  
  surv_multi$n_tags <- surv_multi$n_numbers %>% 
    map(~paste0('total: n = ', .x['total'], 
                ', events: n = ', .x['events']))
  
  ## model formulas with the continuous score
  
  surv_multi$formulas$multi <- 
    list(GSE16560 = Surv(vitality_fup, death) ~ collagen_score + age + gleason_factor, 
         GSE40272 = Surv(relapse_fup, relapse) ~ collagen_score + age + gleason_factor, 
         GSE70768 = Surv(relapse_fup, relapse) ~ collagen_score + age + gleason_factor, 
         GSE70769 = Surv(relapse_fup, relapse) ~ collagen_score + gleason_factor, 
         tcga = Surv(relapse_fup, relapse) ~ collagen_score + age + gleason_factor)
  
  surv_multi$formulas$uni <- 
    list(GSE16560 = Surv(vitality_fup, death) ~ collagen_score, 
         GSE40272 = Surv(relapse_fup, relapse) ~ collagen_score, 
         GSE70768 = Surv(relapse_fup, relapse) ~ collagen_score, 
         GSE70769 = Surv(relapse_fup, relapse) ~ collagen_score, 
         tcga = Surv(relapse_fup, relapse) ~ collagen_score)
  
  surv_multi$formulas$clinical <- 
    list(GSE16560 = Surv(vitality_fup, death) ~ age + gleason_factor, 
         GSE40272 = Surv(relapse_fup, relapse) ~ age + gleason_factor, 
         GSE70768 = Surv(relapse_fup, relapse) ~ age + gleason_factor, 
         GSE70769 = Surv(relapse_fup, relapse) ~ gleason_factor, 
         tcga = Surv(relapse_fup, relapse) ~ age + gleason_factor)
  
  ## model colors and labels
  
  surv_multi$model_labels <- 
    c(uni = 'Collagen Score', 
      clinical = 'age and Gleason', 
      multi = 'Collagen Score,\nage and Gleason')
  
  surv_multi$model_colors <- 
    c(uni = 'steelblue', 
      clinical = 'gray60', 
      multi = 'coral3')

# Construction of the models -----
  
  insert_msg('Building the models')
  
  for(i in names(surv_multi$formulas)) {
    
    surv_multi$models[[i]] <- 
      map2(surv_multi$analysis_tbl, 
           surv_multi$formulas[[i]], 
           ~call2('coxph', 
                  formula = .y, 
                  data = .x, 
                  x = TRUE)) %>% 
      map(eval) %>% 
      map2(., surv_multi$analysis_tbl, 
           as_coxex)
    
    
  }

# Model assumptions, fit stats and inference ------
  
  insert_msg('Model assumptions, fit stats and inference')
  
  for(i in names(surv_multi$models)) {
    
    surv_multi$assumptions[[i]] <- surv_multi$models[[i]] %>% 
      map(summary, 'assumptions')
    
    surv_multi$stats[[i]] <- surv_multi$models[[i]] %>% 
      map(summary, 'fit') %>% 
      compress(names_to = 'cohort') %>% 
      mutate(c_lab = paste0(signif(c_index, 2), ' [', 
                            signif(lower_ci, 2), ' - ', 
                            signif(upper_ci, 2), ']'), 
             response = ifelse(cohort == 'GSE16560', 'OS', 'RFS'),  
             cohort_lab = paste(globals$study_labels[cohort], 
                                response, 
                                sep = ', '), 
             cohort_extended = paste0(cohort_lab, '\ntotal: n = ', 
                                      n_complete, ', events: n = ', 
                                      n_events),
             cohort = factor(cohort, names(surv_multi$analysis_tbl)))
    
    surv_multi$inference[[i]] <- surv_multi$models[[i]] %>% 
      map(summary, 'inference') %>% 
      compress(names_to = 'cohort') %>% 
      mutate(hr_lab = paste0(signif(estimate, 2), 
                             ' [', signif(lower_ci, 2), ' - ', 
                             signif(upper_ci, 2), ']'), 
             cohort = factor(cohort, names(surv_multi$analysis_tbl))) %>% 
      left_join(surv_multi$stats[[i]][c('cohort', 
                                        'n_complete', 
                                        'n_events', 
                                        'response', 
                                        'cohort_lab', 
                                        'cohort_extended')], 
                by = 'cohort')
    
  }
  
# Forest plots of the C-indexes -------
  
  insert_msg('Forest plots of the C-indexes')
  
  ## presenting only the clinical factors 
  ## and the Collagen score with clinical factors
  
  surv_multi$c_index_plot <- surv_multi$stats %>% 
    compress(names_to = 'model') %>% 
    filter(model %in% c('clinical', 'multi')) %>% 
    mutate(model = factor(model, c('clinical', 'multi'))) %>% 
    ggplot(aes(x = c_index, 
               y = reorder(cohort_extended, as.numeric(cohort)), 
               color = model)) + 
    geom_vline(xintercept = 0.5, 
               linetype = 'dashed') + 
    geom_errorbarh(aes(xmin = lower_ci, 
                       xmax = upper_ci), 
                   height = 0, 
                   position = position_dodge(0.8)) + 
    geom_point(shape = 16, 
               size = 2, 
               position = position_dodge(0.8)) + 
    geom_text(aes(label = c_lab), 
              size = 2.5, 
              hjust = 0.3, 
              vjust = -1.2, 
              position = position_dodge(0.8), 
              show.legend = FALSE) + 
    scale_color_manual(labels = surv_multi$model_labels, 
                       values = surv_multi$model_colors, 
                       name = '') +
    facet_grid(cohort_extended ~ ., 
               space = 'free', 
               scales = 'free') +
    globals$common_theme +
    theme(axis.title.y = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank()) + 
    labs(title = 'Predictive performance', 
         subtitle = 'Collagen Score, age and Gleason score, Cox models', 
         x = 'C-index \u00B1 95%CI')
  
# Plotting the R-squares ------
  
  insert_msg('Plotting the R-squares')
  
  surv_multi$rsq_plot <- surv_multi$stats %>% 
    compress(names_to = 'model') %>% 
    mutate(model = factor(model, c('uni', 'clinical', 'multi'))) %>%
    ggplot(aes(x = raw_rsq, 
               y = reorder(cohort_extended, -as.numeric(cohort)), 
               fill = model)) +
    geom_bar(stat = 'identity', 
             color = 'black', 
             position = position_dodge(0.9)) + 
    geom_text(aes(label = signif(raw_rsq, 2), 
                  color = model), 
              size = 2.75, 
              hjust = -0.4, 
              position = position_dodge(0.9),
              show.legend = FALSE)  +
    scale_fill_manual(values = surv_multi$model_colors,
                      labels = surv_multi$model_labels, 
                      name = '') + 
    scale_color_manual(values = surv_multi$model_colors,
                      labels = surv_multi$model_labels, 
                      name = '') + 
    scale_x_continuous(limits = c(0, 0.57), 
                       breaks = seq(0, 0.5, by = 0.1)) + 
    globals$common_theme +
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Explanatory performance', 
         subtitle = 'Collagen Score, age and Gleason score, Cox models', 
         x = 'R\u00B2')
  
# Plotting the C-indexes and IBS as a scatter plot -------
  
  insert_msg('C_index and IBS scatter plot')
  
  surv_multi$c_ibs_plot <- surv_multi$stats %>% 
    compress(names_to = 'model') %>% 
    mutate(model = factor(model, 
                          c('uni', 'clinical', 'multi')), 
           cohort = factor(cohort, names(surv_multi$analysis_tbl))) %>% 
    blast(cohort)
  
  for(i in names(surv_multi$c_ibs_plot)) {
    
    surv_multi$c_ibs_plot[[i]] <- surv_multi$c_ibs_plot[[i]] %>% 
      ggplot(aes(x = c_index, 
                 y = 1 - ibs_model, 
                 color = model)) + 
      geom_hline(yintercept = 1 - 0.5^2, 
                 linetype = 'dashed') +
      geom_vline(xintercept = 0.5, 
                 linetype = 'dashed') +
      geom_point(shape = 16, 
                 size = 2) + 
      scale_color_manual(values = surv_multi$model_colors, 
                         labels = surv_multi$model_labels, 
                         name = '') + 
      globals$common_theme + 
      labs(title = paste('Predictive performance,', 
                         globals$study_labels[surv_multi$c_ibs_plot[[i]]$cohort]),
           subtitle = paste0(surv_multi$c_ibs_plot[[i]]$response, 
                             ', total: n = ', 
                             surv_multi$c_ibs_plot[[i]]$n_complete, 
                             ', events: n = ', 
                             surv_multi$c_ibs_plot[[i]]$n_events), 
           x = 'C-index', 
           y = '1 - IBS')
    
  }

# Plotting the hazard ratios ------
  
  insert_msg('Plotting the hazard ratios')
  
  surv_multi$hr_plot <- surv_multi$inference %>% 
    compress(names_to = 'model') %>% 
    mutate(model = factor(model, c('uni', 'clinical', 'multi'))) %>% 
    filter(variable == 'collagen_score') %>% 
    ggplot(aes(x = estimate, 
               y = reorder(cohort_extended, -as.numeric(cohort)), 
               color = model)) + 
    geom_vline(xintercept = 1, 
               linetype = 'dashed') + 
    geom_errorbarh(aes(xmin = lower_ci, 
                       xmax = upper_ci), 
                   height = 0, 
                   position = position_dodge(0.8)) + 
    geom_point(shape = 16, 
               size = 2, 
               position = position_dodge(0.8)) + 
    geom_text(aes(label = hr_lab), 
              size = 2.75, 
              hjust = 0.2, 
              vjust = -1.2, 
              position = position_dodge(0.8)) +
    scale_color_manual(values = surv_multi$model_colors, 
                       labels = surv_multi$model_labels, 
                       name = '') +
    globals$common_theme +
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Inference', 
         subtitle = 'Uni-variable and multi-variable Collagen Score Cox models', 
         x = 'HR \u00B1 95%CI, Collagen Score')
  
# END -----
  
  rm(i)
  
  insert_tail()