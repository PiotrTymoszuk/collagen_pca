# Multi-variate survival modeling: 
# Collagen score corrected for Gleason score, tumor stage and PSA

  insert_head()
  
# container ------
  
  cs_multi <- list()
  
# modeling data -----
  
  insert_msg('Modeling data')
  
  ## analysis tables
  ## Gleason is coded as up to 5-6, 7 versus 8+
  
  cs_multi$survival <- 
    coll_score$score_tbl[c("gse54460", "gse70768", 
                           "gse70769", "gse220095", 
                           "tcga", "dkfz")] %>% 
    map(select, 
        sample_id, 
        collagen_score, 
        relapse, 
        rfs_months)
  
  cs_multi$clinic <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor')
  
  cs_multi$clinic <- cs_multi$clinic[names(cs_multi$survival)] %>% 
    map(select, 
        sample_id, 
        gleason_simple, 
        psa_diagnosis, 
        pt_stage)
  
  cs_multi$data <- 
    map2(cs_multi$survival, 
         cs_multi$clinic, 
         left_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x)))
  
  ## transformation of PSA with sqrt to improve normality
  # merging stages T3 and T4 together
  ## Z-scores of numeric variables
  
  cs_multi$data <- cs_multi$data %>% 
    map(mutate, 
        psa_diagnosis = sqrt(psa_diagnosis), 
        pt_stage = car::recode(pt_stage, 
                               "'T1' = 'T1/T2'; 'T2' = 'T1/T2'; 
                               'T3' = 'T3/T4'; 'T4' = 'T3/T4'"), 
        pt_stage = factor(pt_stage, c('T1/T2', 'T3/T4')))

  cs_multi$data <- cs_multi$data %>% 
    map(function(x) if(is.numeric(x)) scale(x)[, 1] else x)
    
# N numbers -------
  
  insert_msg('N numbers')
  
  cs_multi$n_numbers <- cs_multi$data %>% 
    map(~tibble(total = nrow(.x), 
                events = sum(.x$relapse))) %>% 
    compress(names_to = 'cohort')
  
  cs_multi$n_tags <- 
    map2(cs_multi$n_numbers$total, 
         cs_multi$n_numbers$events, 
         ~paste0('total: n = ', .x, 
                 ', events: n = ', .y)) %>% 
    set_names(cs_multi$n_numbers$cohort)
  
# Model formulas, labels and colors --------
  
  insert_msg('Model formulas, labels and colors')
  
  cs_multi$formulas <- 
    list(uni = Surv(rfs_months, relapse) ~ collagen_score, 
         clinical = Surv(rfs_months, relapse) ~ psa_diagnosis + pt_stage + gleason_simple, 
         multi = Surv(rfs_months, relapse) ~ psa_diagnosis + pt_stage + gleason_simple + collagen_score)

  ## model colors and labels
  
  cs_multi$model_labels <- 
    c(uni = 'Collagen Score', 
      clinical = 'Cilinical factors', 
      multi = 'Collagen Score and clinical factors')
  
  cs_multi$model_colors <- 
    c(uni = 'steelblue', 
      clinical = 'gray60', 
      multi = 'coral3')

# Construction of the models -----
  
  insert_msg('Building the models')

  for(i in names(cs_multi$formulas)) {
    
    cs_multi$models[[i]] <- 
      map(cs_multi$data, 
           ~call2('coxph', 
                  formula = cs_multi$formulas[[i]], 
                  data = .x, 
                  x = TRUE)) %>% 
      map(eval) %>% 
      map2(., cs_multi$data, 
           as_coxex)
    
  }

# Model assumptions, fit stats and inference ------
  
  insert_msg('Model assumptions, fit stats and inference')
  
  ## assumptions
  
  cs_multi$assumptions <- cs_multi$models %>% 
    map(map, summary, 'assumptions') %>% 
    map(compress, names_to = 'cohort')
  
  ## fit stats
  
  cs_multi$stats <- cs_multi$models %>% 
    map(map, summary, 'fit') %>% 
    map(compress, names_to = 'cohort')
  
  cs_multi$stats <- cs_multi$stats %>% 
    map(mutate, 
        c_lab = paste0(signif(c_index, 2), ' [', 
                       signif(lower_ci, 2), ' - ', 
                       signif(upper_ci, 2), ']'), 
        response = ifelse(cohort == 'GSE16560', 'OS', 'RFS'),  
        cohort_lab = paste(globals$study_labels[cohort], 
                           response, 
                           sep = ', '), 
        cohort_extended = paste0(cohort_lab, '\ntotal: n = ', 
                                 n_complete, ', events: n = ', 
                                 n_events),
        cohort = factor(cohort, names(cs_multi$data)))
  
  ## inference
  
  cs_multi$inference <- cs_multi$models %>% 
    map(map, summary, 'inference') %>% 
    map(compress, names_to = 'cohort')
  
  cs_multi$inference <- 
    map2(map(cs_multi$stats, ~.x[c('cohort', 'n_events')]), 
         cs_multi$inference, 
         left_join, by = 'cohort')
  
  cs_multi$inference <- cs_multi$inference %>% 
    map(mutate, 
        hr_lab = paste0(signif(estimate, 2), 
                        ' [', signif(lower_ci, 2), ' - ', 
                        signif(upper_ci, 2), ']'), 
        response = ifelse(cohort == 'GSE16560', 'OS', 'RFS'),  
        cohort_lab = paste(globals$study_labels[cohort], 
                           response, 
                           sep = ', '), 
        cohort_extended = paste0(cohort_lab, '\ntotal: n = ', 
                                 n_complete, ', events: n = ', 
                                 n_events), 
        cohort = factor(cohort, names(cs_multi$data)))
    
# Forest plots of the C-indexes -------
  
  insert_msg('Forest plots of the C-indexes')
  
  cs_multi$c_index_plot <- cs_multi$stats %>% 
    compress(names_to = 'model') %>% 
    mutate(model = factor(model, c('uni', 'clinical', 'multi'))) %>% 
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
    scale_color_manual(labels = cs_multi$model_labels, 
                       values = cs_multi$model_colors, 
                       name = '') +
    facet_grid(cohort_extended ~ ., 
               space = 'free', 
               scales = 'free') +
    globals$common_theme +
    theme(axis.title.y = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank()) + 
    labs(title = 'Predictive performance', 
         subtitle = 'Collagen Score, PSA, stage, Gleason score, Cox models', 
         x = 'C-index \u00B1 95%CI')
  
# Plotting the R-squares ------
  
  insert_msg('Plotting the R-squares')
  
  cs_multi$rsq_plot <- cs_multi$stats %>% 
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
    scale_fill_manual(values = cs_multi$model_colors,
                      labels = cs_multi$model_labels, 
                      name = '') + 
    scale_color_manual(values = cs_multi$model_colors,
                      labels = cs_multi$model_labels, 
                      name = '') + 
    scale_x_continuous(limits = c(0, 0.7), 
                       breaks = seq(0, 0.7, by = 0.1)) + 
    globals$common_theme +
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Explanatory performance', 
         subtitle = 'Collagen Score, PSA, stage, Gleason score, Cox models', 
         x = 'R\u00B2')
  
# Plotting the IBS ------
  
  insert_msg('Plotting the IBS')
  
  cs_multi$ibs_plot <- cs_multi$stats %>% 
    compress(names_to = 'model') %>% 
    mutate(model = factor(model, c('uni', 'clinical', 'multi'))) %>%
    ggplot(aes(x = 1 - ibs_model, 
               y = reorder(cohort_extended, -as.numeric(cohort)), 
               fill = model)) +
    geom_bar(stat = 'identity', 
             color = 'black', 
             position = position_dodge(0.9)) + 
    geom_text(aes(label = signif(1 - ibs_model, 2), 
                  color = model), 
              size = 2.75, 
              hjust = -0.4, 
              position = position_dodge(0.9),
              show.legend = FALSE)  +
    scale_fill_manual(values = cs_multi$model_colors,
                      labels = cs_multi$model_labels, 
                      name = '') + 
    scale_color_manual(values = cs_multi$model_colors,
                       labels = cs_multi$model_labels, 
                       name = '') + 
    scale_x_continuous(limits = c(0, 1), 
                       breaks = seq(0, 1, by = 0.25)) + 
    globals$common_theme +
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Calibration', 
         subtitle = 'Collagen Score, PSA, stage, Gleason score, Cox models', 
         x = '1 - IBS')
  
  
  
# Plotting the C-indexes and IBS as a scatter plot -------
  
  insert_msg('C_index and IBS scatter plot')
  
  cs_multi$c_ibs_plots <- cs_multi$stats %>% 
    compress(names_to = 'model') %>% 
    mutate(model = factor(model, 
                          c('uni', 'clinical', 'multi')), 
           cohort = factor(cohort, names(cs_multi$data))) %>% 
    blast(cohort)
  
  for(i in names(cs_multi$c_ibs_plots)) {
    
    cs_multi$c_ibs_plots[[i]] <- cs_multi$c_ibs_plots[[i]] %>% 
      ggplot(aes(x = c_index, 
                 y = 1 - ibs_model, 
                 color = model, 
                 size = raw_rsq)) + 
      geom_hline(yintercept = 1 - 0.5^2, 
                 linetype = 'dashed') +
      geom_vline(xintercept = 0.5, 
                 linetype = 'dashed') +
      geom_point(shape = 16) + 
      scale_color_manual(values = cs_multi$model_colors, 
                         labels = cs_multi$model_labels, 
                         name = '') + 
      scale_radius(range = c(0.5, 5), 
                   limits = c(0, 0.7),  
                   name = expression('R'^2)) + 
      globals$common_theme + 
      labs(title = paste('Predictive performance,', 
                         globals$study_labels[i]),
           subtitle = cs_multi$n_tags[[i]], 
           x = 'C-index', 
           y = '1 - IBS')
    
  }

# Plotting the hazard ratios ------
  
  insert_msg('Plotting the hazard ratios')
  
  cs_multi$hr_plot <- cs_multi$inference %>% 
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
    scale_color_manual(values = cs_multi$model_colors, 
                       labels = cs_multi$model_labels, 
                       name = '') +
    globals$common_theme +
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Inference', 
         subtitle = 'Uni-variable and multi-variable Collagen Score Cox models', 
         x = 'HR \u00B1 95%CI, Collagen Score')
  
# END -----
  
  rm(i)
  
  insert_tail()