# Development and evaluation of two GBM models: 
# 1) a model with solely clinical parameters (Gleason score, tumor stage and 
# PSA at diagnosis)
#
# 2) a full model with the collagen-related gene expression levels and the 
# clinical variables listed above

  insert_head()
  
# container -----
  
  surv_multi <- list()
  
# analysis data -------
  
  insert_msg('Analysis data')
  
  ## variables: clinical and expression
  
  surv_multi$clinic_variables <- 
    c('gleason_simple', 
      'pt_stage', 
      'psa_diagnosis')
  
  surv_multi$gene_variables <- globals$genes_interest$gene_symbol
  
  ## the expression data sets
  
  surv_multi$expression <- surv_globals$data
  
  for(i in names(surv_multi$expression)) {
    
    surv_multi$expression[[i]][surv_multi$gene_variables] <- 
      surv_multi$expression[[i]][surv_multi$gene_variables] %>% 
      center_data('mean')
    
  }
  
  ## clinical data
  
  surv_multi$clinic <- surv_globals$clinic
  
  surv_multi$data <- 
    map2(surv_multi$expression, 
         surv_multi$clinic[names(surv_multi$expression)],
         left_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x)))

  ## normalization of the PSA levels 
  ## collapsing the stages: T1/T2 and T3/T4
  
  surv_multi$data <- surv_multi$data %>% 
    map(mutate, 
        psa_diagnosis = scale(psa_diagnosis)[, 1])

  surv_multi$data <- surv_multi$data %>% 
    map(mutate, 
        pt_stage = car::recode(pt_stage, 
                               "'T1' = 'T1/T2'; 'T2' = 'T1/T2';
                               'T3' = 'T3/T4'; 'T4' = 'T3/T4'"), 
        pt_stage = factor(pt_stage, c('T1/T2', 'T3/T4'))) %>% 
    map(column_to_rownames, 'sample_id')
  
  ## GBM scores for the the collagen-only survival models
  
  surv_multi$collagen_scores <- gbm_surv$score_tbl %>% 
    map(select, sample_id, gbm_score) %>% 
    map(set_names, c('sample_id', 'collagen_score'))
  
# Model formulas and tune grids  --------
  
  insert_msg('Model formulas and tune grids')
  
  ## formulas of the training cohort multi-parameter models: 
  ## 1) a Cox model with the traditional clinical markers only
  ## 2) a Cox model with the GBM score and the clinical markers
  
  surv_multi$gbm_formulas$clinic <- 
    paste0('Surv(rfs_months, relapse) ~', 
           paste(surv_multi$clinic_variables, collapse = ' + '))
    
  surv_multi$gbm_formulas$full <- 'Surv(rfs_months, relapse) ~ .'
  
  surv_multi$gbm_formulas <- surv_multi$gbm_formulas %>% 
    map(as.formula)
  
  ## formulas of the univariable Cox models used for evaluation
  
  surv_multi$cox_formulas <- 
    list(clinic = Surv(rfs_months, relapse) ~ clinic_score, 
         collagen = Surv(rfs_months, relapse) ~ collagen_score, 
         full = Surv(rfs_months, relapse) ~ full_score)
  
  ## tune grid
  
  surv_multi$tune_grid <- surv_globals$gbm_grid

# N numbers -------
  
  insert_msg('N numbers of complete cases and events')
  
  surv_multi$n_numbers <- surv_multi$data %>% 
    map(~tibble(n_total = nrow(.x), 
                n_events = sum(.x$relapse))) %>% 
    compress(names_to = 'cohort')
  
# Characteristic and comparison of the cohorts ------
  
  insert_msg('Characteristic and comparison of the clinics')
  
  ## descriptive stats
  
  surv_multi$exploration$stats <- surv_multi$data %>% 
    compress(names_to = 'cohort') %>% 
    explore(variables = surv_multi$clinic_variables, 
            split_factor = 'cohort', 
            what = 'table', 
            pub_styled = TRUE) %>% 
    format_desc
  
  ## comparison: non-parametric tests
  
  surv_multi$exploration$test <- surv_multi$data %>% 
    compress(names_to = 'cohort') %>% 
    compare_variables(variables = surv_multi$clinic_variables, 
                      split_factor = 'cohort', 
                      what = 'eff_size', 
                      types = c('cramer_v', 
                                'cramer_v', 
                                'kruskal_etasq'), 
                      exact = FALSE, 
                      ci = FALSE, 
                      pub_styled = TRUE, 
                      adj_method = 'BH')
  
# CV tuning in the GEO data set -------
  
  insert_msg('CV tuning in the GEO data set')
  
  set.seed(1234)
  
  surv_multi$tuning <- gbm_tune(data = surv_multi$data$geo, 
                                time_variable = 'rfs_months', 
                                event_variable = 'relapse', 
                                n_folds = 10, 
                                tune_grid = surv_multi$tune_grid, 
                                distribution = 'coxph')
  
# Training of the multi-parameter models, pooled GEO cohort -------
  
  insert_msg('Training of the multi-parameter models, pooled GEO cohort')
  
  ## the clinical-only model is trained with the default settings
  
  set.seed(1234)
  
  surv_multi$gbm_models$clinic <- 
    gbm(formula = surv_multi$gbm_formulas$clinic, 
        data = surv_multi$data$geo, 
        distribution = 'coxph', 
        cv.folds = 10)
  
  ## the full model trained with the CV-tuned parameters
  
  set.seed(1234)
  
  surv_multi$gbm_models$full <- 
    gbm(formula = surv_multi$gbm_formulas$full, 
        data = surv_multi$data$geo, 
        distribution = 'coxph', 
        n.trees = surv_multi$tuning$best_tune$n.trees[1],
        shrinkage = surv_multi$tuning$best_tune$shrinkage[1],
        interaction.depth = surv_multi$tuning$best_tune$interaction.depth[1],
        n.minobsinnode = surv_multi$tuning$best_tune$n.minobsinnode[1], 
        cv.folds = 10)
  
  ## best iterations
  
  surv_multi$best_iter <- surv_multi$gbm_models %>% 
    map(gbm.perf, method = 'cv')
  
# Predictions --------
  
  insert_msg('Predictions: LP scores')
  
  for(i in names(surv_multi$data)) {
    
    surv_multi$score_tbl[[i]] <- 
      list(object = surv_multi$gbm_models, 
           n.trees = surv_multi$best_iter) %>% 
      pmap(predict, 
           newdata = surv_multi$data[[i]]) %>% 
      reduce(cbind) %>% 
      set_colnames(paste0(names(surv_multi$gbm_models), '_score')) %>% 
      as_tibble
    
    surv_multi$score_tbl[[i]] <-
      cbind(surv_multi$data[[i]][c('rfs_months', 'relapse')], 
            surv_multi$score_tbl[[i]]) %>% 
      rownames_to_column('sample_id') %>% 
      as_tibble
    
    surv_multi$score_tbl[[i]] <- 
      left_join(surv_multi$score_tbl[[i]], 
                surv_multi$collagen_scores[[i]], 
                by = 'sample_id')
    
  }
  
# Univariable Cox models used for evaluation --------
  
  insert_msg('Univariable Cox model used for evaluation of predictions')
  
  for(i in names(surv_multi$score_tbl)) {
    
    surv_multi$cox_models[[i]] <- surv_multi$cox_formulas %>% 
      map(~call2(.fn = 'coxph', 
                 formula = .x, 
                 data = surv_multi$score_tbl[[i]], 
                 x = TRUE, 
                 y = TRUE)) %>% 
      map(eval) %>% 
      map(as_coxex, data = surv_multi$score_tbl[[i]])
    
  }
  
# Assumptions of the univariable models, fit stats and inference ------
  
  insert_msg('Assumptions, fit stats and inference')
  
  ## there are major violations of the PH assumption, 
  ## one can think of fitting a model with spline or strata of the linear
  ## predictor scores. Because we're interested in the performance only
  ## and the key statistic, C-index, is anyway non-parametric, we're
  ## sticking to the firs-order Cox models
  
  surv_multi$assumptions <- surv_multi$cox_models %>% 
    map(map, summary, 'assumptions')
  
  ## fit stats
  
  surv_multi$stats <- surv_multi$cox_models %>% 
    map(map, summary, 'fit')
  
  ## inference 
  
  surv_multi$inference <- surv_multi$cox_models %>% 
    map(map, summary, 'inference')
  
  ## collapsing by the cohort and model type
  
  surv_multi[c("stats", "inference")] <- 
    surv_multi[c("stats", "inference")] %>% 
    map(map, compress, names_to = 'model') %>% 
    map(compress, names_to = 'cohort') %>% 
    map(mutate, 
        model = factor(model, names(surv_multi$stats[[1]])), 
        cohort = factor(cohort, names(surv_multi$data)), 
        dataset = ifelse(cohort == 'geo', 'training', 'test'))
  
  surv_multi$inference <- surv_multi$inference %>% 
    mutate(estimate = exp(estimate), 
           lower_ci = exp(lower_ci), 
           upper_ci = exp(upper_ci))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  surv_multi$calibration <- surv_multi$cox_models %>% 
    map(map, calibrate, n = 3, labels = c('low', 'int', 'high'))
  
  surv_multi$global_cal <- surv_multi$calibration %>% 
    map(map, summary, type = 'global') %>% 
    map(compress, names_to = 'model') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(model = factor(model, levels(surv_multi$stats$model)), 
           cohort = factor(cohort, levels(surv_multi$stats$cohort)), 
           dataset = ifelse(cohort == 'geo', 'training', 'test'))
  
  surv_multi$brier_scores <- surv_multi$cox_models %>% 
    map(map, surv_brier)
  
# Differences in survival between the collagen score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  surv_multi$tertile_stats <- surv_multi$calibration %>% 
    map(map, ~.$surv_fit) %>% 
    map(map, surv_median)
  
  surv_multi$tertile_test <- surv_multi$calibration %>% 
    map(map, ~.$surv_fit) %>% 
    map(map, surv_pvalue, method = 'S1') %>% 
    map(compress, names_to = 'model') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(model = factor(model, levels(surv_multi$stats$model)), 
           cohort = factor(cohort, levels(surv_multi$stats$cohort)), 
           dataset = ifelse(cohort == 'geo', 'training', 'test')) %>% 
    re_adjust('pval')
  
# Importance stats ------
  
  insert_msg('Variable importance')
  
  surv_multi$importance <- surv_multi$gbm_models %>% 
    map(summary) %>% 
    map(set_names, c('variable', 'rel_influence')) %>% 
    map(as_tibble)
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  surv_multi$survival <- NULL
  surv_multi$clinic <- NULL
  surv_multi$data <- NULL
  surv_multi$collagen_scores <- NULL
  
  surv_multi$gbm_formulas <- NULL
  surv_multi$cox_formulas <- NULL
  
  surv_multi <- compact(surv_multi)
  
  save(surv_multi, file = './cache/surv_multi.RData')
  
# END --------
  
  rm(i)
  
  insert_tail()