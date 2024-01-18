# Performance statistics, Brier scores and survival in score tertiles 
# for multi-parameter
# relapse-free survival models fitted with the Ridge, Elastic Net, LASSO, SVM 
# and Random Forest algorithms.

  insert_head()
  
# container -----
  
  surv_summary <- list()
  
# Data frames with statistics of performance and variable importance -------
  
  insert_msg('Performance stats and variable importance')
  
  ## performance stats
  
  surv_summary$stats <- list(ridge = ridge_surv, 
                             elnet = elnet_surv, 
                             lasso = lasso_surv, 
                             svm = svm_surv, 
                             rf = rf_surv, 
                             gbm = gbm_surv) %>% 
    map(~.x$stats) %>% 
    map(select, 
        dataset, cohort, c_index, ibs_model, ibs_reference)
  
  ## variable importance
  
  surv_summary$importance <- list(ridge = ridge_surv$coefs, 
                                  elnet = elnet_surv$coefs,
                                  lasso = lasso_surv$coefs, 
                                  svm = svm_imp$test, 
                                  rf = rf_surv$importance$test, 
                                  gbm = gbm_surv$importance) %>% 
    map(filter, variable != 'full')
  
# Data frames with Brier scores for unique time points ------
  
  insert_msg('Brier scores for the time points')
  
  surv_summary$brier_scores <- list(ridge = ridge_surv, 
                                    elnet = elnet_surv, 
                                    lasso = lasso_surv, 
                                    svm = svm_surv, 
                                    rf = rf_surv, 
                                    gbm = gbm_surv) %>% 
    map(~.x$brier_scores)
  
  ## compatible format
  
  surv_summary$brier_scores$rf <- 
    surv_summary$brier_scores$rf %>% 
    map(mutate, 
        training = test, 
        test = NA)

# Score tertiles and tests for differences between the score tertiles -------
  
  insert_msg('Tertiles and tertile tests')
  
  ## this data is available only for the Cox-like models and SVM
  ## but not for the Random Forests

  surv_summary$tertile_data <- list(ridge = ridge_surv, 
                                    elnet = elnet_surv, 
                                    lasso = lasso_surv, 
                                    svm = svm_surv, 
                                    gbm = gbm_surv) %>% 
    map(~.x$score_tbl) %>% 
    map2(., c(rep('collagen_score', 3), 'svm_score', 'gbm_score'), 
         function(data_lst, var) data_lst %>% 
           map(mutate, score_cuts = cut_tertiles(.data[[var]]))) %>% 
    map(map, 
        select, 
        sample_id, rfs_months, relapse, score_cuts)
  
  ## tertile N numbers: total and events
  
  surv_summary$tertile_n <- surv_summary$tertile_data %>% 
    map(map, count, score_cuts) %>% 
    map(map, set_names, c('score_cuts', 'n_total'))
  
  surv_summary$tertile_events <- surv_summary$tertile_data %>% 
    map(map, filter, relapse == 1) %>% 
    map(map, count, score_cuts) %>% 
    map(map, set_names, c('score_cuts', 'n_events'))
  
  surv_summary$tertile_n <- 
    map2(surv_summary$tertile_n, 
         surv_summary$tertile_events, 
         function(x, y) map2(x, y, left_join, by = 'score_cuts'))
  
  ## tertile surv fits, median survival times and p values
  
  surv_summary$tertile_fits <- surv_summary$tertile_data %>% 
    map(map, 
        survminer::surv_fit, 
        formula = Surv(rfs_months, relapse) ~ score_cuts)
  
  surv_summary$tertile_stats <- surv_summary$tertile_fits %>% 
    map(map, surv_median) %>% 
    map(compress, names_to = 'cohort')
  
  surv_summary$tertile_test <- surv_summary$tertile_fits %>% 
    map(map, surv_pvalue, methor = 'S1') %>% 
    map(map, re_adjust, 'pval') %>% 
    map(compress, names_to = 'cohort')
  
# END ------
  
  surv_summary$tertile_events <- NULL
  
  surv_summary <- compact(surv_summary)
  
  insert_tail()