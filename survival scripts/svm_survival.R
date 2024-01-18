# Modeling of biochemical relapse-free survival as a function of expression 
# of the collagen-related genes in the cancer tissue with 
# Support Vector Machines.
#
# Pre-processing: normalization, mean centering.
# Tuning and training: in the pooled GEO cohort.
# 

  insert_head()
  
# container ------
  
  svm_surv <- list()
  
# modeling data -------
  
  insert_msg('Modeling data')
  
  ## only first-order terms
  
  svm_surv$variables <- globals$genes_interest$gene_symbol
  
  ## analysis tables: from the globals
  ## normalization and centering of the explanatory variables
  ## min-max normalization of the survival time: avoiding zeros
  ## by a small shift
  
  svm_surv$data <- surv_globals$data
  
  for(i in names(svm_surv$data)) {
    
    svm_surv$data[[i]][svm_surv$variables] <- 
      svm_surv$data[[i]][svm_surv$variables] %>% 
      center_data('mean')
    
    svm_surv$data[[i]]['rfs_months'] <- 
      svm_surv$data[[i]]['rfs_months'] %>% 
      min_max
    
  }
  
  svm_surv$data <- svm_surv$data %>% 
    map(mutate, rfs_months = rfs_months + 0.01) %>% 
    map(column_to_rownames, 'sample_id')
  
# Modeling globals ---------
  
  insert_msg('Modeling globals: CV folds and tuning data frames')
  
  ## folds 
  
  set.seed(1234)
  
  svm_surv$n_rep <- 5
  
  svm_surv$folds <- 1:svm_surv$n_rep %>% 
    map(function(x) createFolds(y = factor(svm_surv$data$geo$relapse), 
                                k = 10, 
                                list = TRUE, 
                                returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:svm_surv$n_rep)) %>% 
    unlist(recursive = FALSE)
  
  ## tune grids: the additive kernel and vanbelle1 method tends to function
  ## the best in terms of C-index as tested per hand. We're tuning
  ## hence just the gamma cost parameter
  
  svm_surv$tune_grid <- 
    expand.grid(gamma.mu = c(0.0001, 0.0025, 0.005, 0.01, 
                             0.015, 0.02, 0.025, 0.05, 
                             seq(0.1, 1, by = 0.1)), 
                type = c('vanbelle1'), 
                diff.meth = c('makediff3'), 
                kernel = c('add_kernel'), 
                stringsAsFactors = FALSE)
  
# Tuning of the GEO model -------
  
  insert_msg('Tuning of the GEO model')
  
  plan('multisession')
  
  svm_surv$tuning <- svm_tune(data = svm_surv$data$geo, 
                              time_variable = 'rfs_months', 
                              event_variable = 'relapse', 
                              folds = svm_surv$folds, 
                              tune_grid = svm_surv$tune_grid)
  
  plan('sequential')

# Training of the SVM model in the pooled GEO cohort -------
  
  insert_msg('Training of the GEO model')
  
  svm_surv$svm_model <- 
    survivalsvm(formula = Surv(rfs_months, relapse) ~ ., 
                data = svm_surv$data$geo, 
                type = svm_surv$tuning$best_tune$type[[1]], 
                diff.meth = svm_surv$tuning$best_tune$diff.meth[[1]], 
                gamma.mu = svm_surv$tuning$best_tune$gamma.mu[[1]], 
                kernel = svm_surv$tuning$best_tune$kernel[[1]])
  
# Predictions: SVM scores --------
  
  insert_msg('Predictions: SVM scores')
  
  svm_surv$predictions <- svm_surv$data %>% 
    map(predict, object = svm_surv$svm_model)
  
  svm_surv$score_tbl <- 
    map2(svm_surv$predictions, 
         svm_surv$data, 
         svm_score)
  
# Univariable Cox models with the tertiles of SVM score as explanatory variable -------
  
  insert_msg('Univariable Cox models')
  
  svm_surv$cox_models <- svm_surv$score_tbl %>% 
    map(~call2(.fn = 'coxph', 
               formula = Surv(rfs_months, relapse) ~ score_cuts, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., svm_surv$score_tbl, as_coxex)
  
# Assumptions, fit stats and inference -------
  
  insert_msg('Assumptions, fit stats and inference')
  
  ## assumptions
  
  svm_surv$assumptions <- svm_surv$cox_models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  svm_surv$stats <- svm_surv$cox_models %>% 
    map(summary, type = 'fit')
  
  ## inference
  
  svm_surv$inference <- svm_surv$cox_models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  svm_surv[c("assumptions", "stats", "inference")] <- 
    svm_surv[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'geo', 'training', 'test'))
  
# Differences in survival between the collagen score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')

  svm_surv$tertile_test <- svm_surv$score_tbl %>% 
    map(survminer::surv_fit, 
        formula = Surv(rfs_months, relapse) ~ score_cuts) %>%  
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo', 'training', 'test')) %>% 
    re_adjust(p_variable = 'pval')
  
# Brier scores for the unique time points -------
  
  insert_msg('Brier scores for the unique time points')
  
  svm_surv$brier_scores <- svm_surv$cox_models %>% 
    map(surv_brier)
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  svm_surv$survival <- NULL
  svm_surv$expression <- NULL
  svm_surv$n_rep <- NULL

  svm_surv <- compact(svm_surv)
  
  save(svm_surv, file = './cache/svm_surv.RData')
  
# END -----
  
  insert_tail()