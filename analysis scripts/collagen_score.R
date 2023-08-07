# Development of the collagen score in the TCGA training cohort
# The procedure: Elastic Net with the initial explanatory variable set of 28
# collagen pathway genes, response: relapse-free survival, family: Cox.
# Pre-processing of the explanatory variables: normalization Z/score
# Lambda finding in 200-repeat 10-fold CV
# The Collagen Score is defined as linear predictor score of the training
# model
#
# I'm working with the Combat-normalized dataset with suppressed cohort effects

  insert_head()
  
# container ------
  
  coll_score <- list()
  
# globals -------
  
  insert_msg('Globals')
  
  ## first-order explanatory variables
  ## inclusion of second order variables does not improve the fit
  ## by increases the computation time by factor 4 - 16
  
  coll_score$variables <- globals$genes_interest$gene_symbol
  
  ## analysis tables: survival information
  
  coll_score$analysis_tbl <- study_data %>% 
    map(~.x$expression) %>% 
    map(extract_tumor_samples) %>% 
    map(select, 
        patient_id, 
        any_of(c('death', 'vitality_fup')), 
        any_of(c('relapse', 'relapse_fup')))
  
  ## ... and batch-normalized expression data of the collagen genes
  
  coll_score$analysis_tbl <- 
    map2(coll_score$analysis_tbl, 
         combatch$adjusted_data, 
         left_join, by = 'patient_id') %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(column_to_rownames, 'patient_id')
  
  ## survival object for the training cohort
  
  coll_score$tcga_surv <- Surv(coll_score$analysis_tbl$tcga$relapse_fup, 
                               coll_score$analysis_tbl$tcga$relapse)
  
  ## matrices of normalized explanatory variables
  
  coll_score$x_mats <- coll_score$analysis_tbl %>% 
    map(~.x[coll_score$variables]) %>% 
    map(center_data, type = 'mean') %>% 
    map(as.matrix)
  
  ## folds 
  
  set.seed(1234)
  
  coll_score$n_rep <- 200
  
  coll_score$folds <- 1:coll_score$n_rep %>% 
    map(function(x) createFolds(y = coll_score$analysis_tbl$tcga[[1]], 
                                k = 10, 
                                list = FALSE, 
                                returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:coll_score$n_rep))
  
# Tuning of the lambda parameter ------
  
  insert_msg('Lambda tuning')
  
  plan('multisession')
  
  coll_score$lambda_tune <- coll_score$folds %>% 
    future_map(~cv.glmnet(x = coll_score$x_mats$tcga, 
                          y = coll_score$tcga_surv, 
                          foldid = .x, 
                          type.measure = 'default', 
                          family = 'cox', 
                          alpha = 0.5), 
               .options = furrr_options(seed = TRUE))
  
  plan('sequential')
  
  coll_score$lambda_tbl <- coll_score$lambda_tune %>% 
    map(~as_tibble(.x[c('lambda', 'cvm', 'cvup', 'cvlo')])) %>% 
    map2_dfr(., coll_score$lambda_tune, 
             ~filter(.x, lambda == .y[['lambda.min']]))
  
  coll_score$opt_lambda <- coll_score$lambda_tbl %>% 
    filter(cvm == min(cvm))
  
# Building the training Elastic Net model ------
  
  insert_msg('Building of the training Elastic Net model')
  
  coll_score$tcga_lasso <- 
    glmnet(x = coll_score$x_mats$tcga, 
           y = coll_score$tcga_surv, 
           family = 'cox', 
           alpha = 0.5, 
           lambda = coll_score$opt_lambda$lambda)
  
# Calculating the linear predictor scores for the training and test cohorts -------
  
  insert_msg('Calculating the collagen scores')
  
  ## prediction
  
  coll_score$score_tbl <- coll_score$x_mats %>% 
    map(~predict(coll_score$tcga_lasso, newx = .x)) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'patient_id') %>% 
    map(set_names, c('patient_id', 'collagen_score')) %>% 
    map(as_tibble)
  
  ## appending with the survival information
  
  coll_score$score_tbl <- coll_score$analysis_tbl %>% 
    map(rownames_to_column, 'patient_id') %>% 
    map(select, 
        patient_id, 
        any_of(c('death', 'vitality_fup')), 
        any_of(c('relapse', 'relapse_fup'))) %>% 
    map2(coll_score$score_tbl, ., 
         left_join, by = 'patient_id')

# Building univariable Cox models ------
  
  insert_msg('Uni-variable Cox models')
  
  ## GSE16560 has nor relapse data!!!
  # working with metaprogramming, to get the entire
  # data sets kept in place with the models
  
  coll_score$models <- 
    map2(coll_score$score_tbl, 
         list(Surv(vitality_fup, death) ~ collagen_score, 
              Surv(relapse_fup, relapse) ~ collagen_score, 
              Surv(relapse_fup, relapse) ~ collagen_score, 
              Surv(relapse_fup, relapse) ~ collagen_score, 
              Surv(relapse_fup, relapse) ~ collagen_score), 
         ~call2('coxph', 
                formula = .y, 
                data = .x, 
                x = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         coll_score$score_tbl, 
         as_coxex)

# Characteristic of the Cox model: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of collagen score in the training cohort')

  ## assumptions
  
  coll_score$assumptions <- coll_score$models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  coll_score$stats <- coll_score$models %>% 
    map(summary, type = 'fit')
  
  ## inference
  
  coll_score$inference <- coll_score$models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  coll_score[c("assumptions", "stats", "inference")] <- 
    coll_score[c("assumptions", "stats", "inference")] %>% 
    map(compress, names_to = 'cohort')
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  coll_score$calibration <- coll_score$models %>% 
    map(calibrate, n = 3, labels = c('low', 'int', 'high'))
  
  coll_score$global_cal <- coll_score$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort')

  coll_score$brier_scores <- coll_score$models %>% 
    map(surv_brier)
    
# Differences in survival between the collagen score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  coll_score$tertile_stats <- coll_score$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  coll_score$tertile_test <- coll_score$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    re_adjust(p_variable = 'pval')

# training model estimates ------
  
  insert_msg('Estimates of the training model')
  
  coll_score$coefs <- coef(coll_score$tcga_lasso) %>% 
    as.matrix %>% 
    as.data.frame %>% 
    rownames_to_column('variable') %>% 
    set_names(c('variable', 'coef')) %>% 
    mutate(exp_coef = exp(coef)) %>% 
    filter(coef != 0) %>% 
    as_tibble

# END -----
  
  insert_tail()