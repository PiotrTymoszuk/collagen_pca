# Development of the collagen score in the TCGA training cohort
# The procedure: Elastic Net with the initial explanatory variable set of 55
# collagen pathway genes, response: relapse-free survival, family: Cox.
# Pre-processing of the explanatory variables: normalization Z/score
# Lambda finding in 200-repeat 10-fold CV
# The Collagen Score is defined as linear predictor score of the training
# model. 
# Pre-processing of the data sets: normalization (Z-scores, mean-centered). 
# The model explanatory factors are the firs- and secong-order log2 exprssion 
# values.

  insert_head()
  
# container ------
  
  elnet_surv <- list()
  
# globals -------
  
  insert_msg('Globals')
  
  ## first-order and second order explanatory variables
  
  elnet_surv$variables <- globals$genes_interest$gene_symbol
  
  ## analysis tables: obtained from globals
  ## inclusion of the second order terms
  
  elnet_surv$survival <- surv_globals$data %>% 
    map(select, sample_id, rfs_months, relapse)
  
  elnet_surv$expression <- surv_globals$data %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(elnet_surv$variables))
  
  elnet_surv$sq_expression <- elnet_surv$expression %>% 
    map(map_dfc, ~.x^2) %>% 
    map(~set_colnames(.x, paste0(names(.x), '_sq')))

  elnet_surv$expression <- 
    map2(elnet_surv$expression, 
         elnet_surv$sq_expression, 
         cbind) %>% 
    map(center_data, 'mean') %>% 
    map(rownames_to_column, 'sample_id')
  
  ## complete cases, adding a month to the survival time to guarantee
  ## thet there are no zero times to events
    
  elnet_surv$data <- 
    map2(elnet_surv$survival, 
         elnet_surv$expression[names(elnet_surv$survival)], 
         left_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(mutate, rfs_months = rfs_months + 1) %>% 
    map(column_to_rownames, 'sample_id')
  
  elnet_surv$variables <- 
    c(elnet_surv$variables, 
      paste0(elnet_surv$variables, '_sq'))

  ## survival objects
  
  elnet_surv$y <- elnet_surv$data %>% 
    map(~Surv(.x$rfs_months, .x$relapse))

  ## matrices of normalized explanatory variables
  
  elnet_surv$x <- elnet_surv$data %>% 
    map(~.x[elnet_surv$variables]) %>% 
    map(as.matrix)
  
# CV folds --------
  
  insert_msg('CV folds')

  elnet_surv$folds <- surv_globals$folds
  
# Tuning of the lambda parameter ------
  
  insert_msg('Lambda tuning')
  
  plan('multisession')
  
  elnet_surv$lambda_tune <- elnet_surv$folds %>% 
    future_map(~cv.glmnet(x = elnet_surv$x$geo, 
                          y = elnet_surv$y$geo, 
                          foldid = .x, 
                          type.measure = 'default', 
                          family = 'cox', 
                          alpha = 0.5), 
               .options = furrr_options(seed = TRUE))
  
  plan('sequential')
  
  elnet_surv$lambda_tbl <- elnet_surv$lambda_tune %>% 
    map(~as_tibble(.x[c('lambda', 'cvm', 'cvup', 'cvlo')])) %>% 
    map2_dfr(., elnet_surv$lambda_tune, 
             ~filter(.x, lambda == .y[['lambda.min']]))
  
  elnet_surv$opt_lambda <- elnet_surv$lambda_tbl %>% 
    filter(cvm == min(cvm))
  
# Building the training Elastic Net model ------
  
  insert_msg('Building of the training Elastic Net model')
  
  elnet_surv$glmnet_model <- 
    glmnet(x = elnet_surv$x$geo, 
           y = elnet_surv$y$geo, 
           family = 'cox', 
           alpha = 0.5, 
           lambda = elnet_surv$opt_lambda$lambda)
  
# Calculating the linear predictor scores for the training and test cohorts -------
  
  insert_msg('Calculating the collagen scores')
  
  ## prediction
  
  elnet_surv$score_tbl <- elnet_surv$x %>% 
    map(~predict(elnet_surv$glmnet_model, newx = .x)) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(set_names, c('sample_id', 'collagen_score')) %>% 
    map(as_tibble)
  
  ## appending with the survival information
  
  elnet_surv$score_tbl <- 
    map2(elnet_surv$score_tbl, 
         elnet_surv$survival, 
         left_join, by = 'sample_id')

# Building univariable Cox models ------
  
  insert_msg('Uni-variable Cox models')
  
  # working with metaprogramming, to get the entire
  # data sets kept in place with the models
  
  elnet_surv$models <- elnet_surv$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(rfs_months, relapse) ~ collagen_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         elnet_surv$score_tbl, 
         as_coxex)

# Characteristic of the Cox model: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of collagen score in the training cohort')

  ## assumptions
  
  elnet_surv$assumptions <- elnet_surv$models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  elnet_surv$stats <- elnet_surv$models %>% 
    map(summary, type = 'fit')
  
  ## inference
  
  elnet_surv$inference <- elnet_surv$models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  elnet_surv[c("assumptions", "stats", "inference")] <- 
    elnet_surv[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'geo', 'training', 'test'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  elnet_surv$calibration <- elnet_surv$models %>% 
    map(calibrate, n = 3, labels = c('low', 'int', 'high'))
  
  elnet_surv$global_cal <- elnet_surv$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo', 'training', 'test'))

  elnet_surv$brier_scores <- elnet_surv$models %>% 
    map(surv_brier)
    
# Differences in survival between the collagen score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  elnet_surv$tertile_stats <- elnet_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  elnet_surv$tertile_test <- elnet_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo', 'training', 'test')) %>% 
    re_adjust(p_variable = 'pval')

# training model estimates ------
  
  insert_msg('Estimates of the training model')
  
  elnet_surv$coefs <- coef(elnet_surv$glmnet_model) %>% 
    as.matrix %>% 
    as.data.frame %>% 
    rownames_to_column('variable') %>% 
    set_names(c('variable', 'coef')) %>% 
    mutate(exp_coef = exp(coef)) %>% 
    filter(coef != 0) %>% 
    as_tibble
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  elnet_surv$survival <- NULL
  elnet_surv$expression <- NULL
  elnet_surv$x <- NULL
  elnet_surv$y <- NULL
  elnet_surv$sq_expression <- NULL
  elnet_surv$data <- NULL
  elnet_surv$n_rep <- NULL
  elnet_surv$folds <- NULL
  
  elnet_surv <- compact(elnet_surv)
  
  save(elnet_surv, file = './cache/elnet_surv.RData')

# END -----
  
  insert_tail()