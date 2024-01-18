# Development of the collagen score in the pooled GEO training cohort
# The procedure: Ridge with the initial explanatory variable set of 55
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
  
  ridge_surv <- list()
  
# globals -------
  
  insert_msg('Globals')
  
  ## first-order and second order explanatory variables
  
  ridge_surv$variables <- globals$genes_interest$gene_symbol
  
  ## analysis tables: obtained from globals
  ## inclusion of the second order terms
  
  ridge_surv$survival <- surv_globals$data %>% 
    map(select, sample_id, rfs_months, relapse)
  
  ridge_surv$expression <- surv_globals$data %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(ridge_surv$variables))
  
  ridge_surv$sq_expression <- ridge_surv$expression %>% 
    map(map_dfc, ~.x^2) %>% 
    map(~set_colnames(.x, paste0(names(.x), '_sq')))

  ridge_surv$expression <- 
    map2(ridge_surv$expression, 
         ridge_surv$sq_expression, 
         cbind) %>% 
    map(center_data, 'mean') %>% 
    map(rownames_to_column, 'sample_id')
  
  ## complete cases, adding a month to the survival time to guarantee
  ## thet there are no zero times to events
    
  ridge_surv$data <- 
    map2(ridge_surv$survival, 
         ridge_surv$expression[names(ridge_surv$survival)], 
         left_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(mutate, rfs_months = rfs_months + 1) %>% 
    map(column_to_rownames, 'sample_id')
  
  ridge_surv$variables <- 
    c(ridge_surv$variables, 
      paste0(ridge_surv$variables, '_sq'))

  ## survival objects
  
  ridge_surv$y <- ridge_surv$data %>% 
    map(~Surv(.x$rfs_months, .x$relapse))

  ## matrices of normalized explanatory variables
  
  ridge_surv$x <- ridge_surv$data %>% 
    map(~.x[ridge_surv$variables]) %>% 
    map(as.matrix)
  
# CV folds --------
  
  insert_msg('CV folds')

  ridge_surv$folds <- surv_globals$folds
  
# Tuning of the lambda parameter ------
  
  insert_msg('Lambda tuning')
  
  plan('multisession')
  
  ridge_surv$lambda_tune <- ridge_surv$folds %>% 
    future_map(~cv.glmnet(x = ridge_surv$x$geo, 
                          y = ridge_surv$y$geo, 
                          foldid = .x, 
                          type.measure = 'default', 
                          family = 'cox', 
                          alpha = 0), 
               .options = furrr_options(seed = TRUE))
  
  plan('sequential')
  
  ridge_surv$lambda_tbl <- ridge_surv$lambda_tune %>% 
    map(~as_tibble(.x[c('lambda', 'cvm', 'cvup', 'cvlo')])) %>% 
    map2_dfr(., ridge_surv$lambda_tune, 
             ~filter(.x, lambda == .y[['lambda.min']]))
  
  ridge_surv$opt_lambda <- ridge_surv$lambda_tbl %>% 
    filter(cvm == min(cvm))
  
# Building the training Elastic Net model ------
  
  insert_msg('Building of the training Elastic Net model')
  
  ridge_surv$glmnet_model <- 
    glmnet(x = ridge_surv$x$geo, 
           y = ridge_surv$y$geo, 
           family = 'cox', 
           alpha = 0.5, 
           lambda = ridge_surv$opt_lambda$lambda)
  
# Calculating the linear predictor scores for the training and test cohorts -------
  
  insert_msg('Calculating the collagen scores')
  
  ## prediction
  
  ridge_surv$score_tbl <- ridge_surv$x %>% 
    map(~predict(ridge_surv$glmnet_model, newx = .x)) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(set_names, c('sample_id', 'collagen_score')) %>% 
    map(as_tibble)
  
  ## appending with the survival information
  
  ridge_surv$score_tbl <- 
    map2(ridge_surv$score_tbl, 
         ridge_surv$survival, 
         left_join, by = 'sample_id')

# Building univariable Cox models ------
  
  insert_msg('Uni-variable Cox models')
  
  # working with metaprogramming, to get the entire
  # data sets kept in place with the models
  
  ridge_surv$models <- ridge_surv$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(rfs_months, relapse) ~ collagen_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         ridge_surv$score_tbl, 
         as_coxex)

# Characteristic of the Cox model: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of collagen score in the training cohort')

  ## assumptions
  
  ridge_surv$assumptions <- ridge_surv$models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  ridge_surv$stats <- ridge_surv$models %>% 
    map(summary, type = 'fit')
  
  ## inference
  
  ridge_surv$inference <- ridge_surv$models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  ridge_surv[c("assumptions", "stats", "inference")] <- 
    ridge_surv[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'geo', 'training', 'test'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  ridge_surv$calibration <- ridge_surv$models %>% 
    map(calibrate, n = 3, labels = c('low', 'int', 'high'))
  
  ridge_surv$global_cal <- ridge_surv$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo', 'training', 'test'))

  ridge_surv$brier_scores <- ridge_surv$models %>% 
    map(surv_brier)
    
# Differences in survival between the collagen score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  ridge_surv$tertile_stats <- ridge_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  ridge_surv$tertile_test <- ridge_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo', 'training', 'test')) %>% 
    re_adjust(p_variable = 'pval')

# training model estimates ------
  
  insert_msg('Estimates of the training model')
  
  ridge_surv$coefs <- coef(ridge_surv$glmnet_model) %>% 
    as.matrix %>% 
    as.data.frame %>% 
    rownames_to_column('variable') %>% 
    set_names(c('variable', 'coef')) %>% 
    mutate(exp_coef = exp(coef)) %>% 
    filter(coef != 0) %>% 
    as_tibble
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  ridge_surv$survival <- NULL
  ridge_surv$expression <- NULL
  ridge_surv$x <- NULL
  ridge_surv$y <- NULL
  ridge_surv$sq_expression <- NULL
  ridge_surv$data <- NULL
  ridge_surv$n_rep <- NULL
  ridge_surv$folds <- NULL
  
  ridge_surv <- compact(ridge_surv)
  
  save(ridge_surv, file = './cache/ridge_surv.RData')

# END -----
  
  insert_tail()