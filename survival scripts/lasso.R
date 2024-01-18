# Development of the collagen score in the pooled GEO training cohort
# The procedure: LASSO with the initial explanatory variable set of 55
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
  
  lasso_surv <- list()
  
# globals -------
  
  insert_msg('Globals')
  
  ## first-order and second order explanatory variables
  
  lasso_surv$variables <- globals$genes_interest$gene_symbol
  
  ## analysis tables: obtained from globals
  ## inclusion of the second order terms
  
  lasso_surv$survival <- surv_globals$data %>% 
    map(select, sample_id, rfs_months, relapse)
  
  lasso_surv$expression <- surv_globals$data %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(lasso_surv$variables))
  
  lasso_surv$sq_expression <- lasso_surv$expression %>% 
    map(map_dfc, ~.x^2) %>% 
    map(~set_colnames(.x, paste0(names(.x), '_sq')))

  lasso_surv$expression <- 
    map2(lasso_surv$expression, 
         lasso_surv$sq_expression, 
         cbind) %>% 
    map(center_data, 'mean') %>% 
    map(rownames_to_column, 'sample_id')
  
  ## complete cases, adding a month to the survival time to guarantee
  ## thet there are no zero times to events
    
  lasso_surv$data <- 
    map2(lasso_surv$survival, 
         lasso_surv$expression[names(lasso_surv$survival)], 
         left_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(mutate, rfs_months = rfs_months + 1) %>% 
    map(column_to_rownames, 'sample_id')
  
  lasso_surv$variables <- 
    c(lasso_surv$variables, 
      paste0(lasso_surv$variables, '_sq'))

  ## survival objects
  
  lasso_surv$y <- lasso_surv$data %>% 
    map(~Surv(.x$rfs_months, .x$relapse))

  ## matrices of normalized explanatory variables
  
  lasso_surv$x <- lasso_surv$data %>% 
    map(~.x[lasso_surv$variables]) %>% 
    map(as.matrix)
  
# CV folds --------
  
  insert_msg('CV folds')

  lasso_surv$folds <- surv_globals$folds
  
# Tuning of the lambda parameter ------
  
  insert_msg('Lambda tuning')
  
  plan('multisession')
  
  lasso_surv$lambda_tune <- lasso_surv$folds %>% 
    future_map(~cv.glmnet(x = lasso_surv$x$geo, 
                          y = lasso_surv$y$geo, 
                          foldid = .x, 
                          type.measure = 'default', 
                          family = 'cox', 
                          alpha = 1), 
               .options = furrr_options(seed = TRUE))
  
  plan('sequential')
  
  lasso_surv$lambda_tbl <- lasso_surv$lambda_tune %>% 
    map(~as_tibble(.x[c('lambda', 'cvm', 'cvup', 'cvlo')])) %>% 
    map2_dfr(., lasso_surv$lambda_tune, 
             ~filter(.x, lambda == .y[['lambda.min']]))
  
  lasso_surv$opt_lambda <- lasso_surv$lambda_tbl %>% 
    filter(cvm == min(cvm))
  
# Building the training LASSO model ------
  
  insert_msg('Building of the training LASSO model')
  
  lasso_surv$glmnet_model <- 
    glmnet(x = lasso_surv$x$geo, 
           y = lasso_surv$y$geo, 
           family = 'cox', 
           alpha = 1, 
           lambda = lasso_surv$opt_lambda$lambda)
  
# Calculating the linear predictor scores for the training and test cohorts -------
  
  insert_msg('Calculating the collagen scores')
  
  ## prediction
  
  lasso_surv$score_tbl <- lasso_surv$x %>% 
    map(~predict(lasso_surv$glmnet_model, newx = .x)) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(set_names, c('sample_id', 'collagen_score')) %>% 
    map(as_tibble)
  
  ## appending with the survival information
  
  lasso_surv$score_tbl <- 
    map2(lasso_surv$score_tbl, 
         lasso_surv$survival, 
         left_join, by = 'sample_id')

# Building univariable Cox models ------
  
  insert_msg('Uni-variable Cox models')
  
  # working with metaprogramming, to get the entire
  # data sets kept in place with the models
  
  lasso_surv$models <- lasso_surv$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(rfs_months, relapse) ~ collagen_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         lasso_surv$score_tbl, 
         as_coxex)

# Characteristic of the Cox model: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of collagen score in the training cohort')

  ## assumptions
  
  lasso_surv$assumptions <- lasso_surv$models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  lasso_surv$stats <- lasso_surv$models %>% 
    map(summary, type = 'fit')
  
  ## inference
  
  lasso_surv$inference <- lasso_surv$models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  lasso_surv[c("assumptions", "stats", "inference")] <- 
    lasso_surv[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'geo', 'training', 'test'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  lasso_surv$calibration <- lasso_surv$models %>% 
    map(calibrate, n = 3, labels = c('low', 'int', 'high'))
  
  lasso_surv$global_cal <- lasso_surv$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo', 'training', 'test'))

  lasso_surv$brier_scores <- lasso_surv$models %>% 
    map(surv_brier)
    
# Differences in survival between the collagen score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  lasso_surv$tertile_stats <- lasso_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  lasso_surv$tertile_test <- lasso_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo', 'training', 'test')) %>% 
    re_adjust(p_variable = 'pval')

# training model estimates ------
  
  insert_msg('Estimates of the training model')
  
  lasso_surv$coefs <- coef(lasso_surv$glmnet_model) %>% 
    as.matrix %>% 
    as.data.frame %>% 
    rownames_to_column('variable') %>% 
    set_names(c('variable', 'coef')) %>% 
    mutate(exp_coef = exp(coef)) %>% 
    filter(coef != 0) %>% 
    as_tibble
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  lasso_surv$survival <- NULL
  lasso_surv$expression <- NULL
  lasso_surv$x <- NULL
  lasso_surv$y <- NULL
  lasso_surv$sq_expression <- NULL
  lasso_surv$data <- NULL
  lasso_surv$n_rep <- NULL
  lasso_surv$folds <- NULL
  
  lasso_surv <- compact(lasso_surv)
  
  save(lasso_surv, file = './cache/lasso_surv.RData')

# END -----
  
  insert_tail()