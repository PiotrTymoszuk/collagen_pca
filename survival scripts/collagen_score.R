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
  
  coll_score <- list()
  
# globals -------
  
  insert_msg('Globals')
  
  ## first-order explanatory variables
  ## inclusion of second order variables does not improve the fit
  ## by increases the computation time by factor 4 - 16
  
  coll_score$variables <- globals$genes_interest$gene_symbol
  
  ## analysis tables: survival information
  ## for curiosity I'm checking the association of the collagen score
  ## with overall survival in the GSE16560 cohort as well (names of the 
  ## survival responses changed for consistency with the other cohorts)
  
  coll_score$survival <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(safely(select), sample_id, relapse, rfs_months) %>% 
    map(~.x$result) %>% 
    compact
  
  coll_score$survival$gse16560 <- gse16560$clinic %>% 
    filter(tissue_type == 'tumor') %>% 
    select(sample_id, death, os_months) %>% 
    set_names(c('sample_id', 'relapse', 'rfs_months'))
  
  coll_score$survival <- 
    coll_score$survival[c("gse16560", "gse54460", 
                          "gse70768", "gse70769", 
                          "gse220095", "tcga", 
                          "dkfz")]

  ## inclusion of the second order terms
  
  coll_score$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(coll_score$variables))
  
  coll_score$sq_expression <- coll_score$expression %>% 
    map(map_dfc, ~.x^2) %>% 
    map(~set_colnames(.x, paste0(names(.x), '_sq')))

  coll_score$expression <- 
    map2(coll_score$expression, 
         coll_score$sq_expression, 
         cbind) %>% 
    map(center_data, 'mean') %>% 
    map(rownames_to_column, 'sample_id')
    
  coll_score$data <- 
    map2(coll_score$survival, 
         coll_score$expression[names(coll_score$survival)], 
         left_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(column_to_rownames, 'sample_id')
  
  coll_score$variables <- 
    c(coll_score$variables, 
      paste0(coll_score$variables, '_sq'))

  ## survival objects
  
  coll_score$y <- coll_score$data %>% 
    map(~Surv(.x$rfs_months, .x$relapse))

  ## matrices of normalized explanatory variables
  
  coll_score$x <- coll_score$data %>% 
    map(~.x[coll_score$variables]) %>% 
    map(as.matrix)
  
  ## folds 
  
  set.seed(1234)
  
  coll_score$n_rep <- 200
  
  coll_score$folds <- 1:coll_score$n_rep %>% 
    map(function(x) createFolds(y = factor(coll_score$data$tcga$relapse), 
                                k = 10, 
                                list = FALSE, 
                                returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:coll_score$n_rep))
  
# Tuning of the lambda parameter ------
  
  insert_msg('Lambda tuning')
  
  plan('multisession')
  
  coll_score$lambda_tune <- coll_score$folds %>% 
    future_map(~cv.glmnet(x = coll_score$x$tcga, 
                          y = coll_score$y$tcga, 
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
  
  coll_score$glmnet_model <- 
    glmnet(x = coll_score$x$tcga, 
           y = coll_score$y$tcga, 
           family = 'cox', 
           alpha = 0.5, 
           lambda = coll_score$opt_lambda$lambda)
  
# Calculating the linear predictor scores for the training and test cohorts -------
  
  insert_msg('Calculating the collagen scores')
  
  ## prediction
  
  coll_score$score_tbl <- coll_score$x %>% 
    map(~predict(coll_score$glmnet_model, newx = .x)) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(set_names, c('sample_id', 'collagen_score')) %>% 
    map(as_tibble)
  
  ## appending with the survival information
  
  coll_score$score_tbl <- coll_score$data %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(select, 
        sample_id, 
        relapse, 
        rfs_months) %>% 
    map2(coll_score$score_tbl, ., 
         left_join, by = 'sample_id')

# Building univariable Cox models ------
  
  insert_msg('Uni-variable Cox models')
  
  # working with metaprogramming, to get the entire
  # data sets kept in place with the models
  
  coll_score$models <- coll_score$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(rfs_months, relapse) ~ collagen_score, 
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
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, response = ifelse(cohort == 'gse16560', 'OS', 'RFS'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  coll_score$calibration <- coll_score$models %>% 
    map(calibrate, n = 3, labels = c('low', 'int', 'high'))
  
  coll_score$global_cal <- coll_score$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(respose = ifelse(cohort == 'gse16560', 'OS', 'RFS'))

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
    mutate(response = ifelse(cohort == 'gse16560', 'OS', 'RFS')) %>% 
    re_adjust(p_variable = 'pval')

# training model estimates ------
  
  insert_msg('Estimates of the training model')
  
  coll_score$coefs <- coef(coll_score$glmnet_model) %>% 
    as.matrix %>% 
    as.data.frame %>% 
    rownames_to_column('variable') %>% 
    set_names(c('variable', 'coef')) %>% 
    mutate(exp_coef = exp(coef)) %>% 
    filter(coef != 0) %>% 
    as_tibble
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  coll_score$survival <- NULL
  coll_score$expression <- NULL
  coll_score$x <- NULL
  coll_score$y <- NULL
  coll_score$sq_expression <- NULL
  coll_score$data <- NULL
  coll_score$n_rep <- NULL
  coll_score$folds <- NULL
  
  coll_score <- compact(coll_score)
  
  save(coll_score, file = './cache/coll_score.RData')

# END -----
  
  insert_tail()