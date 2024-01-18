# Survival random forests. Tuning in the out-of-bag setting

  insert_head()
  
# container -------
  
  rf_surv <- list()
  
# Modeling data -------
  
  insert_msg('Modeling data')
  
  ## only first-order terms
  
  rf_surv$variables <- globals$genes_interest$gene_symbol
  
  ## analysis tables: from the globals
  ## normalization and centering of the explanatory variables
  ## min-max normalization of the survival time: avoiding zeros
  ## by a small shift
  
  rf_surv$data <- surv_globals$data
  
  for(i in names(rf_surv$data)) {
    
    rf_surv$data[[i]][rf_surv$variables] <- 
      rf_surv$data[[i]][rf_surv$variables] %>% 
      center_data('mean')
    
    rf_surv$data[[i]]['rfs_months'] <- 
      rf_surv$data[[i]]['rfs_months'] %>% 
      min_max
    
  }
  
  rf_surv$data <- rf_surv$data %>% 
    map(mutate, rfs_months = rfs_months + 0.01) %>% 
    map(column_to_rownames, 'sample_id')
  
# Dummy survival models -------
  
  insert_msg('Dummy survival models')
  
  rf_surv$null_models <- rf_surv$data %>% 
    map(~call2(.fn = 'coxph', 
               formula = Surv(rfs_months, relapse) ~ 1, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., rf_surv$data, as_coxex)
  
# Tune grid -------
  
  insert_msg('Tune grid')
  
  rf_surv$tune_grid <- 
    expand.grid(mtry = seq(2, 27, by = 2), 
                splitrule = c('logrank', 'bs.gradient'), 
                nsplit = c(1, 2, 5), 
                nodesize = c(5, 10, 15, 20), 
                stringsAsFactors = FALSE)
  
# Construction of the tuning models --------
  
  insert_msg('construction of the tuning models')
  
  set.seed(1234)
  
  rf_surv$tuning$models <- rf_surv$tune_grid %>% 
    pmap(rfsrc, 
         formula = Surv(rfs_months, relapse) ~ ., 
         data = rf_surv$data$geo, 
         save.memory = TRUE, 
         ntree = 500)
  
  rf_surv$tuning$models <- rf_surv$tuning$models %>% 
    set_names(paste0('cond_', 1:length(rf_surv$tuning$models)))
  
  rf_surv$tune_grid$condition <- 
    paste0('cond_', 1:length(rf_surv$tuning$models))
  
# Extraction of the performance stats from the tuning models -------
  
  insert_msg('Extraction of the tuning stats')
  
  ## C-indexes for the OOB predictions: they are stored as 1 - error rate
  
  rf_surv$tuning$summary <- rf_surv$tuning$models %>% 
    map(~.x$err.rate) %>% 
    map(na.omit) %>% 
    map(as.numeric) %>% 
    map_dbl(function(x) 1 - x) %>% 
    compress(names_to = 'condition', 
             values_to = 'c_index')
  
  rf_surv$tuning$summary <- 
    left_join(rf_surv$tune_grid, 
              rf_surv$tuning$summary, 
              by = 'condition')
  
  ## the best tune
  
  rf_surv$tuning$best_tune <- rf_surv$tuning$summary %>% 
    filter(c_index == max(c_index))
  
  rf_surv$tuning$best_tune <- 
    rf_surv$tuning$best_tune[1, ]
  
# Training the GEO model --------
  
  insert_msg('Training the GEO model')
  
  set.seed(1234)
  
  rf_surv$rf_model <- 
    rfsrc(formula = Surv(rfs_months, relapse) ~ ., 
          data = rf_surv$data$geo, 
          save.memory = FALSE, 
          ntree = 500, 
          mtry = rf_surv$tuning$best_tune$mtry[1], 
          splitrule = rf_surv$tuning$best_tune$splitrule[1],
          nodesize = rf_surv$tuning$best_tune$nodesize[1])
  
# Predictions --------
  
  insert_msg('Predictions')
  
  rf_surv$predictions <- 
    list(newdata = rf_surv$data, 
         outcome = c('train', 'test', 'test')) %>% 
    pmap(predict, 
         object = rf_surv$rf_model)
  
# Prediction stats ------
  
  insert_msg('Prediction stats')
  
  ## C-indexes, extracted as above for the tuning models
  
  rf_surv$stats$c_index <- rf_surv$predictions %>% 
    map(~.x$err.rate) %>% 
    map(na.omit) %>% 
    map(as.numeric) %>% 
    map_dbl(function(x) 1 - x)
  
  ## IBS for the models
  
  rf_surv$stats$ibs_model <- rf_surv$predictions %>% 
    map(get.brier.survival, 
        cens.model = 'rfsrc') %>% 
    map_dbl(~.x$crps.std)
  
  ## IBS for the dummy models: reference
  
  rf_surv$stats$ibs_reference <- rf_surv$null_models %>% 
    map(summary, 'fit') %>% 
    map_dbl(~.x$ibs_reference[1])
  
  ## a common table
  
  rf_surv$stats <- rf_surv$stats %>% 
    as_tibble %>%
    mutate(cohort = names(rf_surv$stats[[1]]), 
           dataset = ifelse(cohort == 'geo', 'training', 'test'))
  
# Brier scores for unique time points -------
  
  insert_msg('Brier scores for the unique time points')
  
  ## Brier scores for the model
  
  rf_surv$brier_scores <- rf_surv$predictions %>% 
    map(get.brier.survival, 
        cens.model = 'rfsrc') %>% 
    map(~.x$brier.score) %>% 
    map(as_tibble) %>% 
    map(set_names, c('time', 'test'))
  
  ## Brier scores for the reference, i.e. dummy Cox model 
  ## for the training data
  
  rf_surv$ref_brier_scores <- rf_surv$null_models$geo %>% 
    surv_brier %>% 
    select(time, reference)
  
  rf_surv$brier_scores <- rf_surv$brier_scores %>% 
    map(mutate, 
        training = rf_surv$brier_scores$geo$test) %>% 
    map(left_join, rf_surv$ref_brier_scores, 
        by = 'time') %>% 
    map(as.list) %>% 
    transpose %>% 
    pmap(brier)
  
# variable importance --------
  
  insert_msg('Variable importance')
  
  ## constructing a model with permutation importance measures

  rf_surv$importance$model <- 
    rfsrc(formula = Surv(rfs_months, relapse) ~ ., 
          data = rf_surv$data$geo, 
          save.memory = FALSE, 
          ntree = 500, 
          mtry = rf_surv$tuning$best_tune$mtry[1], 
          splitrule = rf_surv$tuning$best_tune$splitrule[1],
          nodesize = rf_surv$tuning$best_tune$nodesize[1], 
          importance = 'permute')
  
  ## extracting the importance stats, i.e. deltas of C-indexes
  
  rf_surv$importance$test <- rf_surv$importance$model$importance %>% 
    compress(names_to = 'variable', 
             values_to = 'delta_c_index')
  
# Caching the results ------
  
  insert_msg('Caching the results')
  
  rf_surv$data <- NULL
  rf_surv$ref_brier_scores <- NULL
  rf_surv$null_models <- NULL
  
  rf_surv$tuning$models <- NULL
  rf_surv$tuning <- compact(rf_surv$tuning)
  
  rf_surv$importance$model <- NULL
  rf_surv$importance <- compact(rf_surv$importance)

  rf_surv <- compact(rf_surv)
  
  save(rf_surv, file = './cache/rf_surv.RData')
  
# END ------
  
  insert_tail()
  
  