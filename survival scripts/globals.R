# COMBAT expression estimates used for multi-parameter modeling of survival
#
# Generation of a pooled GO cohort.

  insert_head()
  
# container ------
  
  surv_globals <- list()
  
# input data: COMBAT estimates -----
  
  insert_msg('COMBAT expression estimates')
  
  ## survival 
  
  surv_globals$survival <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(function(x) if('relapse' %in% names(x)) x else NULL) %>% 
    compact %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(select, sample_id, relapse, rfs_months)
  
  ## expression: the pooled GEO cohort, TCGA and DKFZ
  
  surv_globals$expression$geo <- 
    combat$adjusted_data[c("gse54460", "gse70768", "gse70769", "gse220095")] 
  
  surv_globals$expression$tcga <- combat$adjusted_data$tcga
  surv_globals$expression$dkfz <- combat$adjusted_data$dkfz
  
  ## merging with the survival data
  
  surv_globals$data$geo <- 
    map2_dfr(surv_globals$survival[names(surv_globals$expression$geo)], 
             surv_globals$expression$geo, 
             left_join, by = 'sample_id')
  
  for(i in c('tcga', 'dkfz')) {
    
    surv_globals$data[[i]] <- 
      left_join(surv_globals$survival[[i]], 
                surv_globals$expression[[i]], 
                by = 'sample_id')
    
  }
  
  surv_globals$data <- surv_globals$data %>% 
    map(~filter(.x, complete.cases(.x)))
  
# N numbers -------
  
  insert_msg('N numbers')
  
  ## n numbers: total and events
  
  surv_globals$n_numbers <- surv_globals$data %>% 
    map(~tibble(n_total = nrow(.x), 
                n_events = sum(.x$relapse))) %>% 
    compress(names_to = 'cohort')
  
  ## ready-to-use plot captions
  
  surv_globals$n_tags <- surv_globals$n_numbers %>% 
    mutate(plot_cap = paste0('total: n = ', n_total, 
                             ', events: n = ', n_events)) %>% 
    .$plot_cap %>% 
    set_names(surv_globals$n_numbers$cohort)
  
# Colors and labels --------
  
  insert_msg('Colors and labels')
  
  surv_globals$tertile_colors <- c(low = 'darkolivegreen', 
                                   int = 'steelblue', 
                                   high = 'firebrick')
  
  surv_globals$model_colors <-
    c(test = 'steelblue', 
      training = 'indianred3')
  
  surv_globals$study_labels <- 
    c(geo = 'pooled GEO', 
      tcga = 'TCGA', 
      dkfz = 'DKFZ')
  
  surv_globals$algo_labels <- 
    c(ridge = 'Ridge Cox', 
      elnet = 'Elastic Net Cox', 
      lasso = 'LASSO Cox', 
      svm = 'SVM', 
      rf = 'Random Forest', 
      gbm = 'GBM')
  
  surv_globals$algo_colors <- 
    c(ridge = 'plum4', 
      elnet = 'steelblue2', 
      lasso = 'indianred3', 
      svm = 'gray60', 
      rf = 'darkolivegreen', 
      gbm = 'orangered2')
  
  surv_globals$algo_xlabs <- 
    c(rep('Relapse-free survival, months', 3), 
      rep('Min/max scaled relapse-free survival', 3)) %>% 
    set_names(c('ridge', 'elnet', 'lasso', 'svm', 'rf', 'gbm'))
  
# CV folds used for tuning of the GLMNET models -----
  
  insert_msg('CV folds')
  
  set.seed(1234)
  
  surv_globals$n_rep <- 200
  
  surv_globals$folds <- 1:surv_globals$n_rep %>% 
    map(function(x) createFolds(y = factor(surv_globals$data$geo$relapse), 
                                k = 10, 
                                list = FALSE, 
                                returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:surv_globals$n_rep))
  
# END -----
  
  rm(i)
  
  surv_globals$survival <- NULL
  surv_globals$expression <- NULL
  
  surv_globals <- compact(surv_globals)
  
  insert_tail()