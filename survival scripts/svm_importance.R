# Permutation importance of explanatory variables in a survival SVM score.
# The importance is measured by comparing the CV out-of-fold C-indexes between 
# the full model and models with particular explanatory variables re-shuffled 
# at random

  insert_head()
  
# container ------
  
  svm_imp <- list()
  
# analysis globals --------
  
  insert_msg('Analysis globals')
  
  ## analysis data set for the training pooled GEO cohort
  ## and tuned parameters of the SVM model
  ## CV folds used also for tuning of the genuine model
  
  svm_imp$data <- svm_surv$data$geo
  
  svm_imp$best_tune <- svm_surv$tuning$best_tune
  
  svm_imp$folds <- svm_surv$folds
  
# Importance testing -------
  
  insert_msg('Importance testing')
  
  plan('multisession')
  
  svm_imp$test <- svm_importance(data = svm_imp$data, 
                                 time_variable = 'rfs_months', 
                                 event_variable = 'relapse', 
                                 folds = svm_imp$folds, 
                                 type = svm_imp$best_tune$type[[1]], 
                                 diff.meth = svm_imp$best_tune$diff.meth[[1]], 
                                 gamma.mu = svm_imp$best_tune$gamma.mu[[1]], 
                                 kernel = svm_imp$best_tune$kernel[[1]])
  
  plan('sequential')
  
# Caching the results --------
  
  insert_msg('Caching the results')
  
  svm_imp <- svm_imp[c("best_tune", "test")]
  
  save(svm_imp, file = './cache/svm_imp.RData')
  
# END ------