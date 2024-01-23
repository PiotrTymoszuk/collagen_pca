# Development of a multi-gene signatures of the Gleason scores (5 - 6, 7 and 8+), 
# with ordinal Random Forest. For this reason I'm re-coding the GS 
# strata into a 
#
# As shown with overall model performance stats, prediction quality is low: 
# expression of the collagen-related genes is only loosely associated with tumor
# aggressiveness gauged by the Gleason score.


  insert_head()
  
# container -------
  
  gs_sign <- list()
  
# gene expression data -------
  
  insert_msg('Gene expression data')
  
  ## variables
  
  gs_sign$variables <- globals$genes_interest$gene_symbol
  
  ## expression of the collagen-related genes
  
  gs_sign$expression <- combat$adjusted_data

  ## Gleason scores
  
  gs_sign$clinic <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(safely(select), 
        sample_id, gleason_simple) %>% 
    map(~.x$result) %>% 
    compact
  
  ## normalization
  
  gs_sign$data <- 
    map2(gs_sign$clinic, 
         gs_sign$expression[names(gs_sign$clinic)], 
         left_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x))) 

  for(i in names(gs_sign$data)) {
    
    gs_sign$data[[i]][gs_sign$variables] <- 
      gs_sign$data[[i]][gs_sign$variables] %>% 
      center_data('mean')
    
  }
  
  ## factor level renaming and complete cases
  
  gs_sign$data <- gs_sign$data %>% 
    map(mutate, 
        gleason_simple = as.numeric(gleason_simple)) %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(column_to_rownames, 'sample_id')
  
# N numbers -------
  
  insert_msg('N numbers')
  
  gs_sign$n_numbers <- gs_sign$data %>% 
    map(count, gleason_simple) %>% 
    map(transmute, 
        .outcome = gleason_simple, 
        n = n)
  
# Tuning of a Ranger RF model in the TCGA cohort --------
  
  insert_msg('Tuning')
  
  registerDoParallel(cores = 7)
  
  gs_sign$caret_model <- 
    train(gleason_simple ~ ., 
          data = gs_sign$data$tcga, 
          method = 'ranger', 
          metric = 'Rsquared', 
          tuneGrid = expand.grid(mtry = 2:27, 
                                 splitrule = c('variance', 
                                               'extratrees', 
                                               'maxstat'), 
                                 min.node.size = c(1, 3, 5, 10)), 
          trControl = trainControl(method = 'cv', 
                                   number = 10, 
                                   savePredictions = 'final', 
                                   returnData = TRUE, 
                                   returnResamp = 'final'), 
          num.trees = 1000) %>% 
    as_caretx
  
  stopImplicitCluster()
  
# Predictions --------
  
  insert_msg('Predictions')
  
  gs_sign$predictions <- gs_sign$data %>% 
    map(~predict(gs_sign$caret_model, newdata = .x)) %>% 
    map(~.x$test)
  
# Overall prediction performance -------
  
  insert_msg('Overall prediction performance')
  
  gs_sign$stats <- gs_sign$predictions %>% 
    map(summary, wide = TRUE) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(type = ifelse(cohort == 'tcga', 
                         'training', 'test'), 
           type = factor(type, c('training', 'test')), 
           cohort = factor(cohort, names(gs_sign$predictions)))
  
# Variable importance -------
  
  insert_msg('Variable importance')
  
  ## done only at request, the models perform quite poor
  ## so it does not make much sense to present variable importances
  
  ## non-zero coeffcients corresponding to the best lambda value

  #gs_sign$coefs <- gs_sign$caret_model$finalModel %>% 
   # coef %>% 
    #compress(names_to = 'parameter', 
     #        values_to = 'coef') %>% 
  #  mutate(variable = stri_replace(parameter, regex = ':.*', replacement = ''), 
   #        level = stri_extract(parameter, regex = ':.*'), 
    #       level = stri_extract(level, regex = '\\d+'), 
     #      level = ifelse(is.na(level), 0, as.numeric(level)), 
      #     level = levels(gs_sign$data[[1]]$gleason_simple)[level + 1]) %>% 
  #  filter(coef != 0)

# Caching the results --------
  
  insert_msg('Caching the results')
  
  gs_sign <- 
    gs_sign[c("variables", 
              "caret_model", "predictions", 
              "stats", "strata_stats", "coefs")]
  
  save(gs_sign, file = './cache/gs_sign.RData')
  
# END --------
  
  rm(i)
  
  insert_tail()
  