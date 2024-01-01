# Correlation of expression of the collagen-related genes of interest with PSA 
# at diagnosis done with Spearman's correlation. 
#
# Significant correlations are defined by pFDR < 0.05 and at least weak effect 
# size (rho >= 0.1). Common effects are shared by at least 5 cohorts except of 
# GSE16560.

  insert_head()
  
# container ------
  
  psa <- list()
  
# Parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# expression and PSA concentrations -------
  
  insert_msg('Expression and PSA levels')
  
  psa$genes <- globals$genes_interest$gene_symbol
  
  psa$expression <- globals$study_exprs %>% 
    eval %>%
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(select, 
        sample_id, 
        all_of(psa$genes))
  
  ## PSA levels
  
  psa$clinic <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(safely(select), sample_id, psa_diagnosis) %>% 
    map(~.x$result) %>% 
    compact
  
  psa$data <- 
    map2(psa$clinic, 
         psa$expression[names(psa$clinic)], 
         left_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x)))
  
# N numbers -------
  
  insert_msg('N numbers')
  
  psa$n_numbers <- psa$data %>% 
    map_dbl(nrow)

# Spearman's correlations --------
  
  insert_msg('Testing')
  
  for(i in names(psa$data)) {
    
    psa$test[[i]] <- list(variables = map(psa$genes, c, 'psa_diagnosis')) %>% 
      future_pmap_dfr(correlate_variables, 
                      psa$data[[i]], 
                      type = 'spearman', 
                      ci = FALSE) %>% 
      re_adjust('p_value', method = 'BH')
    
  }
  
# Formatting the results --------
  
  insert_msg('Formatting the results')
  
  psa$test <- psa$test %>% 
    map(mutate, 
        correlation = ifelse(p_adjusted >= 0.05, 'ns', 
                             ifelse(estimate >= 0.1, 'positive', 
                                    ifelse(estimate <= -0.1, 'negative', 'ns'))), 
        correlation = factor(correlation, c('positive', 'negative', 'ns')))
  
# Significant effects -------
  
  insert_msg('Significant effects')
  
  ## in single cohorts
  
  psa$significant[c('positive', 'negative')] <- 
    c('positive', 'negative') %>% 
    map(function(sig) psa$test %>% 
          map(~filter(.x, correlation == sig))) %>% 
    map(map, ~.x$variable1)
  
  ## shared by at least five cohorts: none found
  
 # psa$common_significant <- psa$significant %>% 
  #  map(shared_features, m = 5)
  
# END ------
  
  psa$expression <- NULL
  psa$clinic <- NULL
  
  psa <- compact(psa)
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()