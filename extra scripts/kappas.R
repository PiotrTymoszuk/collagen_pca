# Cohen's kappa and other ROC stats for the proteome models

  insert_head()
  
# container -----
  
  extra_kappa <- list()
  
# parallel backend --------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# data ------
  
  insert_msg('Data: confusion matrices')
  
  ## input per hand based on the sensitivity and specificity measures
  ## as well as n numbers provided for the proteome SVM models

  ## detection of any PCa with the tissue proteome
  ## whole data set n numbers
  
  extra_kappa$data$svm_tissue_pca <- 
    data_from_sesp(sens = 0.83, 
                   spec = 1, 
                   n_outcome = 90, 
                   n_total = 104)
  
  ## discrimination between ISUP2+ and ISUP1
  ## whole data set n numbers
  
  extra_kappa$data$svm_tissue_isup <- 
    data_from_sesp(sens = 0.8, 
                   spec = 0.88, 
                   n_outcome = 33, 
                   n_total = 57)
  
  ## urinome: detection of any PCa
  ## test subset n numbers
  
  extra_kappa$data$svm_urine_pca <-
    data_from_ppnp(ppv = 85.4/100, 
                   npv = 84.2/100, 
                   n_outcome = 138, 
                   n_total = 138 + 272)
  
  ## urinome: detection of ISUP2+ lesions vs benign PCa
  
  extra_kappa$data$svm_urine_isup <-
    data_from_ppnp(ppv = 58.4/100, 
                   npv = 92.8/100, 
                   n_outcome = 26, 
                   n_total = 26 + 272)
  
  ## urinome: combi models with collagen-related peptides
  ## age and PSA
  
  extra_kappa$data$age_psa_collagen_urine_pca <-
    data_from_ppnp(ppv = 87/100, 
                   npv = 70/100, 
                   n_outcome = 138, 
                   n_total = 138 + 272)
  
  extra_kappa$data$age_psa_collagen_urine_isup <-
    data_from_ppnp(ppv = 49.3/100, 
                   npv = 92.6/100, 
                   n_outcome = 26, 
                   n_total = 26 + 272)
  

# ROC stats ------  

  insert_msg('ROC stats')
  
  ## basic stats
  
  extra_kappa$stats <- extra_kappa$data %>% 
    map(roc_stats) %>% 
    compress(names_to = 'model') %>% 
    relocate(model)
  
  ## confidence intervals for Kappa, sensitivity and specificity
  
  extra_kappa[c('kappa_ci', 'sensitivity_ci', 'specificity_ci')] <- 
    c('Kappa', 'Sensitivity', 'Specificity') %>% 
    map(function(stat_name) extra_kappa$data %>% 
          future_map(~bmap(.x, 
                           FUN = function(x) roc_stats(x)[[stat_name]], 
                           B = 1000, 
                           ci_method = 'bca'), 
                     .options = furrr_options(seed = TRUE)))
  
  extra_kappa[c('kappa_ci', 'sensitivity_ci', 'specificity_ci')] <- 
    extra_kappa[c('kappa_ci', 'sensitivity_ci', 'specificity_ci')] %>% 
    map(function(boot_res) boot_res %>% 
          map(~tibble(statistic_value = .x[['dataset']], 
                      lower_ci = .x[['bootstrap']]$boot_lower_ci, 
                      upper_ci = .x[['bootstrap']]$boot_upper_ci)) %>% 
          compress(names_to = 'model') %>% 
          relocate(names_to = 'model'))
  
# Saving the table as an Excel file ------
  
  insert_msg('Saving an Excel file')
  
  extra_kappa[c("stats", "kappa_ci", "sensitivity_ci", "specificity_ci")] %>% 
    write_xlsx(path = './report/manuscript extras/svm_stats.xlsx')

# END ------
  
  plan('sequential')
  
  insert_tail()