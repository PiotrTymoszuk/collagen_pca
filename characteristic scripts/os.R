# Differences in overall survival between the collagen clusters investigated 
# by the Peto-Peto test

  insert_head()
  
# container ------
  
  ana_os <- list()
  
# analysis data -------
  
  insert_msg('Analysis data')
  
  ## survival
  
  ana_os$data <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(safely(select), sample_id, death, os_months) %>% 
    map(~.x$result) %>% 
    compact
  
  ## cluster assignment
  
  ana_os$data <- 
    map2(ana_os$data, 
         ana_globals$assignment[names(ana_os$data)], 
         left_join, by = 'sample_id')
  
# n numbers ------
  
  insert_msg('N numbers')
  
  ## in the clusters
  
  ana_os$n_numbers <- ana_os$data %>% 
    map(count, clust_id)
  
  ana_os$n_labs <- ana_os$n_numbers %>% 
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = '))
  
  ## total and events
  
  ana_os$plot_caps <- ana_os$data %>% 
    map(~paste0('total: n = ', nrow(.x), 
                ', events: n = ', sum(.x$death)))
  
  
# Survfit objects ------
  
  insert_msg('Surv fit objects')
  
  ana_os$survfit_obj <- ana_os$data %>% 
    map(survminer::surv_fit, 
        formula = Surv(os_months, death) ~ clust_id)
  
# Median survival and testing ------
  
  insert_msg('Survival stats and tests')
  
  ana_os$stats <- ana_os$survfit_obj %>% 
    map(surv_median) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(clust_id = stri_extract(strata, regex = 'Collag.*'), 
           clust_id = factor(clust_id, levels(ana_os$data[[1]]$clust_id))) %>% 
    as_tibble
  
  ana_os$test <- ana_os$survfit_obj %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    re_adjust('pval') %>% 
    as_tibble
  
# Kaplan-Meier plots -------
  
  insert_msg('Kaplan-Meier plots')
  
  ana_os$plots <- 
    list(fit = ana_os$survfit_obj, 
         title = globals$study_labels[names(ana_os$survfit_obj)], 
         legend.labs = ana_os$n_labs, 
         pval = ana_os$test$significance) %>% 
    pmap(ggsurvplot, 
         palette = unname(globals$cluster_colors[c("Collagen low", "Collagen hi")]), 
         legend.title = '', 
         xlab = 'Overall survival, months', 
         pval.size = 2.75) %>% 
    map(~.x$plot)
  
  ## styling and number of events
  
  ana_os$plots <- 
    list(x = ana_os$plots, 
         y = ana_os$plot_caps) %>% 
    pmap(function(x, y) x + 
           labs(subtitle = y) + 
           globals$common_theme)
  
# END -------
  
  ana_os <- ana_os[c("survfit_obj", "stats", "test", "plots")]
  
  insert_tail()