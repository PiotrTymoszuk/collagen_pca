# Differences in relapse-free survival between the collagen clusters 
# investigated by the Peto-Peto test

  insert_head()
  
# container ------
  
  ana_rfs <- list()
  
# analysis data -------
  
  insert_msg('Analysis data')
  
  ## data set labels
  
  ana_rfs$study_labels <- 
    c(geo = 'pooled GEO', 
      globals$study_labels[c("tcga", "dkfz")])
  
  ## survival
  
  ana_rfs$data <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(safely(select), sample_id, relapse, rfs_months) %>% 
    map(~.x$result) %>% 
    compact
  
  ## cluster assignment
  
  ana_rfs$data <- 
    map2(ana_rfs$data, 
         ana_globals$assignment[names(ana_rfs$data)], 
         left_join, by = 'sample_id') %>% 
    map(filter, 
        !is.na(relapse),
        !is.na(rfs_months))
  
  ## the pooled GEO cohort: as discussed in the manuscript text.
  ## single GEO cohorts are likely underpowered
  
  ana_rfs$data$geo <- 
    ana_rfs$data[c("gse54460", "gse70768", "gse70769", "gse220095")] %>% 
    compress(names_to = 'cohort')
  
  ana_rfs$data <- ana_rfs$data[c("geo", "tcga", "dkfz")]
  
# n numbers ------
  
  insert_msg('N numbers')
  
  ## in the clusters
  
  ana_rfs$n_numbers <- ana_rfs$data %>% 
    map(count, clust_id)
  
  ana_rfs$n_labs <- ana_rfs$n_numbers %>% 
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = '))
  
  ## total and events
  
  ana_rfs$plot_caps <- ana_rfs$data %>% 
    map(~paste0('total: n = ', nrow(.x), 
                ', events: n = ', sum(.x$relapse)))
  
  
# Survfit objects ------
  
  insert_msg('Surv fit objects')
  
  ana_rfs$survfit_obj <- ana_rfs$data %>% 
    map(survminer::surv_fit, 
        formula = Surv(rfs_months, relapse) ~ clust_id)
  
# Median survival and testing ------
  
  insert_msg('Survival stats and tests')
  
  ana_rfs$stats <- ana_rfs$survfit_obj %>% 
    map(surv_median) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(clust_id = stri_extract(strata, regex = 'Collag.*'), 
           clust_id = factor(clust_id, levels(ana_rfs$data[[1]]$clust_id))) %>% 
    as_tibble
  
  ana_rfs$test <- ana_rfs$survfit_obj %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    re_adjust('pval') %>% 
    as_tibble
  
# Kaplan-Meier plots -------
  
  insert_msg('Kaplan-Meier plots')
  
  ana_rfs$plots <- 
    list(fit = ana_rfs$survfit_obj, 
         title = ana_rfs$study_labels[names(ana_rfs$survfit_obj)], 
         legend.labs = ana_rfs$n_labs, 
         pval = ana_rfs$test$significance) %>% 
    pmap(ggsurvplot, 
         palette = unname(globals$cluster_colors[c("Collagen low", "Collagen hi")]), 
         legend.title = '', 
         xlab = 'Relapse-free survival, months', 
         pval.size = 2.75) %>% 
    map(~.x$plot)
  
  ## styling and number of events
  
  ana_rfs$plots <- 
    list(x = ana_rfs$plots, 
         y = ana_rfs$plot_caps) %>% 
    pmap(function(x, y) x + 
           labs(subtitle = y) + 
           globals$common_theme)
  
# END -------
  
  ana_rfs <- ana_rfs[c("survfit_obj", "stats", "test", "plots")]
  
  insert_tail()