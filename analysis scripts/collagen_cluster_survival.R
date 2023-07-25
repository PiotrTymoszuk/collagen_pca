# Differences in survival between the collagen clusters

  insert_head()
  
# container -----
  
  clust_surv <- list()
  
# globals -----
  
  insert_msg('Globals setup')
  
  ## analysis tables
  
  clust_surv$analysis_tbl <- study_data %>% 
    map(~.x$expression) %>% 
    map(extract_tumor_samples) %>% 
    map(select, 
        patient_id, 
        any_of(c('death', 'vitality_fup')), 
        any_of(c('relapse', 'relapse_fup'))) %>% 
    map2(coll_clust$assignment, ., 
         left_join, by = 'patient_id') %>% 
    map(~filter(.x, complete.cases(.x)))
  
  ## model formulas
  
  clust_surv$formulas <- list(GSE16560 = Surv(vitality_fup, death) ~ clust_id, 
                              GSE40272 = Surv(relapse_fup, relapse) ~ clust_id, 
                              GSE70768 = Surv(relapse_fup, relapse) ~ clust_id, 
                              GSE70769 = Surv(relapse_fup, relapse) ~ clust_id,
                              tcga = Surv(relapse_fup, relapse) ~ clust_id)
  
  ## n numbers, clusters
  
  clust_surv$clust_n <- clust_surv$analysis_tbl %>% 
    map(count, clust_id)
  
  clust_surv$clust_caps <- clust_surv$clust_n %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = ': n = ') %>% 
          set_names(levels(clust_surv$analysis_tbl[[1]]$clust_id)))
  
  ## n numbers: total and events
  
  clust_surv$n_numbers <- 
    map2(clust_surv$analysis_tbl, 
         c('death', rep('relapse', 4)), 
         ~count(.x, .data[[.y]])) %>% 
    map(~c(total = sum(.x$n), 
           events = .x$n[2]))

  clust_surv$n_tags <- clust_surv$n_numbers %>% 
    map(~paste0('total: n = ', .x['total'], 
                ', events: n = ', .x['events']))
  
# Testing with log-rank tests ------
  
  insert_msg('Log-rank tests')
  
  ## surv-fit objects: overall and pairwise
  
  clust_surv$survfit_obj$overall <- 
    map2(clust_surv$formulas, 
         clust_surv$analysis_tbl, 
         survminer::surv_fit)
  
  clust_surv$survfit_obj$hi_lo <- 
    map2(clust_surv$formulas, 
         clust_surv$analysis_tbl %>% 
           map(filter, clust_id %in% c('Collagen hi', 'Collagen low')) %>% 
           map(mutate, clust_id = droplevels(clust_id)), 
         survminer::surv_fit)

  ## test
  
  clust_surv$test <- clust_surv$survfit_obj %>% 
    map(map, surv_pvalue) %>% 
    map(compress, names_to = 'cohort') %>% 
    compress(names_to = 'comparison') %>% 
    mutate(comparison = factor(comparison, c('overall', 'hi_lo'))) %>% 
    re_adjust('pval') %>% 
    blast(comparison)
  
  ## KM plot labels
  
  clust_surv$km_p_labs <- 
    map2(clust_surv$test$overall$significance, 
         clust_surv$test$hi_lo$significance, 
         ~paste0('overall: ', .x, 
                 '\nhi vs low: ', .y))

# Kaplan-Meier plots -------
  
  insert_msg('KM plots')

  clust_surv$km_plots <- list(fit = clust_surv$survfit_obj$overall, 
                              pval = clust_surv$km_p_labs, 
                              title = paste('Collagen clusters and survival,', 
                                            globals$study_labels[names(clust_surv$analysis_tbl)]), 
                              xlab = c('Overall survival, months', 
                                       rep('Relapse-free survival, months', 4)),  
                              legend.labs = clust_surv$clust_caps) %>% 
    pmap(ggsurvplot, 
         palette = unname(globals$cluster_colors), 
         ggtheme = globals$common_theme, 
         legend = 'right', 
         legend.title = 'Cluster', 
         pval.size = 2.75) %>% 
    map2(., clust_surv$n_tags, 
         ~.x$plot + 
           labs(subtitle = .y))
  
# END ------
  
  insert_tail()