# Identification of the significantly regulated metabolic subsystems
# Enrichment analysis with Fisher's exact test

  insert_head()
  
# container ------
  
  meta_sub <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# enrichment testing -------
  
  insert_msg('Enrichment testing')
  
  meta_sub$test <- meta$models %>% 
    map(map, 
        suba, 
        method = 'fisher', 
        signif_type = 'fdr', 
        .parallel = TRUE)
  
# Numbers of all and significantly regulated reactions ------
  
  insert_msg('Numbers of all and significantly regulated reactions')
  
  meta_sub$react_n <- meta_sub$test %>% 
    map(map, 
        summarise, 
        n_all = sum(n), 
        n_all_total = sum(n_total), 
        .by = status)
  
# Effect sizes for the enrichment ------
  
  insert_msg('Effect sizes for the enrichment')
  
  ## effect sizes of the enrichment: odds ratio
  ## this is computed with the following formula
  ## (n(DE, subsystem)/n(subsystem))/(n(DE)/n(total))
  
  for(i in names(meta_sub$test)) {
    
    meta_sub$test[[i]] <- 
      map2(meta_sub$test[[i]], 
           meta_sub$react_n[[i]], 
           left_join, by = 'status') %>% 
      map(mutate, 
          OR = (n/n_total)/(n_all/n_all_total), 
          eff_size = interpret_oddsratio(OR, rules = 'cohen1988'))
    
  }
  
# Identification of significantly activated and inhibited subsystems ------
  
  insert_msg('Significant results')
  
  ## pFDR < 0.05, OR > 1 and at least small effect size (Cohen1988)
  
  meta_sub$significant <- meta_sub$test %>% 
    map(map, 
        filter, 
        status %in% c('activated', 'inhibited'), 
        p_adjusted < 0.05, 
        OR > 1, 
        eff_size %in% c('small', 'medium', 'large')) %>% 
    map(map, blast, status) %>% 
    map(transpose)
  
# Identification of common activated and inhibited subsystems -----
  
  insert_msg('Common enriched metabolic subsystems')
  
  ## significantly enriched in at least 4 out of 5 investigated cohorts
  
  meta_sub$common_sets <- names(meta_sub$significant[[1]][[1]]) %>% 
    combn(m = 4, simplify = FALSE)
  
  ## common activated and inhibited pathways
  ## all cohorts except of GSE40272
  
  for(i in names(meta_sub$significant)) {
    
    meta_sub$common[[i]] <- meta_sub$significant[[i]] %>% 
      map(function(data) meta_sub$common_sets %>% 
            map(~data[.x]) %>% 
            map(map, ~.x$subsystem) %>% 
            map(reduce, intersect)) %>% 
      map(reduce, union)

  }
  
# Plotting fractions of regulated reactions, common subsystems ------
  
  insert_msg('Bubble plots with percentages of regulated reacitions')
  
  ## plotting data: enrichment and fractions of regulated reactions
  ## for common significantly enriched subsystems
  
  meta_sub$regulation_plots$data <- meta_sub$test %>% 
    map(map, filter, status %in% c('activated', 'inhibited')) %>% 
    map(map, mutate, eff_size = as.character(eff_size)) %>% 
    map(compress, names_to = 'cohort') %>% 
    map(blast, status)
  
  for(i in names(meta_sub$regulation_plots$data)) {
    
    meta_sub$regulation_plots$data[[i]] <- 
      map2_dfr(meta_sub$regulation_plots$data[[i]], 
               meta_sub$common[[i]], 
               ~filter(.x, subsystem %in% .y))
    
  }
  
  ## coding for significance
  
  meta_sub$regulation_plots$data <- meta_sub$regulation_plots$data %>% 
    map(mutate, 
        significant = ifelse(p_adjusted < 0.05 & 
                               OR > 1 & 
                               eff_size %in% c('small', 'medium', 'large'), 
                             'yes', 'no'), 
        regulation = ifelse(significant == 'yes', status, 'ns'))
  
  ## bubble plots
  
  meta_sub$regulation_plots$plots <- 
    list(x = meta_sub$regulation_plots$data, 
         y = c('Common metabolic sybsystems, Collagen int vs low', 
               'Common metabolic sybsystems, Collagen high vs low')) %>% 
    pmap(function(x, y) x %>% 
           ggplot(aes(x = cohort, 
                      y = reorder(subsystem, OR), 
                      size = OR, 
                      fill = regulation)) + 
           geom_point(shape = 21) + 
           geom_text(aes(label = signif(OR, 2), 
                         alpha = significant, 
                         color = regulation), 
                     size = 2.75, 
                     hjust = -1.4, 
                     vjust = 0.5) +
           scale_fill_manual(values = c(activated = 'firebrick', 
                                        inhibited = 'steelblue', 
                                        ns = 'gray60'), 
                             name = 'regulation\nstatus') + 
           scale_color_manual(values = c(activated = 'firebrick', 
                                         inhibited = 'steelblue', 
                                         ns = 'black'), 
                              name = 'regulation\nstatus') + 
           scale_alpha_manual(values = c(no = 0.25, 
                                         yes = 1)) +
           scale_x_discrete(labels = globals$study_labels) + 
           scale_y_discrete(labels = function(x) stri_replace_last(x, 
                                                                   fixed = ' and', 
                                                                   replacement = '\nand')) + 
           guides(alpha = 'none', 
                  color = 'none') + 
           facet_grid(status ~ ., 
                      scales = 'free', 
                      space = 'free') + 
           globals$common_theme + 
           theme(axis.title = element_blank()) + 
           labs(title = y))
  
# caching the results -------
  
  insert_msg('Caching the results')
  
  save(meta_sub, file = './cache/meta_sub.RData')

# END -------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()