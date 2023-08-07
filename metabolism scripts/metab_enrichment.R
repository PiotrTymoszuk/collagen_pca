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
    map(suba, 
        method = 'fisher', 
        signif_type = 'fdr', 
        .parallel = TRUE)
  
# Numbers of all and significantly regulated reactions ------
  
  insert_msg('Numbers of all, and activated and inhibited reactions')
  
  meta_sub$react_n <- meta_sub$test %>% 
    map(summarise, 
        n_all = sum(n), 
        n_all_total = sum(n_total), 
        .by = status)
  
# Effect sizes for the enrichment ------
  
  insert_msg('Effect sizes for the enrichment')
  
  ## effect sizes of the enrichment: odds ratio
  ## this is computed with the following formula
  ## (n(DE, subsystem)/n(subsystem))/(n(DE)/n(total))
  
  meta_sub$test <- 
    map2(meta_sub$test, 
         meta_sub$react_n, 
         left_join, 
         by = 'status') %>% 
    map(mutate, 
        OR = (n/n_total)/(n_all/n_all_total), 
        eff_size = interpret_oddsratio(OR, rules = 'cohen1988'))
  
# Identification of significantly activated and inhibited subsystems ------
  
  insert_msg('Significant results')
  
  ## pFDR < 0.05, OR > 1 and at least small effect size (Cohen1988)
  
  meta_sub$significant <- meta_sub$test %>% 
    map(filter, 
        status %in% c('activated', 'inhibited'), 
        p_adjusted < 0.05, 
        OR > 1, 
        eff_size %in% c('small', 'medium', 'large')) %>% 
    map(blast, status) %>% 
    transpose
  
# Identification of common activated and inhibited subsystems -----
  
  insert_msg('Common enriched metabolic subsystems')
  
  ## significantly enriched in at least 4 out of 5 investigated cohorts
  
  meta_sub$common_sets <- names(meta_sub$test) %>% 
    combn(m = 4, simplify = FALSE)

  for(i in names(meta_sub$significant)) {
    
    meta_sub$common[[i]] <- meta_sub$common_sets %>% 
      map(~meta_sub$significant[[i]][.x]) %>% 
      map(map, ~.x$subsystem) %>% 
      map(reduce, intersect) %>% 
      reduce(union)

  }
  
# Plotting fractions of regulated reactions, common subsystems ------
  
  insert_msg('Bubble plots with percentages of regulated reacitions')
  
  ## plotting data: enrichment and fractions of regulated reactions
  ## for common significantly enriched subsystems
  
  meta_sub$regulation_plots$data <- meta_sub$test %>% 
    map(filter, status %in% c('activated', 'inhibited')) %>% 
    map(mutate, eff_size = as.character(eff_size)) %>% 
    compress(names_to = 'cohort') %>% 
    blast(status)
  
  meta_sub$regulation_plots$data <- 
    map2_dfr(meta_sub$regulation_plots$data, 
             meta_sub$common[names(meta_sub$regulation_plots$data)], 
             ~filter(.x, subsystem %in% .y)) 

  ## coding for significance
  
  meta_sub$regulation_plots$data <- meta_sub$regulation_plots$data %>% 
    mutate(significant = ifelse(p_adjusted < 0.05 & 
                                  OR > 1 & 
                                  eff_size %in% c('small', 'medium', 'large'), 
                                'yes', 'no'), 
           regulation = ifelse(significant == 'yes', status, 'ns'), 
           cohort = factor(cohort, names(meta_sub$test)))

  ## bubble plots
  
  meta_sub$regulation_plots$plot <- meta_sub$regulation_plots$data %>% 
    ggplot(aes(x = cohort, 
               y = reorder(subsystem, OR), 
               size = OR, 
               fill = regulation)) + 
    geom_point(shape = 21) + 
    geom_text(aes(label = signif(OR, 2), 
                  x = as.numeric(cohort) + 0.25, 
                  alpha = significant, 
                  color = regulation), 
              size = 2.5, 
              hjust = 0, 
              vjust = 0.5, 
              show.legend = FALSE) +
    scale_fill_manual(values = c(activated = 'firebrick', 
                                 inhibited = 'steelblue', 
                                 ns = 'gray60'), 
                      name = 'regulation\nstatus') + 
    scale_color_manual(values = c(activated = 'firebrick', 
                                  inhibited = 'steelblue', 
                                  ns = 'black'), 
                       name = 'regulation\nstatus') + 
    scale_alpha_manual(values = c(no = 0.5, 
                                  yes = 1)) +
    scale_radius(limits = c(0, 9.5), 
                 range = c(0.5, 5), 
                 name = 'OR') + 
    scale_x_discrete(labels = globals$study_labels, 
                     limits = names(meta_sub$test)) + 
    guides(alpha = 'none', 
           color = 'none') + 
    facet_grid(status ~ ., 
               scales = 'free', 
               space = 'free') + 
    globals$common_theme + 
    theme(axis.title = element_blank()) + 
    labs(title = 'Metabolic subsystems, hi vs low')

# caching the results -------
  
  insert_msg('Caching the results')
  
  save(meta_sub, file = './cache/meta_sub.RData')

# END -------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()