# Graphs for single subsystems subsystems
# I'm not plotting hyper graphs, since they are not informative 
# for the selected pathways

  insert_head()
  
# container ------
  
  meta_hg <- list()

# Globals -------
  
  insert_msg('Globals setup')
  
  ## regex to exclude irrelevant compounds
  
  meta_hg$excl_regex <- 
    paste0('(^o2_)|(^h_)|(^h2o_)|(^h2o2_)|(^nadph_)|(^nadp_)|(^nh4_)|', 
           '(^accoa_)|(^coa_)|(^co2_)|(^nad_)|(^nadh_)|(^atp_)|(^adp_)|', 
           '(^hco3)|(^fad_)|(^fadh2_)|(^amp_)|(^ppi_)|(^pi_)|(^ac_)|', 
           '(^gdp_)|(^gtp_)')

  ## subsystems of interest, reactions and metabolites
  ## the common ones, citric acid cycle and OxPhos
  
  meta_hg$subsystems <- meta_sub$common[c("activated", "inhibited")] %>% 
    unlist %>% 
    unique
  
  meta_hg$subsystems <- 
    c(meta_hg$subsystems, c('Oxidative phosphorylation'))
  
  meta_hg$reactions <- meta$regulation[[1]] %>% 
    filter(subsystem %in% meta_hg$subsystems) %>% 
    mutate(subsystem = factor(subsystem, meta_hg$subsystems)) %>% 
    blast(subsystem) %>% 
    map(~.$react_id)
  
  meta_hg$metabolites <- meta_hg$reactions %>% 
    map(react_to_metab, 
        annotation_db = reactions, 
        exc_regex = meta_hg$excl_regex) %>% 
    map(~.x$metab_id) %>% 
    map(unlist) %>% 
    map(unique)

# regulation estimates ------
  
  insert_msg('Regulation estimates')
  
  meta_hg$estimates <- meta$regulation %>% 
    map(filter, 
        subsystem %in% meta_hg$subsystems, 
        !is.na(regulation)) %>% 
    map(mutate, 
        react_name = annotate_bigg(react_id), 
        react_lab = ifelse(stri_detect(react_name, regex = '^RE\\d+'), 
                           NA, 
                           paste(react_name, react_id, sep = '\n'))) %>% 
    map(mutate, subsystem = factor(subsystem, meta_hg$subsystems)) %>% 
    map(blast, subsystem) %>% 
    transpose
  
# Numbers of total and regulated reactions -----

  insert_msg('Numbers of total and regulated reactions')
  
  for(i in names(meta_hg$estimates)) {
    
    meta_hg$regulation_numbers[[i]] <- 
      meta_hg$estimates[[i]] %>% 
      map(count, 
          regulation, 
          .drop = FALSE) %>% 
      map(~rbind(.x, 
                 tibble(regulation = 'total', 
                        n = sum(.x$n))))
    
    ## ready-to-use plot captions for the regulation plots
    
    meta_hg$regulation_caps[[i]] <- meta_hg$regulation_numbers[[i]] %>% 
      map(~paste0('activated: n = ', .x$n[1], 
                  ', inhibited: n = ', .x$n[2]))
    
    meta_hg$regulation_x[[i]] <- meta_hg$regulation_numbers[[i]] %>% 
      map(~paste0('Reactions, n = ', .x$n[4]))
    
  }

# Regulation plots ------
  
  insert_msg('Regulation plots')
  
  for(i in names(meta_hg$estimates)) {
    
    meta_hg$regulation_plots[[i]] <- 
      list(data = meta_hg$estimates[[i]], 
           plot_title = paste(i, 
                              globals$study_labels[names(meta_hg$estimates[[i]])], 
                              sep = ', '), 
           plot_subtitle = meta_hg$regulation_caps[[i]], 
           x_lab = meta_hg$regulation_x[[i]]) %>% 
      pmap(plot_sign, 
           regulation_variable = 'log_fold_reg', 
           p_variable = 'p_adjusted', 
           signif_level = 0.05, 
           regulation_level = 0, 
           cust_theme = globals$common_theme, 
           y_lab = expression('log'[2] * ' regulation, hi vs low'), 
           show_trend = FALSE) %>% 
      map(~.x + 
            scale_fill_manual(values = c(upregulated = 'firebrick', 
                                         downregulated = 'steelblue', 
                                         ns = 'gray60'), 
                              labels = c(upregulated = 'activated', 
                                         downregulated = 'inhibited', 
                                         ns = 'ns'), 
                              name = '') + 
            theme(plot.tag = element_blank()))
 
  }
  
# Top regulated reactions per subsystem ------
  
  insert_msg('Top regulated reactions per subsystem')
  
  ## Forest plots for the annotated reactions 
  
  for(i in names(meta_hg$estimates)) {
    
    meta_hg$top_plots[[i]] <- 
      list(data = meta_hg$estimates[[i]] %>% 
             map(filter, 
                 !is.na(react_lab), 
                 regulation != 'ns'), 
           plot_title = paste(i, 
                              globals$study_labels[names(meta_hg$estimates[[i]])], 
                              sep = ', ')) %>% 
      pmap(plot_top, 
           regulation_variable = 'log_fold_reg', 
           label_variable = 'react_lab', 
           p_variable = 'p_adjusted', 
           lower_ci_variable = 'log_lower_ci', 
           upper_ci_variable = 'log_upper_ci', 
           top_regulated = 10, 
           signif_level = 0.05, 
           regulation_level = 0, 
           cust_theme = globals$common_theme, 
           x_lab = expression('log'[2] * ' regulation'), 
           plot_subtitle = 'Top regulated reactions') %>% 
      map(~.x + 
            scale_color_manual(values = c(upregulated = 'firebrick', 
                                          downregulated = 'steelblue', 
                                          ns = 'gray60'), 
                               labels = c(upregulated = 'activated', 
                                          downregulated = 'inhibited', 
                                          ns = 'ns')) + 
            scale_color_manual(values = c(upregulated = 'firebrick', 
                                          downregulated = 'steelblue', 
                                          ns = 'gray60'), 
                               labels = c(upregulated = 'activated', 
                                          downregulated = 'inhibited', 
                                          ns = 'ns')) + 
            theme(axis.title.y = element_blank(), 
                  plot.tag = element_blank()))

  }

# END ------
  
  #rm(i)
  
  insert_tail()