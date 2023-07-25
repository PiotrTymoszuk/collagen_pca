# Graphs for single subsystems subsystems
# I'm not plotting hyper graphs, since they are not informative 
# for the selected pathways

  insert_head()
  
# container ------
  
  meta_hg <- list()

# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
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
  
  meta_hg$subsystems <- meta_sub$common[c("int", "hi")] %>% 
    unlist %>% 
    unique
  
  meta_hg$subsystems <- 
    c(meta_hg$subsystems, c('Citric acid cycle', 'Oxidative phosphorylation'))
  
  meta_hg$reactions <- meta$regulation[[1]][[1]] %>% 
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

  ## plot title suffixes
  
  meta_hg$title_suffix <- 
    c(int = 'Collagen int vs low', 
      hi = 'Collagen high vs low')
  
  meta_hg$short_suffix <- 
    c(int = 'int vs low', 
      hi = 'high vs low')
  
# regulation estimates ------
  
  insert_msg('Regulation estimates')
  
  for(i in names(meta$regulation)) {
    
    ## labeling of the reactions: those annotated with R_*, i.e. 
    ## with no enzymes attached to them, are labeled as NA and skipped
    ## in the Forest plots of the top regulated reactions
    
    meta_hg$estimates[[i]] <- meta$regulation[[i]] %>% 
      map(mutate, 
          log_reg = log2(fold_reg), 
          log_lower = log2(lower_ci), 
          log_upper = log2(upper_ci), 
          react_name = annotate_bigg(react_id), 
          react_lab = ifelse(stri_detect(react_name, regex = '^RE\\d+'), 
                             NA, 
                             paste(react_name, react_id, sep = '\n'))) %>% 
      map(filter, 
          subsystem %in% meta_hg$subsystems, 
          !is.na(regulation)) %>% 
      map(mutate, subsystem = factor(subsystem, meta_hg$subsystems)) %>% 
      map(blast, subsystem) %>% 
      transpose
    
  }

# Numbers of total and regulated reactions -----
  
  insert_msg('Numbers of total and regulated reactions')
  
  for(i in names(meta_hg$estimates)) {
    
    meta_hg$regulation_numbers[[i]] <- 
      meta_hg$estimates[[i]] %>% 
      map(map, 
          count, 
          regulation, 
          .drop = FALSE) %>% 
      map(map, 
          ~rbind(.x, 
                 tibble(regulation = 'total', 
                        n = sum(.x$n))))
    
    ## ready-to-use plot captions for the regulation plots
    
    meta_hg$regulation_caps[[i]] <- meta_hg$regulation_numbers[[i]] %>% 
      map(map, 
          ~paste0('activated: n = ', .x$n[1], 
                  ', inhibited: n = ', .x$n[2]))
    
    meta_hg$regulation_x <- meta_hg$regulation_numbers[[i]] %>% 
      map(map, 
          ~paste0('Reactions, n = ', .x$n[4]))
    
    
  }

# Regulation plots ------
  
  insert_msg('Regulation plots')
  
  for(i in names(meta_hg$estimates)) {
    
    meta_hg$regulation_plots[[i]] <- 
      list(data = meta_hg$estimates[[i]], 
           sub = names(meta_hg$estimates[[i]]), 
           caps = meta_hg$regulation_caps[[i]], 
           x_labs = meta_hg$regulation_x) %>% 
      pmap(function(data, sub, caps, x_labs) list(data = data, 
                                                  plot_title = paste(sub, 
                                                                     globals$study_labels[names(data)], 
                                                                     sep = ', '), 
                                                  plot_subtitle = caps) %>% 
             future_pmap(safely(plot_sign), 
                         regulation_variable = 'log_reg', 
                         p_variable = 'p_adjusted', 
                         signif_level = 0.05, 
                         regulation_level = 0, 
                         cust_theme = globals$common_theme, 
                         x_lab = x_labs, 
                         y_lab = paste0(meta_hg$short_suffix[[i]], 
                                        ', log<sub>2</sub> fold-regulation'), 
                         show_trend = FALSE, 
                         .options = furrr_options(seed = TRUE)) %>% 
             map(~.x$result) %>% 
             compact %>% 
             map(~.x + 
                   scale_fill_manual(values = c(upregulated = 'firebrick', 
                                                downregulated = 'steelblue', 
                                                ns = 'gray60'), 
                                     labels = c(upregulated = 'activated', 
                                                downregulated = 'inhibited', 
                                                ns = 'ns'), 
                                     name = '') + 
                   theme(axis.title.y = element_markdown(), 
                         plot.tag = element_blank())))
    
  }
  
# Top regulated reactions per subsystem ------
  
  insert_msg('Top regulated reactions per subsystem')
  
  ## Forest plots for the annotated reactions 
  
  for(i in names(meta_hg$estimates)) {
    
    meta_hg$top_plots[[i]] <- 
      list(data = meta_hg$estimates[[i]], 
           sub = names(meta_hg$estimates[[i]])) %>% 
      pmap(function(data, sub) list(data = data %>% 
                                      map(filter, 
                                          !is.na(react_lab), 
                                          fold_reg != 1), 
                                    plot_title = paste(sub, 
                                                       globals$study_labels[names(data)], 
                                                       sep = ', ')) %>% 
             future_pmap(plot_top, 
                         regulation_variable = 'log_reg', 
                         label_variable = 'react_lab', 
                         p_variable = 'p_adjusted', 
                         lower_ci_variable = 'log_lower', 
                         upper_ci_variable = 'log_upper', 
                         top_regulated = 10, 
                         signif_level = 0.05, 
                         regulation_level = 0, 
                         cust_theme = globals$common_theme, 
                         x_lab = paste0(meta_hg$title_suffix[[i]], 
                                        ', log<sub>2</sub> fold-regulation'), 
                         plot_subtitle = 'Top regulated reactions', 
                         .options = furrr_options(seed = TRUE)) %>% 
             map(~.x + 
                   theme(axis.title.x = element_markdown(), 
                         axis.title.y = element_blank(), 
                         plot.tag = element_blank())))
    
  }

# END ------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()