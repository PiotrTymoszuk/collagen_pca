# Differences in xCell and MCP counter infiltration 
# between the Collagen Clusters

  insert_head()
  
# container ------
  
  clust_infil <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals -------
  
  insert_msg('Globals')

  ## infiltration variables
  
  clust_infil$variables[c('xcell', 'mcp')] <- 
    infil[c("xcell_types", "mcp_counter_types")]
  
  ## analysis tables
  
  clust_infil$assignment <- coll_clust$assignment %>% 
    map(mutate, 
        clust_id = stri_extract(clust_id, regex = 'low|int|hi'), 
        clust_id = factor(clust_id, c('low', 'int', 'hi')), 
        clust_id = droplevels(clust_id))
  
  clust_infil$analysis_tbl$mcp <- clust_infil$assignment
  clust_infil$analysis_tbl$xcell <- clust_infil$analysis_tbl$mcp
  
  clust_infil$analysis_tbl$xcell <- 
    map2(clust_infil$analysis_tbl$xcell, 
         infil$xcell, 
         left_join, by = 'patient_id')
  
  clust_infil$analysis_tbl$mcp <- 
    map2(clust_infil$analysis_tbl$mcp, 
         infil$mcp_counter, 
         left_join, by = 'patient_id')
  
  ## n numbers 
  
  clust_infil$n_numbers <- clust_infil$analysis_tbl %>% 
    map(map, count, clust_id)
  
  for(i in names(clust_infil$n_numbers)) {
    
    clust_infil$n_caps[[i]] <- clust_infil$n_numbers[[i]] %>% 
      map(~map2_chr(.x[[1]], .x[[2]], 
                    paste, sep = ': n = ')) %>% 
      map(paste, collapse = ', ')
    
    
  }
  
# Descriptive stats ------
  
  insert_msg('Descriptive stats')
  
  for(i in names(clust_infil$analysis_tbl)) {
    
    clust_infil$stats[[i]] <- clust_infil$analysis_tbl[[i]] %>% 
      future_map(~explore(.x, 
                          variables = clust_infil$variables[[i]][clust_infil$variables[[i]] %in% names(.x)], 
                          split_factor = 'clust_id', 
                          what = 'table', 
                          pub_styled = TRUE), 
                 .options = furrr_options(seed = TRUE)) %>% 
      map(reduce, left_join, by = 'variable') %>% 
      map(set_names, 
          c('variable', levels(clust_infil$assignment[[1]]$clust_id)))

  }
  
# Mann-Whitney test ------
  
  insert_msg('Mann-whitney test')
  
  for(i in names(clust_infil$analysis_tbl)) {
    
    clust_infil$test[[i]] <- clust_infil$analysis_tbl[[i]] %>% 
      future_map(~compare_variables(.x, 
                                    variables = clust_infil$variables[[i]][clust_infil$variables[[i]] %in% names(.x)], 
                                    split_factor = 'clust_id', 
                                    what = 'eff_size',
                                    types = 'wilcoxon_r', 
                                    ci = FALSE, 
                                    exact = FALSE, 
                                    pub_styled = TRUE, 
                                    adj_method = 'BH'),
                 .options = furrr_options(seed = TRUE)) %>% 
      map(mutate, 
          plot_cap = paste(eff_size, significance, sep = ', '), 
          plot_lab = paste(variable, plot_cap, sep = '<br>'), 
          plot_lab = ifelse(p_adjusted < 0.05, 
                            paste0('<b>', plot_lab, '</b>'), 
                            plot_lab))
    
    ## significant effects
    
    clust_infil$significant[[i]] <- clust_infil$test[[i]] %>% 
      map(filter, p_adjusted < 0.05) %>% 
      map(~.x$variable)
    
  }
  
# Plots for infiltration variables, MCP counter ------
  
  insert_msg('Plots, MCP counter')
  
  for(i in names(clust_infil$analysis_tbl$mcp)) {
    
    clust_infil$plots$mcp[[i]] <- 
      list(variable = clust_infil$test$mcp[[i]]$variable, 
           plot_title = paste(clust_infil$test$mcp[[i]]$variable, 
                              globals$study_labels[i], 
                              sep = ', '), 
           plot_subtitle = clust_infil$test$mcp[[i]]$plot_cap) %>% 
      future_pmap(safely(plot_variable), 
                  clust_infil$analysis_tbl$mcp[[i]], 
                  split_factor = 'clust_id', 
                  type = 'violin', 
                  cust_theme = globals$common_theme, 
                  y_lab = 'infiltration estimate, MCP counter', 
                  x_n_labs = TRUE, 
                  point_hjitter = 0, 
                  .options = furrr_options(seed = TRUE)) %>% 
      set_names(clust_infil$test$mcp[[i]]$variable) %>% 
      map(~.x$result) %>% 
      map(~.x + 
            scale_fill_manual(values = set_names(globals$cluster_colors, 
                                                 c('hi', 'int', 'low'))))
    
  }
  
# Plots for infiltration variables, xCell ------
  
  insert_msg('Plots, xCell')
  
  for(i in names(clust_infil$analysis_tbl$xcell)) {
    
    clust_infil$plots$xcell[[i]] <- 
      list(variable = clust_infil$test$xcell[[i]]$variable, 
           plot_title = paste(clust_infil$test$xcell[[i]]$variable, 
                              globals$study_labels[i], 
                              sep = ', '), 
           plot_subtitle = clust_infil$test$xcell[[i]]$plot_cap) %>% 
      future_pmap(safely(plot_variable), 
                  clust_infil$analysis_tbl$xcell[[i]], 
                  split_factor = 'clust_id', 
                  type = 'violin', 
                  cust_theme = globals$common_theme, 
                  y_lab = 'infiltration estimate, xCell', 
                  x_n_labs = TRUE, 
                  point_hjitter = 0, 
                  .options = furrr_options(seed = TRUE)) %>% 
      set_names(clust_infil$test$xcell[[i]]$variable) %>% 
      map(~.x$result) %>% 
      map(~.x + 
            scale_fill_manual(values = set_names(globals$cluster_colors, 
                                                 c('hi', 'int', 'low'))))
    
  }
  
# Box panels with the infiltration estimates, MCP counter -------
  
  insert_msg('Box plot panels, MCP counter')
  
  ## I'm displaying Z-scores of the infiltration estimates
  
  for(i in names(clust_infil$analysis_tbl$mcp)) {
    
    clust_infil$box_panels$mcp[[i]] <- 
      draw_violin_panel(data = clust_infil$analysis_tbl$mcp[[i]] %>% 
                          map_dfc(function(x) if(is.numeric(x)) scale(x)[, 1] else x), 
                        variables = clust_infil$variables$mcp[clust_infil$variables$mcp %in% names(clust_infil$analysis_tbl$mcp[[i]])], 
                        split_factor = 'clust_id', 
                        distr_geom = 'box', 
                        cust_theme = globals$common_theme, 
                        plot_title = globals$study_labels[i], 
                        plot_subtitle = clust_infil$n_caps$mcp[[i]], 
                        point_size = 1, 
                        x_lab = 'infiltration estimate Z-score, MCP counter') + 
      scale_fill_manual(values = set_names(globals$cluster_colors, 
                                           c('hi', 'int', 'low'))) + 
      scale_y_discrete(labels = set_names(clust_infil$test$mcp[[i]]$plot_lab, 
                                          clust_infil$test$mcp[[i]]$variable)) + 
      theme(axis.text.y = element_markdown())
    
  }
  
# Box panels with the infiltration estimates, xCell -------
  
  insert_msg('Box plot panels, xCell')
  
  ## I'm displaying Z-scores of the infiltration estimates
  
  for(i in names(clust_infil$analysis_tbl$xcell)) {
    
    clust_infil$box_panels$xcell[[i]] <- 
      draw_violin_panel(data = clust_infil$analysis_tbl$xcell[[i]] %>% 
                          map_dfc(function(x) if(is.numeric(x)) scale(x)[, 1] else x), 
                        variables = clust_infil$variables$xcell[clust_infil$variables$xcell %in% names(clust_infil$analysis_tbl$xcell[[i]])], 
                        split_factor = 'clust_id', 
                        distr_geom = 'box', 
                        cust_theme = globals$common_theme, 
                        plot_title = globals$study_labels[i], 
                        plot_subtitle = clust_infil$n_caps$xcell[[i]], 
                        point_size = 1, 
                        x_lab = 'infiltration estimate Z-score, MCP counter') + 
      scale_fill_manual(values = set_names(globals$cluster_colors, 
                                           c('hi', 'int', 'low'))) + 
      scale_y_discrete(labels = set_names(clust_infil$test$xcell[[i]]$plot_lab, 
                                          clust_infil$test$xcell[[i]]$variable)) + 
      theme(axis.text.y = element_markdown())
    
  }

# Ribbon panels with the infiltration estimates ------
  
  insert_msg('Ribbon panels')
  
  ## with Z-scores
  
  for(i in names(clust_infil$analysis_tbl$mcp)) {
    
    clust_infil$ribbon_panels$mcp[[i]] <- 
      draw_stat_panel(data = clust_infil$analysis_tbl$mcp[[i]] %>% 
                        map_dfc(function(x) if(is.numeric(x)) scale(x)[, 1] else x), 
                      variables = clust_infil$variables$mcp[clust_infil$variables$mcp %in% names(clust_infil$analysis_tbl$mcp[[i]])], 
                      split_factor = 'clust_id', 
                      stat = 'mean', 
                      err_stat = '2se', 
                      form = 'line', 
                      alpha = 0.35,  
                      cust_theme = globals$common_theme, 
                      plot_title = globals$study_labels[i], 
                      plot_subtitle = clust_infil$n_caps$mcp[[i]], 
                      x_lab = 'MCP Counter Z-score, \u00B1 2 \u00D7 SEM') + 
      scale_fill_manual(values = set_names(globals$cluster_colors, 
                                           c('hi', 'int', 'low'))) + 
      scale_color_manual(values = set_names(globals$cluster_colors, 
                                            c('hi', 'int', 'low'))) + 
      scale_y_discrete(labels = set_names(clust_infil$test$mcp[[i]]$plot_lab, 
                                          clust_infil$test$mcp[[i]]$variable)) + 
      theme(axis.text.y = element_markdown(),
            axis.title.y = element_blank())
    
  }
  
  for(i in names(clust_infil$analysis_tbl$xcell)) {
    
    clust_infil$ribbon_panels$xcell[[i]] <- 
      draw_stat_panel(data = clust_infil$analysis_tbl$xcell[[i]] %>% 
                        map_dfc(function(x) if(is.numeric(x)) scale(x)[, 1] else x), 
                      variables = clust_infil$variables$xcell[clust_infil$variables$xcell %in% names(clust_infil$analysis_tbl$xcell[[i]])], 
                      split_factor = 'clust_id', 
                      stat = 'mean', 
                      err_stat = '2se', 
                      form = 'line', 
                      alpha = 0.35,  
                      cust_theme = globals$common_theme, 
                      plot_title = globals$study_labels[i], 
                      plot_subtitle = clust_infil$n_caps$mcp[[i]], 
                      x_lab = 'xCell Z-score, \u00B1 2 \u00D7 SEM') + 
      scale_fill_manual(values = set_names(globals$cluster_colors, 
                                           c('hi', 'int', 'low'))) + 
      scale_color_manual(values = set_names(globals$cluster_colors, 
                                            c('hi', 'int', 'low'))) + 
      scale_y_discrete(labels = set_names(clust_infil$test$mcp[[i]]$plot_lab, 
                                          clust_infil$test$mcp[[i]]$variable)) + 
      theme(axis.text.y = element_markdown(),
            axis.title.y = element_blank())
    
  }
  
# END -----
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()