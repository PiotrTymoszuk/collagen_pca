# Differences in infiltration between the Collagen clusters

  insert_head()
  
# container ------
  
  dge_infil <- list()
  
# globals: analysis tables -----
  
  insert_msg('Analysis tables')
  
  dge_infil[c('analysis_tbl_xcell', 
              'analysis_tbl_mcp')] <- infil[c('xcell', 
                                              'mcp_counter')] %>% 
    map(~map2(.x, 
              map(dge$analysis_tbl, 
                  ~.x[c('patient_id', 'clust_id')]), 
              inner_join, by = 'patient_id') %>% 
          map(~filter(.x, complete.cases(.x))) %>% 
          map(mutate, 
              clust_id = stri_replace(clust_id, 
                                      fixed = 'Collagen ', 
                                      replacement = ''), 
              clust_id = factor(clust_id, c('low', 'int', 'hi'))))

  ## parallel backend
  
  plan('multisession')

# Descriptive stats ------
  
  insert_msg('Descriptive stats')
  
  dge_infil[c('stats_xcell', 
              'stats_mcp')] <- 
    map2(dge_infil[c('analysis_tbl_xcell', 
                     'analysis_tbl_mcp')],
         infil[c('xcell_types', 
                 'mcp_counter_types')], 
         function(data, vars) data %>% 
           future_map(~explore(data = .x, 
                               split_factor = 'clust_id', 
                               variables = vars[vars %in% names(.x)], 
                               what = 'table', 
                               pub_styled = TRUE), 
                      .options = furrr_options(seed = TRUE)) %>% 
           map(reduce, left_join, by = 'variable') %>% 
           map(set_names, c('variable', 
                            levels(dge_infil$analysis_tbl_xcell[[1]]$clust_id))))
  
# Testing: Kruskal-Wallis test -----
  
  insert_msg('Testing for differences between the clusters')
  
  dge_infil[c('test_xcell', 
              'test_mcp')] <- 
    future_map2(dge_infil[c('analysis_tbl_xcell', 
                            'analysis_tbl_mcp')],
                infil[c('xcell_types', 
                        'mcp_counter_types')], 
                function(data, vars) data %>% 
                  map(~compare_variables(.x, 
                                         split_factor = 'clust_id', 
                                         variables = vars[vars %in% names(.x)], 
                                         what = 'eff_size',
                                         types = 'kruskal_eta', 
                                         ci = FALSE, 
                                         exact = FALSE, 
                                         pub_styled = TRUE, 
                                         adj_method = 'BH')) %>% 
                  map(mutate, 
                      plot_cap = paste(eff_size, significance, sep = ', ')), 
                .options = furrr_options(seed = TRUE))
  
# summary plots ------
  
  insert_msg('Summary effect size and significance plots')
  
  ## etest objects
  
  dge_infil[c('etest_xcell', 
              'etest_mcp')] <- 
    future_map2(dge_infil[c('analysis_tbl_xcell', 
                            'analysis_tbl_mcp')],
                infil[c('xcell_types', 
                        'mcp_counter_types')], 
                function(data, vars) data %>% 
                  map(~compare_variables(.x, 
                                         split_factor = 'clust_id', 
                                         variables = vars[vars %in% names(.x)], 
                                         what = 'eff_size',
                                         types = 'kruskal_eta', 
                                         ci = FALSE, 
                                         exact = FALSE, 
                                         pub_styled = FALSE, 
                                         adj_method = 'BH')), 
                .options = furrr_options(seed = TRUE))
  
  ## plotting
  
  dge_infil[c('sum_plots_xcell', 
              'sum_plots_mcp')] <- c('etest_xcell', 
                                     'etest_mcp') %>% 
    map(function(algo) list(x = dge_infil[[algo]], 
                            plot_title = paste('Infiltration in Collagen clusters,', 
                                               globals$study_labels[names(dge_infil[[algo]])])) %>% 
          pmap(plot, 
               cust_theme = globals$common_theme, 
               plot_subtitle = ' ', 
               show_labels = 'signif') %>% 
          map(~.x + 
                geom_hline(yintercept = -log10(0.05), 
                           linetype = 'dashed') + 
                labs(x = 'Effect size, \u03B7\u00B2', 
                     y = expression('-log'[10] * 'pFDR'))))

# Single plots -------
  
  insert_msg('Violin plots')
  
  ## xCell
  
  dge_infil$violin_plots_xcell <- 
    list(data = dge_infil$analysis_tbl_xcell, 
         cohort_lab = globals$study_labels[names(dge_infil$analysis_tbl_xcell)], 
         stats = dge_infil$test_xcell) %>% 
    pmap(function(data, cohort_lab, stats) list(variable = stats$variable, 
                                                plot_title = paste(stats$variable, 
                                                                   cohort_lab, 
                                                                   sep = ', '), 
                                                plot_subtitle = stats$plot_cap) %>% 
           pmap(plot_variable, 
                data, 
                split_factor = 'clust_id', 
                type = 'violin', 
                point_hjitter = 0, 
                y_lab = 'xCell infiltation estimate', 
                cust_theme = globals$common_theme) %>% 
           map(~.x + 
                 scale_fill_manual(values = c('hi' = 'coral3', 
                                              'int' = 'gray60', 
                                              'low' = 'steelblue3')) + 
                 labs(tag = stri_replace_all(.x$labels$tag, 
                                             fixed = '\n', 
                                             replacement = ', ') %>% 
                        paste0('\n', .))) %>% 
           set_names(stats$variable))
  
  ## MCP Counter
  
  dge_infil$violin_plots_mcp <- 
    list(data = dge_infil$analysis_tbl_mcp, 
         cohort_lab = globals$study_labels[names(dge_infil$analysis_tbl_mcp)], 
         stats = dge_infil$test_mcp) %>% 
    pmap(function(data, cohort_lab, stats) list(variable = stats$variable, 
                                                plot_title = paste(stats$variable, 
                                                                   cohort_lab, 
                                                                   sep = ', '), 
                                                plot_subtitle = stats$plot_cap) %>% 
           pmap(plot_variable, 
                data, 
                split_factor = 'clust_id', 
                type = 'violin', 
                point_hjitter = 0, 
                y_lab = 'MCP Counter infiltration estimate', 
                cust_theme = globals$common_theme) %>% 
           map(~.x + 
                 scale_fill_manual(values = c('hi' = 'coral3', 
                                              'int' = 'gray60', 
                                              'low' = 'steelblue3')) + 
                 labs(tag = stri_replace_all(.x$labels$tag, 
                                             fixed = '\n', 
                                             replacement = ', ') %>% 
                        paste0('\n', .))) %>% 
           set_names(stats$variable))
  
# Ribbon plots with the min/max normalized infiltration estimates ------
  
  insert_msg('Ribbon plots with the normalized infiltration estimates')
  
  ## plotting tables

  dge_infil$ribbon_tbl_xcell <- dge_infil$analysis_tbl_xcell %>% 
    map(~map_dfc(.x, function(x) if(is.numeric(x)) (x - min(x))/(max(x) - min(x)) else x))
  
  dge_infil$ribbon_tbl_mcp <- dge_infil$analysis_tbl_mcp %>% 
    map(~map_dfc(.x, function(x) if(is.numeric(x)) (x - min(x))/(max(x) - min(x)) else x))
  
  ## plot labels with the p values
  
  dge_infil$ribbon_labs_xcell <- dge_infil$test_xcell %>% 
    map(mutate, 
        plot_lab = paste(variable, significance, sep = ', '), 
        plot_lab = stri_replace(plot_lab, 
                                fixed = 'Cancer associated fibroblast', 
                                replacement = 'CAF') %>% 
          stri_replace(fixed = 'Endothelial cell', replacement = 'EC') %>% 
          stri_replace(fixed = 'Myeloid dendritic cell', replacement = 'mDC') %>% 
          stri_replace(fixed = 'Macrophage/Monocyte', replacement = 'Macrophage')) %>% 
    map(~set_names(.x$plot_lab, .x$variable))
  
  dge_infil$ribbon_labs_mcp <- dge_infil$test_mcp %>% 
    map(mutate, 
        plot_lab = paste(variable, significance, sep = ', '), 
        plot_lab = stri_replace(plot_lab, 
                                fixed = 'Cancer associated fibroblast', 
                                replacement = 'CAF') %>% 
          stri_replace(fixed = 'Endothelial cell', replacement = 'EC') %>% 
          stri_replace(fixed = 'Myeloid dendritic cell', replacement = 'mDC') %>% 
          stri_replace(fixed = 'Macrophage/Monocyte', replacement = 'Macrophage')) %>% 
    map(~set_names(.x$plot_lab, .x$variable))
  
  ## xCell ribbon plots
  
  dge_infil$ribbon_plots_xcell <- 
    list(data = dge_infil$ribbon_tbl_xcell, 
         plot_title = paste0('Collagen clusters and infiltration, ', 
                             globals$study_labels[names(dge_infil$ribbon_tbl_xcell)]), 
         plot_tag = map(dge_infil$violin_plots_xcell, 
                        ~.x[[1]]$labels$tag)) %>% 
    pmap(draw_stat_panel, 
         variables = infil$xcell_types,
         split_factor = 'clust_id', 
         stat = 'mean', 
         err_stat = '2se', 
         form = 'line', 
         plot_subtitle = 'min/max normalized xCell estimates', 
         x_lab = 'mean \u00B1 2\u00D7SEM', 
         cust_theme = globals$common_theme) %>% 
    map2(., dge_infil$ribbon_labs_xcell, 
         ~.x + 
          theme(axis.title = element_blank()) + 
          scale_fill_manual(values = c('hi' = 'coral3', 
                                       'int' = 'gray60', 
                                       'low' = 'steelblue3'), 
                            labels = c('hi' = 'Collagen hi', 
                                       'int' = 'Collagen int', 
                                       'low' = 'Collagen low'), 
                            name = 'Cluster') + 
          scale_color_manual(values = c('hi' = 'coral3', 
                                        'int' = 'gray60', 
                                        'low' = 'steelblue3'), 
                             labels = c('hi' = 'Collagen hi', 
                                        'int' = 'Collagen int', 
                                        'low' = 'Collagen low'), 
                             name = 'Cluster') + 
          scale_y_discrete(limits = rev(c('Cancer associated fibroblast', 
                                          'Endothelial cell', 
                                          'T cell CD8+', 
                                          'T cell CD4+ (non-regulatory)', 
                                          'T cell regulatory (T regs)', 
                                          'T cell gamma delta', 
                                          'T cell NK', 
                                          'NK cell', 
                                          'Macrophage M1', 
                                          'Macrophage M2', 
                                          'Monocyte', 
                                          'Neutrophil')), 
                           labels = .y))
  
  ## MCP counter ribbon plots
  
  dge_infil$ribbon_plots_mcp <- 
    list(data = dge_infil$ribbon_tbl_mcp, 
         variables = map(dge_infil$ribbon_tbl_mcp, 
                         ~names(.x)[!names(.x) %in% c('patient_id', 'clust_id')]), 
         plot_title = paste0('Collagen clusters and infiltration, ', 
                             globals$study_labels[names(dge_infil$ribbon_tbl_mcp)]), 
         plot_tag = map(dge_infil$violin_plots_mcp, 
                        ~.x[[1]]$labels$tag)) %>% 
    pmap(draw_stat_panel, 
         split_factor = 'clust_id', 
         stat = 'mean', 
         err_stat = '2se', 
         form = 'line', 
         plot_subtitle = 'min/max normalized MCP counter estimates', 
         x_lab = 'mean \u00B1 2\u00D7SEM', 
         cust_theme = globals$common_theme) %>% 
    map2(., dge_infil$ribbon_labs_mcp, 
        ~.x + 
          theme(axis.title.y = element_blank()) + 
          scale_fill_manual(values = c('hi' = 'coral3', 
                                       'int' = 'gray60', 
                                       'low' = 'steelblue3'), 
                            labels = c('hi' = 'Collagen hi', 
                                       'int' = 'Collagen int', 
                                       'low' = 'Collagen low'), 
                            name = 'Cluster') + 
          scale_color_manual(values = c('hi' = 'coral3', 
                                        'int' = 'gray60', 
                                        'low' = 'steelblue3'), 
                             labels = c('hi' = 'Collagen hi', 
                                        'int' = 'Collagen int', 
                                        'low' = 'Collagen low'), 
                             name = 'Cluster') + 
          scale_y_discrete(limits = rev(c('Cancer associated fibroblast', 
                                          'Endothelial cell', 
                                          'T cell', 
                                          'T cell CD8+', 
                                          'NK cell', 
                                          'B cell', 
                                          'Macrophage/Monocyte', 
                                          'Monocyte', 
                                          'Neutrophil', 
                                          'Myeloid dendritic cell')), 
                           labels = .y))
  
# END -----
  
  plan('sequential')
  
  insert_tail()