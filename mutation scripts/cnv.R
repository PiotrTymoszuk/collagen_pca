# Comparison of CNV frequencies between the clusters
# with Chi-squared test with Cramer's V effect size statistic

  insert_head()
  
# container ------
  
  mut_cnv <- list()
  
# analysis globals --------
  
  insert_msg('Analysis globals')
  
  ## variable lexicons, top CNV
  
  mut_cnv$variables <- mut_tables$top[c("amplifications", "deletions")]
  
  mut_cnv$lexicons <- mut_cnv$variables %>% 
    map(~tibble(variable = .x, 
                label = .x, 
                plot_title = paste0('<em><b>', .x, '</em>, TCGA</b>')))
  
  ## analysis tables, removal of the outliers
  
  mut_cnv$analysis_tbl <- 
    map2(mut_tables[c("ampl_tbl", "del_tbl")], 
         mut_cnv$variables, 
         ~.x[c('patient_id', 'clust_id', .y)]) %>% 
    map(filter, 
        !patient_id %in% mut_tables$outliers)
  
  ## numbers of observations in the clusters
  
  mut_cnv$n_numbers <- mut_cnv$analysis_tbl %>% 
    map(count, clust_id)
  
  ## colors
  
  mut_cnv$colors <- 
    list(ampl_tbl = c('amplified' = 'firebrick', 
                      'non-amplified' = 'cornsilk'), 
         del_tbl = c('deleted' = 'steelblue', 
                     'non-deleted' = 'cornsilk'))
  
# Frequency of CNVs in the clusters ------
  
  insert_msg('Frequency of CNVs in the clusters')
  
  mut_cnv$stats <- 
    list(data = mut_cnv$analysis_tbl, 
         vars = mut_cnv$variables, 
         n = mut_cnv$n_numbers) %>% 
    pmap(function(data, vars, n)  data %>% 
           blast(clust_id) %>% 
           map(~.x[vars]) %>% 
           map(map_dfc, ~as.numeric(.x) - 1) %>% 
           map(colSums) %>% 
           map(compress, 
               names_to = 'variable', 
               values_to = 'n') %>% 
           map2(., n$n, 
                ~mutate(.x, 
                        n_complete = .y, 
                        percent = n/n_complete * 100)))

# Testing for differences -------
  
  insert_msg('Testing')
  
  mut_cnv$test <- 
    map2(mut_cnv$analysis_tbl, 
         mut_cnv$variables, 
         ~compare_variables(.x, 
                            variables = .y, 
                            split_factor = 'clust_id', 
                            what = 'eff_size', 
                            types = 'cramer_v', 
                            exact = FALSE, 
                            ci = FALSE, 
                            pub_styled = FALSE, 
                            adj_method = 'BH', 
                            .parallel = TRUE, 
                            .paropts = furrr_options(seed = TRUE, 
                                                     globals = c('mut_cnv'))))
  
  mut_cnv$test <- mut_cnv$test %>%
    map(mutate, 
        eff_size = paste(estimate_name, 
                         signif(estimate, 2), 
                         sep = ' = '), 
        plot_cap = paste(eff_size, significance, sep = ', '))
  
# Top alterations ------
  
  insert_msg('Top alterations')
  
  ## raw significance and at least weak effect size
  ## just for my curiosity
  
  mut_cnv$top_cnv <- mut_cnv$test %>% 
    map(filter, 
        estimate > 0.1, 
        p_value < 0.05) %>% 
    map(mutate, 
        plot_title = paste0('<b><em>', variable, '</em>, TCGA</b>'))
  
# Plots ------
  
  insert_msg('Plots for single mutations')
  
  ## for the top ones
  
  plan('multisession')
  
  for(i in names(mut_cnv$top_cnv)) {
    
    mut_cnv$top_plots[[i]] <- 
      list(variable = mut_cnv$top_cnv[[i]]$variable, 
           plot_title = mut_cnv$top_cnv[[i]]$plot_title, 
           plot_subtitle = mut_cnv$top_cnv[[i]]$plot_cap) %>% 
      future_pmap(plot_variable, 
                  mut_cnv$analysis_tbl[[i]], 
                  split_factor = 'clust_id', 
                  type = 'stack', 
                  scale = 'percent', 
                  cust_theme = globals$common_theme, 
                  y_lab = '% of cluster', 
                  x_n_labs = TRUE, 
                  .options = furrr_options(seed = TRUE)) %>% 
      map(~.x + 
            scale_fill_manual(values = mut_cnv$colors[[i]]) + 
            theme(axis.title.x = element_blank(), 
                  plot.title = element_markdown())) %>% 
      set_names(mut_cnv$top_cnv[[i]]$variable)
    
    
  }
  
  plan('sequential')
  
# END ------
  
  rm(i)
  
  insert_tail()