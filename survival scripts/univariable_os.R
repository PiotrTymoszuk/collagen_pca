# Univariable analysis of overall survival in dichotomous gene strata

  insert_head()
  
# container -------
  
  os_cut <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis data frames --------
  
  insert_msg('Analysis data frames')
  
  ## variables
  
  os_cut$variables <- globals$genes_interest$gene_symbol
  
  ## normalized expression
  
  os_cut$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(os_cut$variables)) %>% 
    map(center_data, 'mean') %>% 
    map(rownames_to_column, 'sample_id')
  
  ## overall survival
  
  os_cut$survival <- globals$study_exprs %>% 
    eval %>%
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(safely(select), sample_id, death, os_months) %>% 
    map(~.x$result) %>% 
    compact
  
  ## merging the analysis data frames
  
  os_cut$data <- 
    map2(os_cut$survival, 
         os_cut$expression[names(os_cut$survival)], 
         left_join, by = 'sample_id') %>% 
    map(as_tibble) %>% 
    map(~filter(.x, complete.cases(.x)))
  
# number of samples and minimal number of observations per strata -------
  
  insert_msg('N numbers')
  
  os_cut$n_total <- os_cut$data %>% 
    map_dbl(nrow)
  
  os_cut$n_events <- os_cut$data %>% 
    map_dbl(~sum(.x$death))
  
  ## minimum 25% of the samples in the strata
  
  os_cut$n_min <- floor(0.25 * os_cut$n_total)

# Cutoff finding --------
  
  insert_msg('Cutoff finding')
  
  os_cut$cutoff_obj <- 
    list(x = os_cut$data, 
         y = os_cut$n_min) %>% 
    pmap(function(x, y) os_cut$variables %>% 
           future_map(~find_cutoff(data = x, 
                                   time = 'os_months', 
                                   event = 'death', 
                                   variable = .x, 
                                   min_n = y, 
                                   .parallel = FALSE), 
                      .options = furrr_options(seed = TRUE)) %>% 
           set_names(os_cut$variables))
  
# Testing summary -------
  
  insert_msg('Tetin summary and significant genes')
  
  ## in case there are more cutoff, the first is selected
  
  os_cut$test <- os_cut$cutoff_obj %>% 
    map(map, summary) %>% 
    map(map, ~.x[1, ]) %>% 
    map(compress, names_to = 'gene_symbol') %>% 
    map(re_adjust, 'p_value')
  
  ## appending the test results with the total n numbers
  ## numbers of events, and ready-to-use plot captions and p value labs
  
  os_cut$test <- 
    list(x = os_cut$test, 
         y = names(os_cut$test), 
         v = os_cut$n_total, 
         z = os_cut$n_events) %>% 
    pmap(function(x, y, v, z) x %>% 
           mutate(cohort = globals$study_labels[y], 
                  n_total = v, 
                  n_events = z, 
                  plot_title = html_italic(gene_symbol), 
                  plot_title = paste(plot_title, cohort, sep = ', '), 
                  plot_title = html_bold(plot_title), 
                  plot_subtitle = paste0('total: n = ', n_total, 
                                         ', events: n = ', n_events), 
                  strata_n = paste0('low: n = ', n_low, 
                                   ', high: n = ', n_high), 
                  plot_subtitle = paste0(plot_subtitle, '\n', 
                                         'cutoff = ', signif(cutoff, 2), 
                                         ', ', strata_n)))
  
  
  ## significant effects and significant effect shared by both cohorts
  
  os_cut$significant <- os_cut$test %>% 
    map(filter, p_adjusted < 0.05) %>% 
    map(~.x$gene_symbol)
  
  os_cut$near_significant <- os_cut$test %>% 
    map(filter, p_value < 0.05) %>% 
    map(~.x$gene_symbol)
  
  ## common significant genes: there are no shared genes following
  ## the FDR adjustment!!! I'm working with the 'near-significant' ones
  ## (raw p < 0.05)
  
  os_cut$common_significant <- os_cut$significant %>% 
    reduce(intersect)
  
  os_cut$common_near <- os_cut$near_significant %>% 
    reduce(intersect)
  
# Kaplan-Meier plots for the common near-significant genes -------
  
  insert_msg('Kaplan-Meier plots')
  
  os_cut$km_plots <- os_cut$cutoff_obj %>% 
    map(~.x[os_cut$common_near]) %>% 
    map(map, plot) %>% 
    map(map, ~.x$plot)
  
  ## testing results solely for the significant genes
  
  os_cut$km_test <- os_cut$test %>% 
    map(filter, 
        gene_symbol %in% os_cut$common_near) %>% 
    map(mutate, 
        gene_symbol = factor(gene_symbol, os_cut$common_near)) %>% 
    map(arrange, gene_symbol)
  
  ## plot adjustment 
  
  for(i in names(os_cut$km_plots)) {
    
    os_cut$km_plots[[i]] <- 
      list(x = os_cut$km_plots[[i]], 
           y = os_cut$km_test[[i]]$plot_title, 
           v = os_cut$km_test[[i]]$plot_subtitle, 
           z = os_cut$km_test[[i]]$significance) %>% 
      pmap(function(x, y, v, z) x + 
             labs(title = y, 
                  subtitle = v, 
                  x = 'Overall survival, months',
                  fill = 'Expression strata') + 
             annotate('text', 
                      label = z, 
                      x = 0, 
                      y = 0, 
                      hjust = 0, 
                      vjust = 0, 
                      size = 2.75) + 
             globals$common_theme + 
             theme(plot.tag = element_blank(), 
                   plot.title = element_markdown()))
    
  }
  
  os_cut$km_plots <- transpose(os_cut$km_plots)
  
# Caching the results ------
  
  insert_msg('Cacheing the results')
  
  os_cut <- os_cut[c("cutoff_obj", "test", "significant", 
                     "common_significant", "near_significant", 
                     "common_near", "km_plots")]
  
  save(os_cut, file = './cache/os_cut.RData')
  
# END -------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()