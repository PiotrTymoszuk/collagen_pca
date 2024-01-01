# Univariable analysis of overall survival in dichotomous gene strata

  insert_head()
  
# container -------
  
  rfs_cut <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis data frames --------
  
  insert_msg('Analysis data frames')
  
  ## variables
  
  rfs_cut$variables <- globals$genes_interest$gene_symbol
  
  ## normalized expression
  
  rfs_cut$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(rfs_cut$variables)) %>% 
    map(center_data, 'mean') %>% 
    map(rownames_to_column, 'sample_id')
  
  ## overall survival
  
  rfs_cut$survival <- globals$study_exprs %>% 
    eval %>%
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(safely(select), sample_id, relapse, rfs_months) %>% 
    map(~.x$result) %>% 
    compact
  
  ## merging the analysis data frames
  
  rfs_cut$data <- 
    map2(rfs_cut$survival, 
         rfs_cut$expression[names(rfs_cut$survival)], 
         left_join, by = 'sample_id') %>% 
    map(as_tibble) %>% 
    map(~filter(.x, complete.cases(.x)))
  
# number of samples and minimal number of observations per strata -------
  
  insert_msg('N numbers')
  
  rfs_cut$n_total <- rfs_cut$data %>% 
    map_dbl(nrow)
  
  rfs_cut$n_events <- rfs_cut$data %>% 
    map_dbl(~sum(.x$relapse))
  
  ## minimum 25% of the samples in the strata
  
  rfs_cut$n_min <- floor(0.25 * rfs_cut$n_total)

# Cutoff finding --------
  
  insert_msg('Cutoff finding')
  
  ## working in a safely mode, since some genes do not yield
  ## any sound cutoffs, e.g. due to constant expression
  
  rfs_cut$cutoff_obj <- 
    list(x = rfs_cut$data, 
         y = rfs_cut$n_min) %>% 
    pmap(function(x, y) rfs_cut$variables %>% 
           future_map(~safely(find_cutoff)(data = x, 
                                           time = 'rfs_months', 
                                           event = 'relapse', 
                                           variable = .x, 
                                           min_n = y, 
                                           .parallel = FALSE), 
                      .options = furrr_options(seed = TRUE)) %>% 
           set_names(rfs_cut$variables))
  
  rfs_cut$cutoff_obj <- rfs_cut$cutoff_obj %>% 
    map(map, ~.x$result) %>% 
    map(compact)
  
# Testing summary -------
  
  insert_msg('Tetin summary and significant genes')
  
  ## in case there are more cutoff, the first is selected
  
  rfs_cut$test <- rfs_cut$cutoff_obj %>% 
    map(map, summary) %>% 
    map(map, ~.x[1, ]) %>% 
    map(compress, names_to = 'gene_symbol') %>% 
    map(re_adjust, 'p_value')
  
  ## appending the test results with the total n numbers
  ## numbers of events, and ready-to-use plot captions and p value labs
  
  rfs_cut$test <- 
    list(x = rfs_cut$test, 
         y = names(rfs_cut$test), 
         v = rfs_cut$n_total, 
         z = rfs_cut$n_events) %>% 
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
  
  rfs_cut$significant <- rfs_cut$test %>% 
    map(filter, p_adjusted < 0.05) %>% 
    map(~.x$gene_symbol)

  ## common significant genes: shared by at least four cohorts
  
  rfs_cut$common_significant <- rfs_cut$significant %>% 
    shared_features(m = 4)
  
# Kaplan-Meier plots for the common near-significant genes -------
  
  insert_msg('Kaplan-Meier plots')
  
  rfs_cut$km_plots <- rfs_cut$cutoff_obj %>% 
    map(~.x[rfs_cut$common_significant]) %>% 
    map(map, plot) %>% 
    map(map, ~.x$plot)
  
  ## testing results solely for the significant genes
  
  rfs_cut$km_test <- rfs_cut$test %>% 
    map(filter, 
        gene_symbol %in% rfs_cut$common_significant) %>% 
    map(mutate, 
        gene_symbol = factor(gene_symbol, rfs_cut$common_significant)) %>% 
    map(arrange, gene_symbol)
  
  ## plot adjustment 
  
  for(i in names(rfs_cut$km_plots)) {
    
    rfs_cut$km_plots[[i]] <- 
      list(x = rfs_cut$km_plots[[i]], 
           y = rfs_cut$km_test[[i]]$plot_title, 
           v = rfs_cut$km_test[[i]]$plot_subtitle, 
           z = rfs_cut$km_test[[i]]$significance) %>% 
      pmap(function(x, y, v, z) x + 
             labs(title = y, 
                  subtitle = v, 
                  x = 'Biochemical relapse-free survival, months',
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
  
  rfs_cut$km_plots <- transpose(rfs_cut$km_plots)
  
# Caching the results ------
  
  insert_msg('Cacheing the results')
  
  rfs_cut <- rfs_cut[c("cutoff_obj", "test", "significant", 
                       "common_significant", "near_significant", 
                       "common_near", "km_plots")]
  
  save(rfs_cut, file = './cache/rfs_cut.RData')
  
# END -------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()