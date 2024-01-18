# Univariable analysis of overall survival in dichotomous gene strata. 
# significant association with survival is assumed for: 
# pFDR < 0.05 and Cox's HR > 0 (unfavorable marker), and pFDR < 0.05 
# and Cox's HR < 0 (favorable marker).

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

  ## normalized expression of collagen genes
  
  rfs_cut$data <- surv_globals$data
  
  for(i in names(rfs_cut$data)) {
    
    rfs_cut$data[[i]][rfs_cut$variables] <- 
      rfs_cut$data[[i]][rfs_cut$variables] %>% 
      center_data('mean')
    
  }
  
# number of samples and minimal number of observations per strata -------

  rfs_cut$n_min <- floor(0.25 * surv_globals$n_numbers$n_total)

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
  
  insert_msg('Testing summary and significant genes')
  
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
         v = surv_globals$n_numbers$n_total, 
         z = surv_globals$n_numbers$n_events) %>% 
    pmap(function(x, y, v, z) x %>% 
           mutate(cohort = surv_globals$study_labels[y], 
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
  
# Favorable vs unfavorable markers: Cox hazard ratios-----
  
  insert_msg('Cox HR, favorable and unfavorable markers')
  
  rfs_cut$cox_hr <- rfs_cut$cutoff_obj %>% 
    map(map, hr_from_cut) %>% 
    map(compress, names_to = 'gene_symbol')
  
  rfs_cut$test <- map2(rfs_cut$test, 
                       rfs_cut$cox_hr, 
                       left_join, by = 'gene_symbol')
  
  rfs_cut$test <- rfs_cut$test %>% 
    map(mutate, 
        marker = ifelse(p_adjusted < 0.05, 
                        as.character(marker), 'ns'), 
        marker = factor(marker, c('unfavorable', 'favorable', 'ns')), 
        volcano_lab = ifelse(marker != 'ns', 
                             gene_symbol, NA), 
        hr_lab = paste0('HR = ', 
                        signif(hr, 2), 
                        ' [', signif(hr_lower_ci, 2), ' to ', 
                        signif(hr_upper_ci, 2), ']'), 
        km_lab = paste(hr_lab, significance, sep = ', '))
  
# Significant effects --------
  
  insert_msg('Significant effects')
  
  ## significant effects in single cohorts
  
  rfs_cut$significant <- rfs_cut$test %>% 
    map(filter, 
        p_adjusted < 0.05, 
        marker %in% c('unfavorable', 'favorable')) %>%
    map(blast, marker) %>% 
    transpose %>% 
    map(map, ~.x$gene_symbol)

  ## common significant genes: shared by all cohorts
  
  rfs_cut$common_significant <- rfs_cut$significant %>% 
    map(reduce, intersect)

# Kaplan-Meier plots for the common significant genes -------
  
  insert_msg('Kaplan-Meier plots')
  
  rfs_cut$km_plots <- rfs_cut$cutoff_obj %>% 
    map(~.x[reduce(rfs_cut$common_significant, union)]) %>% 
    map(map, plot) %>% 
    map(map, ~.x$plot)
  
  ## testing results solely for the significant genes
  
  rfs_cut$km_test <- rfs_cut$test %>% 
    map(filter, 
        gene_symbol %in% reduce(rfs_cut$common_significant, union)) %>% 
    map(mutate, 
        gene_symbol = factor(gene_symbol, 
                             reduce(rfs_cut$common_significant, union))) %>% 
    map(arrange, gene_symbol) 
  
  ## plot adjustment 
  
  for(i in names(rfs_cut$km_plots)) {
    
    rfs_cut$km_plots[[i]] <- 
      list(x = rfs_cut$km_plots[[i]], 
           y = rfs_cut$km_test[[i]]$plot_title, 
           v = rfs_cut$km_test[[i]]$plot_subtitle, 
           z = rfs_cut$km_test[[i]]$km_lab) %>% 
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
  
# Volcano plots with the HR and significance --------
  
  insert_msg('Volcano plots with the HR and significance')
  
  rfs_cut$volcano_plots <- 
    list(x = rfs_cut$test, 
         y = surv_globals$study_labels[names(rfs_cut$test)]) %>% 
    pmap(function(x, y) x %>% 
           ggplot(aes(x = beta, 
                      y = -log10(p_adjusted), 
                      fill = marker, 
                      color = marker)) + 
           geom_vline(xintercept = 0, 
                      linetype = 'dashed') + 
           geom_hline(yintercept = -log10(0.05),
                      linetype = 'dashed') + 
           geom_point(shape = 21, 
                      size = 2, 
                      color = 'black') + 
           geom_text_repel(aes(label = volcano_lab), 
                           size = 2.5, 
                           fontface = 'italic', 
                           show.legend = FALSE) + 
           scale_fill_manual(values = c(unfavorable = 'firebrick', 
                                        favorable = 'steelblue', 
                                        ns = 'gray70'), 
                             name = 'survival marker') + 
           scale_color_manual(values = c(unfavorable = 'firebrick', 
                                         favorable = 'steelblue', 
                                         ns = 'gray70'), 
                              name = 'survival marker') + 
           globals$common_theme + 
           labs(title = y, 
                x = 'log HR, high vs low expressors', 
                y = expression('-log'[10] * ' pFDR')))
  
# Forest plots: for the common significant genes --------
  
  insert_msg('Forest plots for the common significant genes')
  
  rfs_cut$common_forest_plot <- rfs_cut$km_test %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = factor(cohort, rev(names(rfs_cut$test)))) %>% 
    ggplot(aes(x = beta, 
               y = cohort, 
               color = marker)) + 
    geom_vline(xintercept = 0, 
               linetype = 'dashed') + 
    geom_errorbarh(aes(xmin = beta_lower_ci, 
                       xmax = beta_upper_ci), 
                   height = 0, 
                   position = position_dodge(0.9)) +
    geom_point(shape = 16, 
               size = 2, 
               position = position_dodge(0.9)) + 
    geom_text(aes(label = hr_lab), 
              size = 2.3, 
              hjust = 0.5, 
              vjust = -0.8, 
              position = position_dodge(0.9)) + 
    facet_grid(gene_symbol ~ .) + 
    scale_color_manual(values = c(unfavorable = 'firebrick', 
                                  favorable = 'steelblue', 
                                  ns = 'gray70'), 
                       name = 'survival marker') + 
    scale_y_discrete(labels = surv_globals$study_labels) + 
    globals$common_theme + 
    theme(strip.text.y = element_text(angle = 0, 
                                    hjust = 0, 
                                    face = 'italic')) + 
    labs(title = 'Transcriptional markers of survival', 
         subtitle = 'Biochemical relapse-free survival, shared significant effects', 
         x = 'log HR \u00B1 95% CI, high vs low expressors', 
         y = 'Cohort')
    
    
  
# Caching the results ------
  
  insert_msg('Caching the results')
  
  rfs_cut <- rfs_cut[c("cutoff_obj", "test", "significant", 
                       "common_significant", "km_plots", 
                       "volcano_plots", "common_forest_plot")]
  
  rfs_cut <- compact(rfs_cut)
  
  save(rfs_cut, file = './cache/rfs_cut.RData')
  
# END -------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()