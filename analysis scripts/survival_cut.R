# Survival in dichotomous strata of collagen pathway genes 
# and of the Collagen Score

  insert_head()
  
# container ------
  
  surv_cut <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')

  plan('multisession')
    
# globals ------
  
  insert_msg('Globals')
  
  ## variables and their labels
  
  surv_cut$variables <- c(globals$genes_interest$gene_symbol, 
                          'collagen_score')
  
  surv_cut$var_labs <- ifelse(surv_cut$variables == 'collagen_score', 
                              'Collagen Score', 
                              surv_cut$variables) %>% 
    set_names(surv_cut$variables) %>% 
    compress(names_to = 'variable', 
             values_to = 'label')
  
  surv_cut$var_labs <- surv_cut$var_labs %>% 
    mutate(html_label = ifelse(variable != 'collagen_score', 
                               paste0('<em>', label, '</em>'), 
                               label))
  
  ## analysis tables with the complete survival information
  
  surv_cut$analysis_tbl <- study_data %>% 
    map(~.x$expression) %>% 
    map(extract_tumor_samples) %>% 
    map(select, 
        patient_id, 
        any_of(c('death', 'vitality_fup')), 
        any_of(c('relapse', 'relapse_fup')), 
        all_of(globals$genes_interest$gene_symbol))
  
  surv_cut$analysis_tbl <- 
    map2(surv_cut$analysis_tbl, 
         map(coll_score$score_tbl, ~.x[c('patient_id', 'collagen_score')]), 
         right_join, by = 'patient_id')
  
  surv_cut$analysis_tbl <- 
    map2(surv_cut$analysis_tbl, 
         list(c('death', 'vitality_fup'), 
              c('relapse', 'relapse_fup'), 
              c('relapse', 'relapse_fup'), 
              c('relapse', 'relapse_fup'), 
              c('relapse', 'relapse_fup')), 
         ~filter(.x, 
                 !is.na(.data[[.y[1]]]), 
                 !is.na(.data[[.y[2]]])))

  ## n numbers: total and events
  
  surv_cut$n_numbers <- 
    map2(surv_cut$analysis_tbl, 
         c('death', rep('relapse', 4)), 
         ~count(.x, .data[[.y]])) %>% 
    map(~c(total = sum(.x$n), 
           events = .x$n[2]))

  surv_cut$n_tags <- surv_cut$n_numbers %>% 
    map(~paste0('total: n = ', .x['total'], 
                ', events: n = ', .x['events']))
  
  ## minimal n numbers, 25% of the cohort size
  
  surv_cut$min_n <- surv_cut$analysis_tbl %>% 
    map(~floor(nrow(.x) * 0.25))
  
  ## KM plot axis titles with survival types
  
  surv_cut$km_axis_title <- 
    c('Overall survival, months', 
      rep('Relapse-free survival, months', 4)) %>% 
    set_names(names(surv_cut$analysis_tbl))

# Serial finding of the optimal cutoffs -------
  
  insert_msg('Cutoff finding')
  
  surv_cut$cutoff_obj <- 
    list(data = surv_cut$analysis_tbl, 
         time_var = c('vitality_fup', rep('relapse_fup', 4)), 
         event_var = c('death', rep('relapse', 4)), 
         min_n = surv_cut$min_n) %>% 
    pmap(function(data, time_var, event_var, min_n) surv_cut$variables %>% 
           future_map(~find_cutoff(data = data, 
                                   time = time_var, 
                                   event = event_var, 
                                   variable = .x, 
                                   min_n = min_n, 
                                   .parallel = FALSE), 
                      .options = furrr_options(seed = TRUE)) %>% 
           set_names(surv_cut$variables))
  
# Significance summary -------
  
  insert_msg('Significance summary')
  
  surv_cut$summary_tbl <- surv_cut$cutoff_obj %>% 
    map(map, summary) %>% 
    map(map, function(x) if(is_tibble(x)) x[1, ] else NULL) %>% 
    map(compress, names_to = 'variable')
  
  ## FDR adjustment
  
  surv_cut$summary_tbl <- surv_cut$summary_tbl %>% 
    map(re_adjust) %>% 
    map(mutate, 
        significant = ifelse(p_adjusted < 0.05, 'yes', 'no'), 
        plot_cap = paste0('cutoff = ', signif(cutoff, 2), '\n', significance))
  
  ## numbers of cases in the expression strata
  ## to be shown in the legends of KM plots
  
  surv_cut$km_legends <- surv_cut$summary_tbl %>% 
    map(~map2(.x[['n_low']], 
              .x[['n_high']], 
              ~c(paste('low: n =', .x), 
                 paste('high: n =', .y))))

# Identification of significant markers ------
  
  insert_msg('Identification of significant markers')
  
  surv_cut$significant_fdr <- surv_cut$summary_tbl %>% 
    map(filter, p_adjusted < 0.05) %>% 
    map(~.x$variable)
  
  surv_cut$significant_raw <- surv_cut$summary_tbl %>% 
    map(filter, p_value < 0.05) %>% 
    map(~.x$variable)
  
# Upset plots with the numbers of significant genes -----
  
  insert_msg('Upset plots with the numbers of significant genes')
  
  surv_cut$upset_fdr <- surv_cut$significant_fdr %>% 
    set_names(globals$study_labels[names(surv_cut$significant_fdr)]) %>% 
    map(exchange, surv_cut$var_labs) %>% 
    plot_upset(label_common = FALSE, 
               plot_title = 'Significant survival markers', 
               plot_subtitle = 'Optimal stratification, FDR-corrected', 
               rel_widths = c(0.95, 0.05))
  
  surv_cut$upset_raw <- surv_cut$significant_raw %>% 
    set_names(globals$study_labels[names(surv_cut$significant_raw)]) %>% 
    map(exchange, surv_cut$var_labs) %>% 
    plot_upset(label_common = TRUE, 
               plot_title = 'Significant survival markers', 
               plot_subtitle = 'Optimal stratification, raw p value', 
               fct_per_line = 1)
  
# Plots of the p values -------
  
  insert_msg('Plots of the p values')
  
  surv_cut[c('p_plots_fdr', 'p_plots_raw')] <- 
    list(x = c('p_adjusted', 'p_value'), 
         y = list(expression('-log'[10] * ' pFDR, Mentel-Henszel test'), 
                  expression('-log'[10] * ' p raw, Mentel-Henszel test'))) %>% 
    pmap(function(x, y) list(data = surv_cut$summary_tbl, 
                             plot_title = paste('Survival markers', 
                                                globals$study_labels[names(surv_cut$summary_tbl)], 
                                                sep = ', '), 
                             plot_subtitle = map2(c('OS', rep('RFS', 4)), 
                                                  surv_cut$n_tags, 
                                                  paste, sep = ', ')) %>% 
           future_pmap(plot_signifcant, 
                       p_variable = x, 
                       label_variable = 'variable', 
                       top_significant = 30, 
                       x_lab = y, 
                       cust_theme = globals$common_theme, 
                       .options = furrr_options(seed = TRUE)) %>% 
           map(~.x + 
                 theme(axis.text.y = element_markdown(), 
                       plot.tag = element_blank()) + 
                 scale_y_discrete(labels = function(x) exchange(x, 
                                                                surv_cut$var_labs, 
                                                                value = 'html_label')) + 
                 scale_fill_manual(values = c(ns = 'gray60', 
                                              significant = 'coral3'), 
                                   labels = c(ns = 'ns', 
                                              significant = 'p < 0.05'))))
  
# KM plots ------
  
  insert_msg('KM plots')
  
  ## not done at the moment, memory sparing!
  
 # for(i in names(surv_cut$cutoff_obj)) {
    
    ## base plots
    
  #  surv_cut$km_plots[[i]] <- 
   #   list(x = surv_cut$cutoff_obj[[i]], 
    #       title = names(surv_cut$cutoff_obj[[i]]) %>% 
     #        exchange(surv_cut$var_labs, 
      #                value = 'html_label') %>% 
       #      paste(globals$study_labels[[1]], sep = ', ') %>% 
        #     paste0('<b>', ., '</b>'), 
         #  xlab = surv_cut$km_axis_title[[i]]) %>% 
    #  future_pmap(safely(plot.survcut), 
     #             type = 'km', 
      #            .options = furrr_options(seed = TRUE)) %>% 
    #  map(~.x$result) %>% 
     # compact
    
    ## appending them with the number of cases and events
    ## cutoff and p value is displayed in the plot
    ## numbers of cases in the strata are displayed in the legend
    
  #  surv_cut$km_plots[[i]] <- 
   #   list(x = surv_cut$km_plots[[i]], 
    #       y = names(surv_cut$km_plots[[i]]) %>% 
     #        exchange(surv_cut$var_labs, 
      #                value = 'html_label'), 
       #    z = names(surv_cut$km_plots[[i]]) %>% 
        #     exchange(surv_cut$summary_tbl[[i]], 
         #             key = 'variable', 
          #            value = 'plot_cap'), 
          # v = surv_cut$km_legends[[i]]) %>% 
    #  pmap(function(x, y, z, v) x$plot + 
     #        labs(subtitle = surv_cut$n_tags[[i]], 
      #            color = y) + 
       #      annotate('text', 
        #              label = z, 
         #             x = 0, 
          #            y = 0.1, 
           #           size = 2.75, 
            #          hjust = 0, 
             #         vjust = 0) + 
          #   scale_color_manual(values = c('steelblue', 'firebrick'), 
           #                     labels = v) + 
            # globals$common_theme + 
          #   theme(plot.tag = element_blank(), 
           #        plot.title = element_markdown(), 
            #       legend.title = element_markdown()))
    
  #}

# END -----
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()