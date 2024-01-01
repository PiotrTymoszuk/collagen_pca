# Graphs for single subsystems subsystems. Essentially these are cascade 
# plots for estimates of metabolic reaction activity.

  insert_head()
  
# container --------
  
  ana_metaplots <- list()

# significantly regulated subsystems, reactions and their activities -------
  
  insert_msg('Significantly reguated subsystems')
  
  # we're working with reactions of the common enriched Recon subsystems
  
  ana_metaplots$subsystems <- ana_meta$common_significant_subs %>% 
    reduce(union)
  
  ana_metaplots$regulation <- ana_meta$models %>% 
    map(components, 'regulation') %>% 
    map(filter, subsystem %in% ana_metaplots$subsystems) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = factor(cohort, names(ana_meta$models)), 
           subsystem = factor(subsystem, ana_metaplots$subsystems)) %>% 
    blast(subsystem) %>% 
    map(blast, cohort)
  
  ana_metaplots$regulation <- ana_metaplots$regulation %>% 
    map(map, filter, !is.na(fold_reg)) %>% 
    map(map, 
        mutate, 
        fold_reg = log2(fold_reg), 
        lower_ci = log2(lower_ci), 
        upper_ci = log2(upper_ci))
  
# Numbers of significantly regulate reactions ---------
  
  insert_msg('Numbers of significantly regulated reactions')
  
  ## numbers of regulated activated and inhibited reactions
  
  ana_metaplots$react_numbers <- ana_meta$models %>% 
    map(count) %>% 
    map(filter, 
        status %in% c('activated', 'inhibited')) %>% 
    map(mutate, 
        percent = n/n_total * 100) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = factor(cohort, names(ana_meta$models)))
  
  ## ready to use captions and X axis titles for cascade plots with numbers of 
  ## significantly
  ## activated and inhibited reactions for the subsystems of interest
  
  ana_metaplots$regulation_caps <- ana_metaplots$react_numbers %>% 
    filter(subsystem %in% ana_metaplots$subsystems) %>% 
    mutate(subsystem = factor(subsystem, ana_metaplots$subsystems)) %>% 
    blast(subsystem) %>% 
    map(blast, cohort)
  
  ana_metaplots$x_titles <- ana_metaplots$regulation_caps %>% 
    map(map, 
        ~paste(.x$subsystem[[1]], .x$n_total[[1]], sep = ', n = '))
  
  ana_metaplots$regulation_caps <- ana_metaplots$regulation_caps %>% 
    map(map, 
        ~paste0('activated: n = ', .x$n[1], 
                ', inhibited: n = ', .x$n[2]))
  
# Plots of percentages of all regulated reactions --------
  
  insert_msg('Plots of percentages of regulated reactions')

  ## bar plot of percentages of regulated reactions
  
  ana_metaplots$react_number_plot <- ana_metaplots$react_numbers %>% 
    filter(subsystem == 'All reactions') %>% 
    mutate(percent = ifelse(status == 'inhibited', -percent, percent)) %>% 
    ggplot(aes(x = percent, 
               y = cohort, 
               fill = status)) + 
    geom_bar(stat = 'identity', 
             color = 'black') + 
    geom_text(aes(label = signif(abs(percent), 2), 
                  x = percent * 0.8), 
              size = 2.5, 
              color = 'white') + 
    scale_x_continuous(labels = function(x) abs(x)) +
    scale_y_discrete(labels = globals$study_labels) + 
    scale_fill_manual(values = c(activated = 'firebrick', 
                                 inhibited = 'steelblue'), 
                      name = '') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Metabolic reaction activity, collagen hi vs low', 
         subtitle = paste('total reactions: n =', 
                          ana_metaplots$react_numbers$n_total[[1]]))

# Cascade plots ------
  
  insert_msg('Cascade plots')
  
  for(i in names(ana_metaplots$regulation)) {
    
    ana_metaplots$regulation_plots[[i]] <- 
      list(data = ana_metaplots$regulation[[i]], 
           plot_title = globals$study_labels[names(ana_metaplots$regulation[[i]])], 
           plot_subtitle = ana_metaplots$regulation_caps[[i]], 
           x_lab = ana_metaplots$x_titles[[i]]) %>% 
      pmap(plot_sign, 
           regulation_variable = 'fold_reg', 
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
  
# Subsystem enrichment ---------
  
  insert_msg('Subsystem enrichment')
  
  ## results of subsystem enrichment analysis for the common significant ones
  
  ana_metaplots$sub_enrichment <- ana_meta$subsystems %>% 
    map(filter, 
        status %in% c('activated', 'inhibited')) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = factor(cohort, names(ana_meta$subsystems))) %>% 
    blast(status)
  
  ana_metaplots$sub_enrichment <- 
    map2_dfr(ana_metaplots$sub_enrichment, 
             ana_meta$common_significant_subs, 
             ~filter(.x, subsystem %in% .y)) %>% 
    mutate(est_plot = ifelse(status == 'inhibited', -OR, OR), 
           plot_status = ifelse(p_adjusted < 0.05 & OR >= 1.44, status, 'ns'), 
           plot_status = factor(plot_status, c('activated', 'inhibited', 'ns'))) %>% 
    blast(cohort)
  
# Bubble plots of the  enrichment OR --------
  
  insert_msg('Plots of the enrichment odds ratios')
  
  ana_metaplots$cmm_plot <-
    regulation_bubble(data_lst = ana_metaplots$sub_enrichment, 
                      variables = ana_metaplots$subsystems, 
                      label_variable = 'subsystem', 
                      regulation_variable = 'est_plot', 
                      status_variable = 'plot_status', 
                      abs_estimate = TRUE, 
                      plot_title = 'Metabolic subsystem enrichment', 
                      plot_subtitle = 'Subsystems enriched in at least 5 cohorts', 
                      fill_lab = 'Regulation status', 
                      size_lab = 'OR') + 
    facet_grid(status ~ ., 
               scales = 'free', 
               space = 'free') + 
    theme(strip.background = element_blank(), 
          strip.text = element_blank())
  
# END ------
  
  insert_tail()