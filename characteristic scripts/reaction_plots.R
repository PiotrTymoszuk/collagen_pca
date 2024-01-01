# Plotting of the metabolic reaction modeling results

  insert_head()
  
# container -----
  
  meta_plots <- list()

# Diagnostic plots: errors ------
  
  insert_msg('Diagnostic plots: errors')
  
  meta_plots$errors <- meta$models %>% 
    map(plot, 
        type = 'errors', 
        cust_theme = globals$common_theme)

# General regulation plots -------
  
  insert_msg('General regulation plots')
  
  meta_plots$regulation <- meta$models %>% 
    map(plot, 
        type = 'regulation', 
        show_fit = FALSE, 
        cust_theme = globals$common_theme) %>% 
    map2(., 
         globals$study_labels[names(meta$models)], 
         ~.x + 
           theme(panel.grid.major.x = element_blank()) + 
           labs(title = .y,, 
                y = expression('log'[2] * 'fold-regulation')))

# Top regulated reactions ------
  
  insert_msg('Top regulated reactions')
  
  meta_plots$top_plots <- meta$models %>% 
    map(plot, 
        type = 'top', 
        n_top = 10, 
        cust_theme = globals$common_theme) %>% 
    map2(., 
         paste('Top regulated reactions,', 
               globals$study_labels[names(meta$models)]), 
         ~.x + 
           labs(title = .y, 
                x = expression('log'[2] * 'fold-regulation')))

# Reaction numbers ------
  
  insert_msg('Reaction numbers')
  
  ## data
  
  meta_plots$counts$data <- meta$models %>% 
    map(count) %>% 
    map(filter, 
        subsystem == 'All reactions', 
        status %in% c('activated', 'inhibited')) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(percent = n/n_total * 100)
  
  ## plots
  
  meta_plots$counts$plot <- meta_plots$counts$data %>% 
    ggplot(aes(x = percent, 
               y = cohort, 
               fill = status)) + 
    geom_bar(stat = 'identity', 
             color = 'black', 
             position = position_dodge(0.9)) +
    scale_fill_manual(values = c(activated = 'firebrick', 
                                 inhibited = 'steelblue'), 
                      name = '') + 
    scale_y_discrete(labels = globals$study_labels) + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Differentially modulated reactions, hi vs low', 
         subtitle = paste('total reactions: n =', 
                          meta_plots$counts$data$n_total[[1]]), 
         x = '% of reactions')
  
# END ------

  insert_tail()