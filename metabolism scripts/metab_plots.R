# Plotting of the metabolic reaction modeling results

  insert_head()
  
# container -----
  
  meta_plots <- list()
  
# globals -----
  
  insert_msg('Globals')
  
  ## common regulated reactions
  
  meta_plots$common_reactions <- 
    meta$common_reactions %>% 
    map(reduce, union)
  
  ## labels
  
  meta_plots$clust_labels <- 
    c(int = 'Collagen int vs Collagen low', 
      hi = 'Collagen high vs Collagen low')
  
# Diagnostic plots: errors ------
  
  insert_msg('Diagnostic plots: errors')
  
  meta_plots$errors <- meta$models %>% 
    map(map, 
        plot, 
        type = 'errors', 
        cust_theme = globals$common_theme)

  ## plot_titles
  
  for(i in names(meta_plots$errors)) {
    
    meta_plots$errors[[i]] <- 
      map2(meta_plots$errors[[i]], 
           paste('Errors and estimates,', globals$study_labels), 
           ~.x + 
             labs(title = .y, 
                  subtitle = meta_plots$clust_labels))
    
  }

# General regulation plots -------
  
  insert_msg('General regulation plots')
  
  meta_plots$regulation <- meta$models %>% 
    map(map, 
        plot, 
        type = 'regulation', 
        show_fit = FALSE, 
        cust_theme = globals$common_theme)
  
  ## plot titles
  
  for(i in names(meta_plots$regulation)) {
    
    meta_plots$regulation[[i]] <- 
      map2(meta_plots$regulation[[i]], 
           paste(meta_plots$clust_labels[[i]], 
                 globals$study_labels, 
                 sep = ', '), 
           ~.x + 
             theme(panel.grid.major.x = element_blank()) + 
             labs(title = .y,, 
                  y = expression('log'[2] * 'fold-regulation')))
    
  }

# Top regulated reactions ------
  
  insert_msg('Top regulated reactions')
  
  meta_plots$top_plots <- meta$models %>% 
    map(map, 
        plot, 
        type = 'top', 
        n_top = 10, 
        cust_theme = globals$common_theme)
  
  ## annotation with titles, subtitles and axis titles
  
  for(i in names(meta_plots$top_plots)) {
    
    meta_plots$top_plots[[i]] <- 
      map2(meta_plots$top_plots[[i]] , 
           paste('Top regulated reactions, ', globals$study_labels), 
           ~.x + 
             labs(title = .y, 
                  subtitle = meta_plots$clust_labels[[i]], 
                  x = expression('log'[2] * 'fold-regulation')))
    
  }

# Common regulated reactions, Forest plots ------
  
  insert_msg('Common regulated reactions, Forest plots')
  
  ## intermediate vs low
  
  for(i in names(meta$models)) {
    
    meta_plots$common_forest[[i]] <- meta$models[[i]] %>% 
      map(plot, 
          type = 'forest', 
          relevant.reaction = meta_plots$common_reactions[[i]], 
          cust_theme = globals$common_theme, 
          show_txt = FALSE, 
          show_ci_txt = FALSE, 
          txt_size = 2.5) %>% 
      map2(., 
           paste(meta_plots$clust_labels, globals$study_labels, sep = ', '), 
           ~.x + 
             labs(title = .y, 
                  subtitle = 'Common significantly regulated reactions', 
                  x = expression('log2'[2] * ' fold-regulation')))
    
    ## faceting by metabolic subsystems
    
    meta_plots$common_forest[[i]] <- 
      meta_plots$common_forest[[i]] %>% 
      map(~.x + 
            facet_grid(subsystem ~ ., 
                       scales = 'free', 
                       space = 'free'))
    
  }

# Reaction numbers ------
  
  insert_msg('Reaction numbers')
  
  ## data
  
  meta_plots$counts$data <- meta$models %>% 
    map(map, count) %>% 
    map(map, 
        filter, 
        subsystem == 'All reactions', 
        status %in% c('activated', 'inhibited')) %>% 
    map(compress, names_to = 'cohort') %>% 
    map(mutate, 
        percent = n/n_total * 100)
  
  ## plots
  
  meta_plots$counts$plots <- meta_plots$counts$data %>% 
    map2(., c('Collagen int vs low', 'Collagen hi vs low'), 
         ~.x %>% 
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
           labs(title = .y, 
                subtitle = paste('total reactions: n =', 
                                 meta_plots$counts$data[[1]]$n_total[[1]]), 
                x = '% of reactions'))
  
# END ------

  rm(i)

  insert_tail()