# Modulation of signaling pathways with SPIA

  insert_head()
  
# container -------
  
  dge_spia <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## plot title prefixes
  
  dge_spia$prefixes <- 
    c(int = 'Collagen int vs low', 
      hi = 'Collagen high vs low')
  
# serial analysis -------
  
  insert_msg('Serial analysis with SPIA')
  
  plan('multisession')
  
  dge_spia$test$int <- list(de = dge$regulation_int, 
                            all = dge$all_vectors) %>% 
    future_pmap(spia, 
                verbose = FALSE, 
                .options = furrr_options(seed = TRUE)) %>% 
    map(as_tibble)
  
  dge_spia$test$hi <- list(de = dge$regulation_hi, 
                            all = dge$all_vectors) %>% 
    future_pmap(spia, 
                verbose = FALSE, 
                .options = furrr_options(seed = TRUE)) %>% 
    map(as_tibble)

  plan('sequential')
  
# identification of the significantly regulated pathways and common ones -----
  
  insert_msg('Significant pathways')
  
  dge_spia$significant <- dge_spia$test %>% 
    map(map, filter, pGFdr < 0.05)

# Common regulated pathways -------
  
  insert_msg('Commmon regulated pathways')
  
  ## in at least 4 out of five cohorts
  
  dge_spia$cmm_sets <- names(dge_spia$significant[[1]]) %>% 
    combn(m = 4, simplify = FALSE) %>% 
    set_names(paste0('set_', 1:length(dge_spia$significant[[1]])))
  
  for(i in names(dge_spia$significant)) {
    
    dge_spia$common[[i]] <- dge_spia$cmm_sets %>% 
      map(function(set) dge_spia$signif[[i]][set] %>% 
            map(blast, Status) %>% 
            map(map, ~.x$Name) %>% 
            transpose %>% 
            map(reduce, intersect)) %>% 
      transpose %>% 
      map(reduce, union)
    
  }

# Volcano plots with the pathways ------
  
  insert_msg('Volcano plots')
  
  for(i in names(dge_spia$test)) {
    
    dge_spia$volcano_plots[[i]] <- 
      list(data = dge_spia$test[[i]], 
           plot_title = paste(dge_spia$prefixes[[i]], 
                               globals$study_labels[names(dge_spia$test[[i]])], 
                              sep = ', ')) %>% 
      pmap(plot_volcano, 
           regulation_variable = 'tA', 
           p_variable = 'pGFdr', 
           signif_level = 0.05, 
           regulation_level = 0, 
           txt_size = 2.5, 
           fill_title = 'Pathway regulation\nvs Collagen low', 
           top_significant = 120, 
           label_type = 'text', 
           label_variable = 'Name', 
           x_lab = paste('Pathway regulation (tA)', 
                         dge_spia$prefixes[[i]], 
                         sep = ', '), 
           y_lab = expression('-log'[10] * 'p FDR'), 
           cust_theme = globals$common_theme) %>% 
      map( ~.x + 
             geom_vline(xintercept = 0, 
                        linetype = 'dashed') + 
             scale_fill_manual(values = c('firebrick', 
                                          'steelblue', 
                                          'gray60'), 
                               labels = c('activated', 
                                          'inhibited', 
                                          'ns'), 
                               name = 'Pathway regulation') + 
             labs(subtitle = .x$labels$tag %>% 
                    stri_replace(fixed = 'upregulated', 
                                 replacement = 'activated') %>% 
                    stri_replace(fixed = 'downregulated', 
                                 replacement = 'inhibited')) + 
             theme(plot.tag = element_blank()))
    
  }

# Bubble heat plot with the common regulated pathways ------
  
  insert_msg('Bubble plot with the common regulated pathways')
  
  ## plotting tables
  
  dge_spia$cmm_plot_tbl <- 
    map2(dge_spia$test, 
         dge_spia$common %>% 
           map(reduce, union), 
         function(data, pathway) data %>% 
           map(filter, Name %in% pathway) %>% 
           compress(names_to = 'cohort')) %>% 
    map(mutate, 
        cohort = factor(cohort, names(dge_spia$test[[1]])), 
        significant = ifelse(pGFdr < 0.05, 'yes', 'no'), 
        fontface = ifelse(significant == 'yes', 'bold', 'plain'))
  
  ## plots
  
  dge_spia$cmm_plots <- 
    map2(dge_spia$cmm_plot_tbl,
         paste('Common regulated pathways', 
               dge_spia$prefixes, 
               sep = ', '), 
         ~ggplot(.x, 
                 aes(x = cohort, 
                     y = reorder(Name, tA), 
                     fill = tA, 
                     size = tA)) + 
           geom_point(shape = 21) + 
           geom_text(aes(label = signif(tA, 2),
                         alpha = significant, 
                         fontface = fontface, 
                         x = as.numeric(cohort) + 0.25), 
                     size = 2.75, 
                     hjust = 0, 
                     vjust = 0.5) + 
           scale_fill_gradient(low = 'white', 
                               high = 'firebrick', 
                               limits = c(2, 175), 
                               name = 'Pathway regulation\ntA', 
                               breaks = seq(0, 175, by = 25)) + 
           scale_radius(limits = c(2, 175), 
                        breaks = seq(0, 175, 
                                     by = 25), 
                        name = 'Pathway regulation\ntA') + 
           scale_x_discrete(labels = globals$study_labels) + 
           scale_alpha_manual(values = c(no = 0.5, 
                                         yes = 1)) + 
           guides(size = 'legend', 
                  fill = 'legend', 
                  alpha = 'none') + 
           globals$common_theme + 
           theme(axis.title = element_blank()) + 
           labs(title = .y))
  
# Caching the results ------
  
  insert_msg('Caching the results')
  
  save(dge_spia, file = './cache/dge_spia.RData')
  
# END -----
  
  rm(i)
  
  insert_tail()