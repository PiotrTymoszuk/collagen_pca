# Modulation of signaling pathways with SPIA
# Significantly regulated pathways are identified as pFDR < 0.05 and magnitude
# of regulation expressed as tA > 2

  insert_head()
  
# container -------
  
  dge_spia <- list()

# serial analysis -------
  
  insert_msg('Serial analysis with SPIA')
  
  plan('multisession')
  
  dge_spia$test <- list(de = dge$regulation_vectors, 
                        all = dge$all_vectors) %>% 
    future_pmap(spia, 
                verbose = FALSE, 
                .options = furrr_options(seed = TRUE)) %>% 
    map(as_tibble)
  
  plan('sequential')

# identification of the significantly regulated pathways and common ones -----
  
  insert_msg('Significant pathways')
  
  dge_spia$test <- dge_spia$test %>% 
    map(mutate, 
        significant = ifelse(pGFdr < 0.05 & abs(tA) > 2, 
                             'yes', 'no'), 
        regulation = ifelse(significant == 'yes', 
                            ifelse(tA > 0, 'activated', 'inhibited'), 
                            'ns'), 
        regulation = factor(regulation, c('activated', 'inhibited', 'ns')))
  
  dge_spia$significant <- dge_spia$test %>% 
    map(filter, significant == 'yes')

# Common regulated pathways -------
  
  insert_msg('Commmon regulated pathways')
  
  ## in at least 4 out of five cohorts
  
  dge_spia$cmm_sets <- names(dge_spia$significant) %>% 
    combn(m = 4, simplify = FALSE)
  
  dge_spia$signif <- dge_spia$significant %>% 
    map(blast, Status) %>% 
    transpose %>% 
    map(map, ~.x$Name)
  
  for(i in names(dge_spia$signif)) {
    
    dge_spia$common[[i]] <- dge_spia$cmm_sets %>% 
      map(~dge_spia$signif[[i]][.x]) %>% 
      map(reduce, intersect) %>% 
      reduce(union)

  }
  
  dge_spia$signif <- NULL
  dge_spia <- compact(dge_spia)

# Volcano plots with the pathways ------
  
  insert_msg('Volcano plots')
  
  dge_spia$volcano_plots <- 
    list(data = dge_spia$test, 
         plot_title = globals$study_labels[names(dge_spia$test)]) %>% 
    pmap(plot_volcano, 
         regulation_variable = 'tA', 
         p_variable = 'pGFdr', 
         signif_level = 0.05, 
         regulation_level = 2, 
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
  
# Bubble heat plot with the common regulated pathways ------
  
  insert_msg('Bubble plot with the common regulated pathways')
  
  ## plotting tables
  
  dge_spia$cmm_plot_tbl <- dge_spia$test %>% 
    compress(names_to = 'cohort') %>% 
    filter(Name %in% unique(unlist(dge_spia$common))) %>% 
    mutate(cohort = factor(cohort, names(dge_spia$test)), 
           fontface = ifelse(significant == 'yes', 'bold', 'plain'))

  ## plots
  
  dge_spia$cmm_plot <- dge_spia$cmm_plot_tbl %>% 
    ggplot(aes(x = cohort, 
               y = reorder(Name, tA), 
               fill = regulation, 
               color = regulation, 
               size = tA)) + 
    geom_point(shape = 21, 
               color = 'black') + 
    geom_text(aes(label = signif(tA, 2),
                  alpha = significant, 
                  x = as.numeric(cohort) + 0.25), 
              size = 2.5, 
              hjust = 0, 
              vjust = 0.5, 
              show.legend = FALSE) + 
    scale_color_manual(values = c(activated = 'firebrick', 
                                  inhibited = 'steelblue', 
                                  ns = 'black'), 
                       drop = FALSE, 
                       name = '') + 
    scale_fill_manual(values = c(activated = 'firebrick', 
                                 inhibited = 'steelblue', 
                                 ns = 'gray60'), 
                      drop = FALSE, 
                      name = '') + 
    scale_radius(limits = c(1, 90), 
                 breaks = seq(0, 90, 
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
    labs(title = 'Signaling, collagen hi vs low')

# Caching the results ------
  
  insert_msg('Caching the results')
  
  save(dge_spia, file = './cache/dge_spia.RData')
  
# END -----
  
  rm(i)
  
  insert_tail()