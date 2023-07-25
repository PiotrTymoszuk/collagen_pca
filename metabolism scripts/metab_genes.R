# Regulation of genes linked to the common enriched metabolic subsystems 
# between the collagen clusters.

  insert_head()
  
# container -----
  
  meta_gene <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals ------
  
  insert_msg('Globals setup')
  
  ## common regulated subsystemts
  
  meta_gene$subsystems <- meta_sub$common %>% 
    map(reduce, union) %>% 
    reduce(union)
  
  ## reactions of the common regulated subsystems
  
  meta_gene$reactions <- meta$regulation[[1]][[1]] %>% 
    filter(subsystem %in% meta_gene$subsystems) %>% 
    mutate(subsystem = factor(subsystem, meta_gene$subsystems)) %>% 
    blast(subsystem) %>% 
    map(~.x$react_id)
  
  ## genes
  
  meta_gene$entrez_id <- meta_gene$reactions %>% 
    map(react_to_gene, meta$models[[1]][[1]]) %>% 
    map(~.x$entrez_id) %>% 
    map(reduce, union)
  
  ## plot title suffixes
  
  meta_gene$title_suffix <- 
    c(int = 'Collagen int vs low', 
      hi = 'Collagen high vs low')
  
# Extraction of regulation estimates -------
  
  insert_msg('Extraction of gene regulation estimates')
  
  meta_gene$estimates <- dge$test_results %>% 
    map(~.x$lm) %>% 
    map(mutate, 
        regulation = ifelse(significant == 'no', 
                            'ns', 
                            ifelse(estimate > log2(1.25), 
                                   'upregulated', 
                                   ifelse(estimate < log2(1.25), 
                                          'downregulated', 'ns'))), 
        regulation = factor(regulation, c('upregulated', 'downregulated', 'ns')))
  
  meta_gene$estimates <- meta_gene$estimates %>% 
    map(function(data) meta_gene$entrez_id %>% 
          map(~filter(data, 
                      entrez_id %in% .x, 
                      level %in% c('Collagen hi', 'Collagen int'))) %>% 
          map(select, -counfounder) %>% 
          map(blast, level) %>% 
          map(set_names, c('hi', 'int'))) %>% 
    map(transpose) %>% 
    transpose %>% 
    map(transpose)
  
# Numbers of genes -------
  
  insert_msg('Numbers of genes')
  
  for(i in names(meta_gene$estimates)) {
    
    meta_gene$gene_counts[[i]] <- meta_gene$estimates[[i]] %>% 
      map(map, 
          count, 
          regulation, 
          .drop = FALSE) %>% 
      map(map, 
          ~rbind(.x, 
                 tibble(regulation = 'total', 
                        n = sum(.x$n))))
    
    ## ready-to-uspe captions with total numbers of genes
    ## upregulated and downregulated ones
    
    meta_gene$gene_caps[[i]] <- meta_gene$gene_counts[[i]] %>% 
      map(map, 
          ~paste0('total: n = ', .x$n[4], 
                  ', upregulated: n = ', .x$n[1], 
                  ', downregulated: n = ', .x$n[2]))
    
  }
  
# Volcano plots ------
  
  insert_msg('Volcano plots')
  
  for(i in names(meta_gene$estimates)) {
    
    meta_gene$volcano_plots[[i]] <- 
      list(data = meta_gene$estimates[[i]], 
           sub = names(meta_gene$estimates[[i]]), 
           caps = meta_gene$gene_caps[[i]]) %>% 
      pmap(function(data, sub, caps) list(data = data, 
                                          plot_title = paste(sub, 
                                                             globals$study_labels[names(data)], 
                                                             sep = ', '), 
                                          plot_subtitle = caps) %>% 
             future_pmap(plot_volcano, 
                         regulation_variable = 'estimate', 
                         p_variable = 'p_adjusted', 
                         signif_level = 0.05, 
                         regulation_level = log2(1.25), 
                         top_significant = 20, 
                         label_variable = 'gene_symbol', 
                         label_type = 'text', 
                         txt_size = 2.5, 
                         txt_face = 'italic', 
                         cust_theme = globals$common_theme, 
                         x_lab = paste0(meta_gene$title_suffix[[i]], 
                                        ', log<sub>2</sub> fold-regulation'), 
                         y_lab = expression('-log'[10] * ' pFDR'), 
                         .options = furrr_options(seed = TRUE)) %>% 
             map(~.x + 
                   theme(plot.tag = element_blank(), 
                         axis.title.x = element_markdown())))
    
    
    
  }
  
# Common up- and downregulated genes --------
  
  insert_msg('Common up and downregulated genes')
  
  ## up- or downregulated in at least 4 out of 5 cohorts
  
  for(i in names(meta_gene$estimates)) {
    
    meta_gene$common[[i]] <- meta_gene$estimates[[i]] %>% 
      map(extract_cmm_subs, m = 4)
    
  }
  
# Heat maps of regulation estimates --------
  
  insert_msg('Heat maps of regulation estimates')

  for(i in names(meta_gene$estimates[[1]])) {
    
    meta_gene$heat_maps[[i]] <- 
      plot_est_heat_map(map2(meta_gene$estimates$hi[[i]], 
                             meta_gene$estimates$int[[i]], 
                             rbind) %>% 
                          map(mutate, 
                              level = factor(level, 
                                             c('Collagen int', 
                                               'Collagen hi'))), 
                        gene_symbols = reduce(meta_gene$common$hi[[i]], c), 
                        plot_title = i, 
                        plot_subtitle = 'common regulated genes',
                        limits = c(-1.5, 1.5), 
                        oob = scales::squish)
    
  }

# END ------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()