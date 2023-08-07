# Regulation of genes linked to the common enriched metabolic subsystems 
# between the collagen clusters.

  insert_head()
  
# container -----
  
  meta_gene <- list()
  
# globals ------
  
  insert_msg('Globals setup')
  
  ## common regulated subsystemts
  
  meta_gene$subsystems <- meta_sub$common %>% 
    unlist %>% 
    unique
  
  meta_gene$subsystems <- 
    c(meta_gene$subsystems, 'Oxidative phosphorylation')
  
  ## reactions of the common regulated subsystems
  
  meta_gene$reactions <- meta$regulation[[1]] %>% 
    filter(subsystem %in% meta_gene$subsystems) %>% 
    mutate(subsystem = factor(subsystem, meta_gene$subsystems)) %>% 
    blast(subsystem) %>% 
    map(~.x$react_id)
  
  ## genes
  
  meta_gene$entrez_id <- meta_gene$reactions %>% 
    map(react_to_gene, meta$models[[1]]) %>% 
    map(~.x$entrez_id) %>% 
    map(reduce, union)

# Extraction of regulation estimates -------
  
  insert_msg('Extraction of gene regulation estimates')
  
  meta_gene$estimates <- meta_gene$entrez_id %>% 
    map(function(sub) dge$test_results %>% 
          map(filter, entrez_id %in% sub) %>% 
          map(mutate, plot_lab = ifelse(regulation != 'ns', 
                                        gene_symbol, NA)))

# Numbers of genes -------
  
  insert_msg('Numbers of genes')
  
  for(i in names(meta_gene$estimates)) {
    
    meta_gene$gene_counts[[i]] <- meta_gene$estimates[[i]] %>% 
      map(count, 
          regulation, 
          .drop = FALSE) %>% 
      map(~rbind(.x, 
                 tibble(regulation = 'total', 
                        n = sum(.x$n))))
    
    ## ready-to-use captions with total numbers of genes
    ## upregulated and downregulated ones
    
    meta_gene$gene_caps[[i]] <- meta_gene$gene_counts[[i]] %>% 
      map(~paste0('total: n = ', .x$n[4], 
                  ', upregulated: n = ', .x$n[1], 
                  ', downregulated: n = ', .x$n[2]))
    
  }
  
# Volcano plots ------
  
  insert_msg('Volcano plots')
  
  for(i in names(meta_gene$estimates)) {
    
    meta_gene$volcano_plots[[i]] <- 
      list(data = meta_gene$estimates[[i]], 
           plot_title = paste(i, 
                              globals$study_labels[names(meta_gene$estimates[[i]])], 
                              sep = ', '), 
           plot_subtitle = meta_gene$gene_caps[[i]]) %>% 
      pmap(function(data, plot_title, plot_subtitle) data %>% 
             ggplot(aes(x = estimate, 
                        y = -log10(p_adjusted), 
                        color = regulation)) + 
             geom_vline(xintercept = 0, 
                        linetype = 'dashed') + 
             geom_hline(yintercept = -log10(0.05), 
                        linetype = 'dashed') + 
             geom_point(shape = 16, 
                        size = 2) + 
             geom_text_repel(aes(label = plot_lab),
                             size = 2.5, 
                             fontface = 'italic') + 
             scale_color_manual(values = c(upregulated = 'firebrick', 
                                           downregulated = 'steelblue', 
                                           ns = 'gray60'), 
                                name = '') + 
             globals$common_theme + 
             labs(title = plot_title, 
                  subtitle = plot_subtitle, 
                  x = expression('log'[2] * ' fold-regulation'), 
                  y = expression('-log'[10] * ' pFDR')))

  }
  
# Common up- and downregulated genes --------
  
  insert_msg('Common up and downregulated genes')
  
  ## up- or downregulated in at least 4 out of 5 cohorts
  
  meta_gene$common <- meta_gene$estimates %>% 
    map(extract_cmm_subs, m = 4)

# Heat maps of regulation estimates --------
  
  insert_msg('Heat maps of regulation estimates')
  
  meta_gene$heat_maps <- 
    list(estimate_lst = meta_gene$estimates, 
         gene_symbols = map(meta_gene$common, reduce, union), 
         plot_title = names(meta_gene$estimates)) %>% 
    pmap(plot_est_heat_map, 
         plot_subtitle = 'common regulated genes',
         limits = c(-1.5, 1.5), 
         oob = scales::squish)

# END ------
  
  rm(i)
  
  insert_tail()