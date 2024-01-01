# Regulation of genes linked to the common enriched metabolic subsystems 
# between the collagen clusters.

  insert_head()
  
# container -----
  
  ana_metagen <- list()
  
# subsystem - gene mapping ------
  
  insert_msg('Subsystem - gene map')

  ## common regulated subsystemts
  
  ana_metagen$subsystems <- ana_meta$common_significant_subs %>% 
    unlist %>% 
    unique
  
  ana_metagen$subsystems <- 
    c(ana_metagen$subsystems, 'Oxidative phosphorylation')
  
  ## gene - subsystem annotation
  
  ana_metagen$gene_map <- recon$lexicon %>% 
    transmute(subsystem = label, 
              gene_symbol = genes) %>% 
    filter(subsystem %in% ana_metagen$subsystems) %>% 
    blast(subsystem) %>% 
    map(~unname(unlist(.x$gene_symbol))) %>% 
    map(~.x[!is.na(.x)])

# Extraction of regulation estimates -------
  
  insert_msg('Extraction of gene regulation estimates')
  
  ana_metagen$estimates <- ana_metagen$gene_map %>% 
    map(function(sub) ana_dge$test %>% 
          map(filter, gene_symbol %in% sub) %>% 
          map(mutate, 
              plot_lab = ifelse(regulation != 'ns', 
                                gene_symbol, NA)))

# Numbers of genes -------
  
  insert_msg('Numbers of genes')
  
  for(i in names(ana_metagen$estimates)) {
    
    ana_metagen$gene_counts[[i]] <- ana_metagen$estimates[[i]] %>% 
      map(count, 
          regulation, 
          .drop = FALSE) %>% 
      map(~rbind(.x, 
                 tibble(regulation = 'total', 
                        n = sum(.x$n))))
    
    ## ready-to-use captions with total numbers of genes
    ## upregulated and downregulated ones
    
    ana_metagen$gene_caps[[i]] <- ana_metagen$gene_counts[[i]] %>% 
      map(~paste0('total: n = ', .x$n[4], 
                  ', upregulated: n = ', .x$n[1], 
                  ', downregulated: n = ', .x$n[2]))
    
  }
  
# Volcano plots ------
  
  insert_msg('Volcano plots')
  
  for(i in names(ana_metagen$estimates)) {
    
    ana_metagen$volcano_plots[[i]] <- 
      list(data = ana_metagen$estimates[[i]], 
           plot_title = paste(i, 
                              globals$study_labels[names(ana_metagen$estimates[[i]])], 
                              sep = ', '), 
           plot_subtitle = ana_metagen$gene_caps[[i]]) %>% 
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

# END ------
  
  rm(i)
  
  insert_tail()