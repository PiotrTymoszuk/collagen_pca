# Analysis of differentially regulated signaling pathways in the collagen 
# clusters with SPIA.
# Significantly regulated pathways are identified as pFDR < 0.05 and magnitude
# of regulation expressed as tA > 2. Common regulated pathways are defined as 
# pathways significantly activated or inhibited in at least 5 cohorts except of 
# GSE16560

  insert_head()
  
# container -------
  
  ana_spia <- list()
  
# parallel backend --------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# Vectors of regulation estimates and universe vectors --------
  
  insert_msg('Vectors of regulation estimates and all genes')
  
  ## all genes 
  
  ana_spia$all <- ana_globals$genes %>% 
    map(names) %>% 
    map(~.x[!duplicated(.x)]) %>% 
    map(~.x[!is.na(.x)])
  
  ## differentially regulated genes in each cohort
  
  ana_spia$dge <- ana_dge$significant %>% 
    transpose %>% 
    map(reduce, union)
  
  ana_spia$estimates <- 
    map2(ana_dge$test, ana_spia$dge, 
         ~filter(.x, gene_symbol %in% .y)) %>%
    map(filter, 
        !is.na(entrez_id), 
        entrez_id != '') %>% 
    map(filter, 
        !duplicated(entrez_id)) %>% 
    map(~set_names(.x$estimate, .x$entrez_id))
  
# Testing ---------
  
  insert_msg('Testing')
  
  ana_spia$test <- list(de = ana_spia$estimates, 
                        all = ana_spia$all) %>% 
    future_pmap(spia, 
                verbose = FALSE, 
                .options = furrr_options(seed = TRUE)) %>% 
    map(as_tibble)
  
# Significantly activated and inhibited pathways --------
  
  insert_msg('Significantly activated and inhibited pathways')
  
  ana_spia$test <- ana_spia$test %>% 
    map(mutate, 
        significant = ifelse(pGFdr < 0.05, 'yes', 'no'), 
        regulation = ifelse(significant == 'no',
                            'ns', 
                            ifelse(tA > 2, 
                                   'activated', 
                                   ifelse(tA < -2, 'inhibited', 'ns'))), 
        regulation = factor(regulation, c('activated', 'inhibited', 'ns')))
  
  ## significantly activated and inhibited pathways
  
  ana_spia$significant <- ana_spia$test %>% 
    map(filter, regulation %in% c('activated', 'inhibited')) %>% 
    map(blast, regulation) %>% 
    transpose %>% 
    map(map, ~.x$Name)
  
  ana_spia$common_significant <- ana_spia$significant %>% 
    map(~.x[names(.x) != 'gse16560']) %>% 
    map(shared_features, m = 5)
  
# volcano plots of regulation estimates and significance -------
  
  insert_msg('Volcano plots')
  
  ana_spia$volcano_plots <- 
    list(data = ana_spia$test, 
         plot_title = globals$study_labels[names(ana_spia$test)]) %>% 
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
         x_lab = paste('Pathway regulation (tA)'), 
         y_lab = expression('-log'[10] * 'pFDR'), 
         cust_theme = globals$common_theme)
  
  ana_spia$volcano_plots <- ana_spia$volcano_plots %>% 
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
  
  ana_spia$cmm_plot <- 
    regulation_bubble(data_lst = ana_spia$test,
                      variables = reduce(ana_spia$common_significant, union), 
                      label_variable = 'Name', 
                      regulation_variable = 'tA', 
                      status_variable = 'regulation', 
                      plot_title = 'Signaling, collagen hi vs low', 
                      plot_subtitle = paste('Pathways differentially modulated', 
                                            'in at least 5 cohorts'), 
                      limits = c(0, 125), 
                      breaks = seq(0, 125, 
                                   by = 25), 
                      max_size = 5)

# Caching the results -------
  
  insert_msg('Caching the results')
  
  ana_spia <- ana_spia[c("test", "significant", 
                         "common_significant", "volcano_plots", 
                         "cmm_plot")]
  
  save(ana_spia, file = './cache/ana_spia.RData')
  
# END -----
  
  plan('sequential')
  
  insert_tail()