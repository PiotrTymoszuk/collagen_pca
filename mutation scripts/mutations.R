# comparing frequency of the most frequent mutations (>= 2% of participants)
# between the collagen clusters

  insert_head()
  
# container -------
  
  mut_som <- list()
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## variable lexicon
  
  mut_som$lexicon <- 
    tibble(variable = mut_tables$top$mutations, 
           label = mut_tables$top$mutations, 
           plot_title = paste0('<b><em>', 
                               mut_tables$top$mutations, 
                               '</em>, TCGA</b>'))
  
  ## analysis table, exclusion of the outliers
  
  mut_som$analysis_tbl <- mut_tables$mut_tbl %>% 
    select(patient_id, clust_id, all_of(mut_som$lexicon$variable)) %>% 
    filter(!patient_id %in% mut_tables$outliers)
  
  mut_som$analysis_tbl[mut_som$lexicon$variable] <- 
    mut_som$analysis_tbl[mut_som$lexicon$variable] %>% 
    map_dfc(car::recode, "0 = 'WT'; 1 = 'mutated'") %>% 
    map_dfc(factor, c('WT', 'mutated'))
    
# descriptive stats -----
  
  insert_msg('Frequency of mutations in the clusters')
  
  mut_som$stats <- mut_som$analysis_tbl %>% 
    explore(split_factor = 'clust_id', 
            variables = mut_som$lexicon$variable, 
            what = 'table', 
            pub_styled = TRUE) %>% 
    reduce(left_join, by = 'variable') %>% 
    set_names(c('variable', levels(mut_som$analysis_tbl$clust_id)))
  
# serial testing ------
  
  insert_msg('Serial testing')
  
  mut_som$test <- mut_som$analysis_tbl %>% 
    compare_variables(variables = mut_som$lexicon$variable, 
                      split_factor = 'clust_id', 
                      what = 'eff_size', 
                      types = 'cramer_v', 
                      exact = FALSE, 
                      ci = FALSE, 
                      pub_styled = FALSE, 
                      adj_method = 'BH', 
                      .parallel = TRUE, 
                      .paropts = furrr_options(seed = TRUE, 
                                               globals =  c('mut_som'))) %>% 
    mutate(eff_size = paste(estimate_name, 
                            signif(estimate, 2), 
                            sep = ' = '), 
           plot_cap = paste(eff_size, significance, sep = ', '), 
           plot_lab = ifelse(p_adjusted < 0.05, variable, NA))
  
# Plot of the effect sizes -------
  
  insert_msg('Plot of the effect sizes')
  
  mut_som$eff_size_plot <- mut_som$test %>% 
    plot(cust_theme = globals$common_theme, 
         plot_title = 'somatic mutations in the clusters') + 
    geom_hline(yintercept = -log10(0.05), 
               linetype = 'dashed') +
    geom_vline(xintercept = 0.1, 
               linetype = 'dashed') + 
    geom_text_repel(aes(label = plot_lab),
                    size = 2.75, 
                    fontface = 'italic')
  
  mut_som$eff_size_plot <- mut_som$eff_size_plot + 
    labs(x = 'Effect size, Cramer V', 
         y = expression(chi^2 * ' test, -log'[10] * ' pFDR'), 
         subtitle = mut_som$eff_size_plot$labels$tag) + 
    theme(plot.tag = element_blank())
  
# Plots for single mutations ------
  
  insert_msg('Plots for single mutations')
  
  mut_som$plots <- 
    list(variable = mut_som$lexicon$variable, 
         plot_title = mut_som$lexicon$plot_title, 
         plot_subtitle = mut_som$test$plot_cap) %>% 
    pmap(plot_variable, 
         mut_som$analysis_tbl, 
         split_factor = 'clust_id', 
         type = 'stack', 
         scale = 'percent', 
         cust_theme = globals$common_theme, 
         y_lab = '% of cluster', 
         x_n_labs = TRUE) %>% 
    map(~.x + 
          theme(axis.title.x = element_blank(), 
                plot.title = element_markdown()) +
          scale_fill_manual(values = c(WT = 'cornsilk', 
                                       mutated = 'coral3'))) %>% 
    set_names(mut_som$lexicon$variable)

# END ------
  
  insert_tail()