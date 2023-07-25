# Comparing counts of mutations, CNV, amplifications and deletions
# between the collagen clusters with Kruskal-Wallis test
#
# Participants with >6 sigma alteration counts are likely outliers 
# and removed from the analysis

  insert_head()
  
# container ------
  
  mut_counts <- list()
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## variable lexicon
  
  mut_counts$lexicon <- 
    c('mut_count' = 'Mutation counts', 
      'cnv_count' = 'CNV counts', 
      'ampl_count' = 'Amplification counts', 
      'del_count' = 'Deletion count') %>% 
    compress(names_to = 'variable', 
             values_to = 'label') %>% 
    mutate(axis_lab = 'number of genes',
           table_lab = paste(label, axis_lab, sep = ', '))
  
  ## analysis table
  
  mut_counts$analysis_tbl <- mut_tables$count_tbl %>% 
    filter(!patient_id %in% mut_tables$outliers)
  
# Descriptive stats ------
  
  insert_msg('Descriptive stats')
  
  mut_counts$stats <- mut_counts$analysis_tbl %>% 
    explore(split_factor = 'clust_id', 
            variables = mut_counts$lexicon$variable, 
            what = 'table', 
            pub_styled = TRUE) %>% 
    reduce(left_join, by = 'variable') %>% 
    set_names(c('variable', levels(mut_counts$analysis_tbl$clust_id)))
  
# Serial testing -----
  
  insert_msg('Testing')
  
  mut_counts$test <- mut_counts$analysis_tbl %>% 
    compare_variables(variables = mut_counts$lexicon$variable, 
                      split_factor = 'clust_id', 
                      what = 'eff_size', 
                      types = 'kruskal_etasq', 
                      exact = FALSE, 
                      ci = FALSE, 
                      pub_styled = TRUE, 
                      adj_method = 'BH') %>% 
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

# Plots -----
  
  insert_msg('Plots')
  
  mut_counts$plots <- 
    list(variable = mut_counts$lexicon$variable, 
         plot_title = paste0(mut_counts$lexicon$label, ', TCGA'), 
         plot_subtitle = mut_counts$test$plot_cap) %>% 
    pmap(plot_variable, 
         mut_counts$analysis_tbl,
         split_factor = 'clust_id', 
         type = 'violin', 
         cust_theme = globals$common_theme,
         y_lab = 'number of genes', 
         point_hjitter = 0, 
         x_n_labs = TRUE) %>% 
    map(~.x + 
          scale_y_continuous(trans = 'sqrt') + 
          scale_fill_manual(values = globals$cluster_colors)) %>% 
    set_names(mut_counts$lexicon$variable)
  
# END ------
  
  insert_tail()