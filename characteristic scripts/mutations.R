# Comparison of mutation frequency between the collagen clusters.
#
# We are analysis all somatic mutations present in at least 2.5% of the samples
# and shared by all cohorts with genetic data available. 
#
# Statistical significance is assessed by Chi-square test with Cramer's V 
# effect size statistic. Significant effects are defined by pFDR < 0.05 and 
# at least weak effect size (V >= 0.1).

  insert_head()
  
# container -------
  
  ana_mut <- list()
  
# analysis data -------
  
  insert_msg('Analysis data')
  
  ana_mut$data <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$mutations) %>% 
    compact
  
  ## cluster assignment
  
  ana_mut$data <- 
    map2(ana_globals$assignment[names(ana_mut$data)], 
         ana_mut$data, 
         left_join, by = 'sample_id') %>% 
    map(filter, !is.na(clust_id))
  
# General mutation frequency in the cohorts ---------
  
  insert_msg('General mutation frequency')
  
  ana_mut$general_freq <- ana_mut$data %>% 
    map(select, -sample_id, -clust_id) %>% 
    map(map_dfc, ~as.numeric(.x) - 1) %>% 
    map(colMeans, na.rm = TRUE) %>% 
    map(compress, 
        names_to = 'gene_symbol', 
        values_to = 'fraction') %>% 
    map(mutate, percent = fraction * 100)
  
  ## top mutations: present in at least 5% of samples in both cohorts
  
  ana_mut$top_mutations <- ana_mut$general_freq %>% 
    map(filter, percent >= 2.5) %>% 
    map(~.x$gene_symbol) %>%
    reduce(intersect)
  
# Stats ---------
  
  insert_msg('Stats')
  
  ## frequencies in the clusters
  
  ana_mut$stats <- ana_mut$data %>% 
    map(select, clust_id, all_of(ana_mut$top_mutations)) %>% 
    map(blast, clust_id, .skip = TRUE) %>% 
    map(map, map_dfc, ~as.numeric(.x) - 1) %>% 
    map(map, colMeans, na.rm = TRUE) %>% 
    map(map, 
        compress, 
        names_to = 'gene_symbol', 
        values_to = 'fraction') %>% 
    map(compress, names_to = 'clust_id')
  
# Testing --------
  
  insert_msg('Testing')
  
  ana_mut$test <- ana_mut$data %>% 
    map(compare_variables, 
        variables = ana_mut$top_mutations, 
        split_factor = 'clust_id', 
        what = 'eff_size', 
        types = 'cramer_v', 
        ci = FALSE, 
        pub_styled = TRUE, 
        adj_method = 'BH') %>% 
    map(mutate, 
        plot_cap = paste(eff_size, significance, sep = ', '))
  
# Significant effects --------
  
  insert_msg('Significant effects')
  
  ## there are no significant differences between the clusters
  
  ana_mut$significant <- ana_mut$test %>% 
    map(filter, p_adjusted < 0.05) %>% 
    map(~.x$variable)
  
# Stack plots --------
  
  insert_msg('Stack plots')
  
  ## for all top mutations
  
  for(i in names(ana_mut$data)) {
    
    ana_mut$plots[[i]] <- 
      list(variable = ana_mut$test[[i]]$variable, 
           plot_title = paste(html_italic(ana_mut$test[[i]]$variable), 
                              globals$study_labels[i], 
                              sep = ', ') %>% 
             html_bold, 
           plot_subtitle = ana_mut$test[[i]]$plot_cap) %>% 
      pmap(plot_variable, 
           ana_mut$data[[i]], 
           split_factor = 'clust_id', 
           scale = 'percent', 
           type = 'stack', 
           cust_theme = globals$common_theme, 
           y_lab = '% of cluster', 
           x_n_labs = TRUE) %>% 
      set_names(ana_mut$test[[i]]$variable)
    
    ana_mut$plots[[i]] <- ana_mut$plots[[i]] %>% 
      map(~.x + 
            scale_fill_manual(values = c(WT = 'cornsilk', 
                                         mutated = 'coral4')) + 
            theme(plot.title = element_markdown()))
    
    
  }

# END -------
  
  rm(i)
  
  ana_mut$data <- NULL
  
  ana_mut <- compact(ana_mut)
  
  insert_tail()