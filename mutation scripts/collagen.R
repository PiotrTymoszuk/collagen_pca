# Genetic alterations of the collagen pathway genes

  insert_head()
  
# container -------
  
  mut_col <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## variables and variable lexicon
  
  mut_col$lexicon <- globals$genes_interest %>% 
    mutate(variable = gene_symbol, 
           plot_title = paste0('<b><em>', variable, '</em>, TCGA</b>'))
  
  ## analysis tables
  
  mut_col$analysis_tbl <- mut_tables[c("mut_tbl", "ampl_tbl", "del_tbl")] %>% 
    map(select, patient_id, clust_id, any_of(mut_col$lexicon$variable))
  
  mut_col$variables <- mut_col$analysis_tbl %>% 
    map(select, any_of(mut_col$lexicon$variable)) %>% 
    map(names)
  
  mut_col$analysis_tbl$mut_tbl[mut_col$variables$mut_tbl] <- 
    mut_col$analysis_tbl$mut_tbl[mut_col$variables$mut_tbl] %>% 
    map_dfc(car::recode, "0 = 'WT'; 1 = 'mutated'") %>% 
    map_dfc(factor, c('WT', 'mutated'))
  
  ## numbers of observations in the clusters
  
  mut_col$n_numbers <- mut_col$analysis_tbl %>% 
    map(count, clust_id)
  
  ## colors 
  
  mut_col$colors <- 
    list(mut_tbl = c('WT' = 'cornsilk', 
                     'mutated' = 'coral3'), 
         ampl_tbl = c('non-amplified' = 'cornsilk', 
                      'amplified' = 'firebrick'), 
         del_tbl = c('non-deleted' = 'cornsilk', 
                     'deleted' = 'steelblue'))
  
# General frequency of alterations ------
  
  insert_msg('General frequency of alterations')
  
  ## frequencies
  
  mut_col$frequency$stats <- mut_col$analysis_tbl %>% 
    map(select, -clust_id, -patient_id) %>% 
    map(map_dfc, function(x) if(is.factor(x)) as.numeric(x) - 1 else x) %>% 
    map(colSums) %>% 
    map(compress, 
        names_to = 'variable', 
        values_to = 'n') %>% 
    map2(., mut_col$analysis_tbl, 
         ~mutate(.x, 
                 n_total = nrow(.y), 
                 percent = n/n_total * 100))
  
  ## plots
  
  mut_col$frequency$plots <- 
    list(x = mut_col$frequency$stats, 
         y = paste0('Collagen pathway, ', 
                    c('gene mutations', 
                      'gene amplifications', 
                      'gene deletions'), 
                    ', TCGA'),
         z = c('coral3', 'firebrick', 'steelblue')) %>% 
    pmap(function(x, y, z) x %>% 
           ggplot(aes(x = percent, 
                      y = reorder(variable, percent))) + 
           geom_bar(stat = 'identity', 
                    color = 'black', 
                    fill = z) + 
           globals$common_theme + 
           theme(axis.title.y = element_blank(), 
                 axis.text.y = element_text(face = 'italic')) + 
           labs(title = y, 
                subtitle = paste('n =', x$n_total[[1]]), 
                x = 'Frequency of alteration, % of cohort'))
  
# Alterations in the clusters ------
  
  insert_msg('Genetic alterations in the clusters')
  
  ## mutation frequencies
  
  mut_col$clusters$stats <- mut_col$analysis_tbl %>% 
    map(select, -patient_id) %>% 
    map(blast, clust_id) %>% 
    map(map, select, any_of(mut_col$lexicon$variable)) %>% 
    map(map, map_dfc, function(x) if(is.factor(x)) as.numeric(x) - 1 else x) %>% 
    map(map, colSums) %>% 
    map(map, 
        compress, 
        names_to = 'variable', 
        values_to = 'n')
  
  mut_col$clusters$stats <- 
    map2(mut_col$clusters$stats, 
         mut_col$n_numbers, 
         function(data, n) map2(data, n$n, 
                                ~mutate(.x, 
                                        n_total = .y, 
                                        percent = n/n_total * 100)))
  
  ## testing for differences between the clusters
  
  mut_col$clusters$test <- 
    future_map2(mut_col$analysis_tbl, 
                mut_col$variables, 
                ~compare_variables(.x, 
                                   variables = .y, 
                                   split_factor = 'clust_id', 
                                   what = 'eff_size', 
                                   types = 'cramer_v', 
                                   exact = FALSE, 
                                   ci = FALSE, 
                                   pub_styled = TRUE, 
                                   adj_method = 'BH'), 
                .options = furrr_options(seed = TRUE)) %>% 
    map(mutate, 
        plot_cap = paste(eff_size, significance, sep = ', '), 
        plot_title = paste0('<b><em>', variable, '</em>, TCGA</b>'))
  
  ## stack plots
  
  for(i in names(mut_col$clusters$test)) {
    
    mut_col$clusters$plots[[i]] <- 
      list(variable = mut_col$clusters$test[[i]]$variable, 
           plot_title = mut_col$clusters$test[[i]]$plot_title,
           plot_subtitle = mut_col$clusters$test[[i]]$plot_cap) %>% 
      future_pmap(plot_variable, 
                  mut_col$analysis_tbl[[i]], 
                  split_factor = 'clust_id', 
                  type = 'stack', 
                  scale = 'percent', 
                  cust_theme = globals$common_theme, 
                  y_lab = '% of cluster', 
                  x_n_labs = TRUE, 
                  .options = furrr_options(seed = TRUE)) %>% 
      map(~.x + 
            scale_fill_manual(values = mut_col$colors[[i]], 
                              name = '') + 
            theme(plot.title = element_markdown(), 
                  axis.title.x = element_blank())) %>% 
      set_names(mut_col$clusters$test[[i]]$variable)
    
  }

# END ------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()