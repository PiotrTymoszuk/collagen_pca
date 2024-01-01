# Comparison of non-malignant cell fractions in the tumor microenvironment 
# predicted by the xCell algorithm between the collagen clusters. 
# Statistical significance is assessed by Mann-Whitney test and r 
# effect size statistic. Differences are considered significant for pFDR < 0.05 
# and at least weak effect size (r >= 0.1). Common significantly 
# regulated cell populations are shared by at least five cohorts except of 
# GSE16560

  insert_head()
  
# container ------
  
  ana_xcell <- list()
  
# parallel backend --------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# cell infiltration data -------
  
  insert_msg('Cell infiltration data')
  
  ## data
  
  ana_xcell$data <- 
    map2(ana_globals$assignment[names(xcell)], 
         xcell, 
         left_join, by = 'sample_id')
  
  ## variables
  
  ana_xcell$variables <- names(xcell[[1]])[-2:-1]
  
  ## n numbers
  
  ana_xcell$n_numbers <- ana_xcell$data %>% 
    map(extract_n_numbers)
    
# descriptive stats -------
  
  insert_msg('Descriptive stats')
  
  ana_xcell$stats <- ana_xcell$data %>% 
    future_map(explore, 
               variables = ana_xcell$variables, 
               split_factor = 'clust_id', 
               what = 'table', 
               .options = furrr_options(seed = TRUE)) %>% 
    map(format_desc)
  
# Testing -------
  
  insert_msg('Testing')
  
  ana_xcell$test <- ana_xcell$data %>% 
    future_map(test_two_groups, 
               split_fct = 'clust_id', 
               variables = ana_xcell$variables, 
               type = 'wilcox', 
               adj_method = 'BH', 
               .options = furrr_options(seed = TRUE)) %>% 
    map(format_two_test)
  
# Significant effects -------
  
  insert_msg('Significant effects')
  
  ana_xcell$test <- ana_xcell$test %>%
    map(mutate, 
        regulation = ifelse(p_adjusted >= 0.05, 
                            'ns', 
                            ifelse(effect_size >= 0.1, 
                                   'upregulated', 
                                   ifelse(effect_size <= -0.1, 
                                          'downregulated', 'ns'))), 
        regulation = factor(regulation, 
                            c('upregulated', 'downregulated', 'ns')))
  
  ## in single cohorts
  
  ana_xcell$significant <- ana_xcell$test %>% 
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>% 
    map(blast, regulation) %>% 
    transpose %>% 
    map(map, ~.x$response)
  
  ## shared significant populations
  
  ana_xcell$common_significant <- ana_xcell$significant %>% 
    map(~.x[names(.x) != 'gse16560']) %>% 
    map(shared_features, m = 5)
  
# Result table --------
  
  insert_msg('Result table')
  
  ana_xcell$result_tbl <- 
    map2(ana_xcell$stats, 
         map(ana_xcell$test, ~.x[c('variable', 'significance', 'eff_lab')]), 
         left_join, by = 'variable') %>% 
    map2(ana_xcell$n_numbers, ., full_rbind) %>% 
    compress(names_to = 'cohort') %>% 
    format_summ_tbl(rm_n = TRUE) %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    relocate(cohort)
  
# Box plot panels for the shared significant populations -----
  
  insert_msg('Box plots for the significant populations')
  
  ana_xcell$box_plots <- 
    list(variable = unname(unlist(ana_xcell$common_significant)), 
         plot_title = paste0(unname(unlist(ana_xcell$common_significant)), 's')) %>% 
    pmap(plot_box_panel, 
         data_lst = ana_xcell$data, 
         test_lst = ana_xcell$test, 
         point_alpha = 0.5, 
         point_size = 1, 
         x_lab = 'fraction of tumor sample, xCell') %>% 
    set_names(unname(unlist(ana_xcell$common_significant)))

# END --------
  
  ana_xcell$data <- NULL
  
  plan('sequential')
  
  insert_tail()