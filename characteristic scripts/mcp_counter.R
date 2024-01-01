# Comparison of non-malignant cell counts in the tumor microenvironment 
# predicted by the MCP Counter algorithm between the collagen clusters. 
# Statistical significance is assessed by Mann-Whitney test and r 
# effect size statistic. Differences are considered significant for pFDR < 0.05 
# and at least weak effect size (r >= 0.1). Common significantly 
# regulated cell populations are shared by at least five cohorts except of 
# GSE16560

  insert_head()
  
# container ------
  
  ana_mcp <- list()
  
# parallel backend --------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# cell infiltration data -------
  
  insert_msg('Cell infiltration data')
  
  ## data
  
  ana_mcp$data <- 
    map2(ana_globals$assignment[names(mcp)], 
         mcp, 
         left_join, by = 'sample_id')
  
  ## variables
  
  ana_mcp$variables <- names(mcp[[1]])[-2:-1]
  
  ## n numbers
  
  ana_mcp$n_numbers <- ana_mcp$data %>% 
    map(extract_n_numbers)
  
# descriptive stats -------
  
  insert_msg('Descriptive stats')
  
  ana_mcp$stats <- ana_mcp$data %>% 
    future_map(explore, 
               variables = ana_mcp$variables, 
               split_factor = 'clust_id', 
               what = 'table', 
               .options = furrr_options(seed = TRUE)) %>% 
    map(format_desc)
  
# Testing -------
  
  insert_msg('Testing')
  
  ana_mcp$test <- ana_mcp$data %>% 
    future_map(test_two_groups, 
               split_fct = 'clust_id', 
               variables = ana_mcp$variables, 
               type = 'wilcox', 
               adj_method = 'BH', 
               .options = furrr_options(seed = TRUE)) %>% 
    map(format_two_test)
  
# Significant effects -------
  
  insert_msg('Significant effects')
  
  ana_mcp$test <- ana_mcp$test %>%
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
  
  ana_mcp$significant <- ana_mcp$test %>% 
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>% 
    map(blast, regulation) %>% 
    transpose %>% 
    map(map, ~.x$response)
  
  ## shared significant populations
  
  ana_mcp$common_significant <- ana_mcp$significant %>% 
    map(~.x[names(.x) != 'gse16560']) %>% 
    map(shared_features, m = 5)
  
# Result table --------
  
  insert_msg('Result table')
  
  ana_mcp$result_tbl <- 
    map2(ana_mcp$stats, 
         map(ana_mcp$test, ~.x[c('variable', 'significance', 'eff_lab')]), 
         left_join, by = 'variable') %>% 
    map2(ana_mcp$n_numbers, ., full_rbind) %>% 
    compress(names_to = 'cohort') %>% 
    format_summ_tbl(rm_n = TRUE) %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    relocate(cohort)
  
# Box plot panels for the shared significant populations -----
  
  insert_msg('Box plots for the significant populations')
  
  ana_mcp$box_plots <- 
    list(variable = unname(unlist(ana_mcp$common_significant)), 
         plot_title = paste0(unname(unlist(ana_mcp$common_significant)), 's')) %>% 
    pmap(plot_box_panel, 
         data_lst = ana_mcp$data, 
         test_lst = ana_mcp$test, 
         normalize = TRUE, 
         norm_center = 'median', 
         point_alpha = 0.5, 
         point_size = 1, 
         x_lab = 'Z scores of cell counts, MCP counter') %>% 
    set_names(unname(unlist(ana_mcp$common_significant)))

# END --------
  
  #ana_mcp$data <- NULL
  
  plan('sequential')
  
  insert_tail()