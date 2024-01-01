# Gene set variation analysis for the Recon metabolic subsystem gene signatures
#
# Statistical significance is assessed by two-tailed test and Cohen's d 
# effect size statistic. Differences are considered significant for pFDR < 0.05 
# and at least weak effect size (d >= 0.2). Common significantly 
# regulated cell populations are shared by at least five cohorts except of 
# GSE16560

  insert_head()
  
# container -------
  
  ana_recon <- list()
  
# parallel backend --------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# gene signature scores -------
  
  insert_msg('Signature scores')
  
  ## variables
  
  ana_recon$lexicon <- recon$lexicon
  
  ana_recon$variables <- ana_recon$lexicon$variable
  
  ## data
  
  ana_recon$data <- 
    map2(ana_globals$assignment[names(recon$signatures)], 
         recon$signatures, 
         left_join, by = 'sample_id')
  
  ## n numbers
  
  ana_recon$n_numbers <- ana_recon$data %>% 
    map(extract_n_numbers)
  
# testing --------
  
  insert_msg('Testing')
  
  ana_recon$test <- ana_recon$data %>% 
    future_map(test_two_groups, 
               split_fct = 'clust_id', 
               variables = ana_recon$variables, 
               type = 't', 
               adj_method = 'BH', 
               .options = furrr_options(seed = TRUE)) %>% 
    map(format_two_test)
  
# Significant effects -------
  
  insert_msg('Significant effects')
  
  ana_recon$test <- ana_recon$test %>%
    map(mutate, 
        regulation = ifelse(p_adjusted >= 0.05, 
                            'ns', 
                            ifelse(effect_size >= 0.2, 
                                   'upregulated', 
                                   ifelse(effect_size <= -0.2, 
                                          'downregulated', 'ns'))), 
        regulation = factor(regulation, 
                            c('upregulated', 'downregulated', 'ns')))
  
  ## in single cohorts
  
  ana_recon$significant <- ana_recon$test %>% 
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>% 
    map(blast, regulation) %>% 
    transpose %>% 
    map(map, ~.x$response)
  
  ## shared significant populations
  
  ana_recon$common_significant <- ana_recon$significant %>% 
    map(~.x[names(.x) != 'gse16560']) %>% 
    map(shared_features, m = 5)
  
# Heat maps of the common signatures for single cohorts ---------
  
  insert_msg('Heat maps for single cohorts')
  
  ## signatures significantly regulated in at least five cohorts
  ## are presented
  
  ana_recon$heat_maps <- 
    list(data = ana_recon$data, 
         plot_title = globals$study_labels[names(ana_recon$data)]) %>% 
    pmap(heat_map, 
         variables = reduce(ana_recon$common_significant, union), 
         split_fct = 'clust_id', 
         hide_x_axis_text = TRUE, 
         x_lab = 'Cancer sample', 
         y_lab = 'Recon subsystem gene signature', 
         cust_theme = globals$common_theme) %>% 
    map(~.x + 
          scale_fill_gradient2(low = 'steelblue', 
                               mid = 'black', 
                               high = 'firebrick', 
                               midpoint = 0,
                               limits = c(-3, 3), 
                               oob = scales::squish, 
                               name = 'Z score, ssGSEA') + 
          theme(axis.text.y = element_blank(), 
                axis.ticks.y = element_blank(), 
                axis.line.y = element_blank()))
  
# Heat map of average ssGSEA scores --------
  
  insert_msg('Heat map of average ssGSEA scores')
  
  ana_recon$mean_heat_map <- 
    cohort_heat_map(data_lst = ana_recon$data, 
                    variables = reduce(ana_recon$common_significant, union), 
                    split_fct = 'clust_id', 
                    normalize = FALSE, 
                    plot_title = 'GSVA, Recon subsystem signatures', 
                    plot_subtitle = 'Regulated in at least five cohorts', 
                    y_lab = NULL, 
                    fill_lab = 'mean ssGSEA\nscore') + 
    scale_y_discrete(labels = function(x) exchange(x, ana_recon$lexicon)) + 
    theme(axis.title.x = element_blank()) + 
    guides(x = guide_axis(angle = 90))
  
# END -------
  
  ana_recon <- 
    ana_recon[c("lexicon", "test", "significant", 
                "common_significant", "heat_maps", 
                "mean_heat_map", "top_signatures")]
  
  plan('sequential')
  
  insert_tail()