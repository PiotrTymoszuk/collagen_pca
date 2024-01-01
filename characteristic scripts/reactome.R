# Comparison of the Reactome pathway gene signatures between the collagen 
# clusters. 
#
# Statistical significance is assessed by two-tailed test and Cohen's d 
# effect size statistic. Differences are considered significant for pFDR < 0.05 
# and at least weak effect size (d >= 0.2). Common significantly 
# regulated cell populations are shared by at least five cohorts except of 
# GSE16560

  insert_head()
  
# container -------
  
  ana_reactome <- list()

# parallel backend --------

  insert_msg('Parallel backend')
  
  plan('multisession')

# gene signature scores -------

  insert_msg('Signature scores')
  
  ## variables
  
  ana_reactome$lexicon <- reactome$lexicon
  
  ana_reactome$variables <- ana_reactome$lexicon$variable
  
  ## data
  
  ana_reactome$data <- 
    map2(ana_globals$assignment[names(reactome$signatures)], 
         reactome$signatures, 
         left_join, by = 'sample_id')

  ## n numbers
  
  ana_reactome$n_numbers <- ana_reactome$data %>% 
    map(extract_n_numbers)
  
# testing --------
  
  insert_msg('Testing')
  
  ana_reactome$test <- ana_reactome$data %>% 
    future_map(test_two_groups, 
               split_fct = 'clust_id', 
               variables = ana_reactome$variables, 
               type = 't', 
               adj_method = 'BH', 
               .options = furrr_options(seed = TRUE)) %>% 
    map(format_two_test)
  
# Significant effects -------
  
  insert_msg('Significant effects')
  
  ana_reactome$test <- ana_reactome$test %>%
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
  
  ana_reactome$significant <- ana_reactome$test %>% 
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>% 
    map(blast, regulation) %>% 
    transpose %>% 
    map(map, ~.x$response)
  
  ## shared significant populations
  
  ana_reactome$common_significant <- ana_reactome$significant %>% 
    map(~.x[names(.x) != 'gse16560']) %>% 
    map(shared_features, m = 5)
  
# Hierarchical clustering of the common regulated signatures --------
  
  insert_msg('Hierarchical clustering of the sigantures')
  
  ## the signatures clusters based on similarity in expression patterns are 
  ## defined in the TCGA cohort and applied to visualization of common 
  ## regulated signatures in all other data sets
  
  ana_reactome$signature_clusters <- 
    clust_signature(data = ana_reactome$data$tcga, 
                    variables = reduce(ana_reactome$common_significant, 
                                       union), 
                    mds_dim = 3, 
                    fun = hcluster, 
                    distance_method = 'cosine', 
                    k = 4)
  
  ## cluster labels and a modified assignment table to be used 
  ## in heat maps
  
  ana_reactome$cluster_labels <-
    c('1' = 'GF\nsignaling', 
      '2' = 'ECM\ncollagen\ncytoskeleton\nGPCR', 
      '3' = 'inflammatory\nsignaling', 
      '4' = 'aminoacid\nenergy\nmetabolism')
  
  ana_reactome$sign_classification <- ana_reactome$signature_clusters %>% 
    extract('assignment') %>% 
    mutate(clust_id = ana_reactome$cluster_labels[as.character(clust_id)]) %>% 
    set_names(c('variable', 'clust_id'))
  
# Heat maps of the common signatures for single cohorts ---------
  
  insert_msg('Heat maps for single cohorts')
  
  ## signatures significantly regulated in at least five cohorts
  ## are presented
  
  ana_reactome$heat_maps <- 
    list(data = ana_reactome$data, 
         plot_title = globals$study_labels[names(ana_reactome$data)]) %>% 
    pmap(heat_map, 
         variables = reduce(ana_reactome$common_significant, union), 
         split_fct = 'clust_id', 
         variable_classification = ana_reactome$sign_classification, 
         hide_x_axis_text = TRUE, 
         x_lab = 'Cancer sample', 
         y_lab = 'Reactome pathway gene signature', 
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
  
  ana_reactome$mean_heat_map <- 
    cohort_heat_map(data_lst = ana_reactome$data, 
                    variables = reduce(ana_reactome$common_significant, union), 
                    split_fct = 'clust_id', 
                    normalize = FALSE, 
                    variable_classification = ana_reactome$sign_classification, 
                    plot_title = 'GSVA, Reactome pathway signatures', 
                    plot_subtitle = 'Regulated in at least five cohorts', 
                    y_lab = 'Reactome pathway gene signature', 
                    fill_lab = 'mean ssGSEA\nscore') + 
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.line.y = element_blank(), 
          axis.title.x = element_blank()) + 
    guides(x = guide_axis(angle = 90))
  
# Top signatures per category -----
  
  insert_msg('Top signatures for category')
  
  ## signatures differentiating between the clusters from 
  ## the common regulated ones
  
  ## such signatures will be subsequently presented in detailed plots
  
  ## mean effect sizes
  
  ana_reactome$top_signatures <- ana_reactome$test %>% 
    map(filter, 
        variable %in% reduce(ana_reactome$common_significant, union)) %>% 
    map(select, variable, effect_size) %>% 
    map(left_join, ana_reactome$sign_classification, by = 'variable') %>% 
    compress(names_to = 'cohort') %>% 
    group_by(clust_id, variable) %>% 
    summarise(effect_size = mean(effect_size)) %>% 
    ungroup
  
  ana_reactome$top_signatures <- ana_reactome$top_signatures %>% 
    group_by(clust_id) %>% 
    top_n(n = 5, abs(effect_size)) %>% 
    ungroup

# END -------
  
  ana_reactome <- 
    ana_reactome[c("lexicon", "test", "significant", 
                   "common_significant", "signature_clusters", 
                   "cluster_labels", "sign_classification", 
                   "heat_maps", "mean_heat_map", 
                   "top_signatures")]
  
  plan('sequential')
  
  insert_tail()