# Gene set variation analysis for the collagen clusters with the Reactome
# pathway ssGSEA scores

  insert_head()
  
# container -----
  
  dge_gsva <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## variables
  
  dge_gsva$lexicon <- reactome$lexicon
  
  ## cluster assignment
  
  dge_gsva$assignment <- coll_clust$assignment %>% 
    map(~mutate(.x, 
                clust_id = factor(clust_id, rev(levels(.x$clust_id)))))
  
  ## signature scores
  
  dge_gsva$analysis_tbl <- 
    map2(dge_gsva$assignment, 
         reactome$signatures, 
         inner_join, by = 'patient_id') %>% 
    map(select, patient_id, clust_id, all_of(dge_gsva$lexicon$variable))
  
# Testing -------
  
  insert_msg('Testing')
  
  ## two-tailed T test: the ssGSVA scores have usually
  ## a nice normal distribution
  ## computation of Cohen's d used later for selection of common significant
  ## signatures and plot ordering
  
  dge_gsva$test <- dge_gsva$analysis_tbl %>% 
    future_map(test_two_groups,
               type = 't', 
               split_fct = 'clust_id', 
               variables = dge_gsva$lexicon$variable, 
               adj_method = 'BH', 
               .parallel = FALSE, 
               .options = furrr_options(seed = TRUE)) %>% 
    map2(., dge_gsva$analysis_tbl, 
         ~mutate(.x,  
                 n = nrow(.y), 
                 se = abs(estimate/stat), 
                 sd = se * sqrt(n), 
                 eff_size = abs(estimate/sd), 
                 eff_size_desc = interpret_cohens_d(eff_size), 
                 regulation = ifelse(significant == 'yes' & 
                                       eff_size > 0.3, 
                                     ifelse(estimate > 0, 
                                            'upregulated', 'downregulated'), 
                                     'ns'), 
                 regulation = factor(regulation, 
                                     c('upregulated', 'downregulated', 'ns'))))
         
# Significant and common signatures ------
  
  insert_msg('Significant and common significant signatures')
  
  ## pFDR < 0.05 and at least small effect size
  
  dge_gsva$significant <- dge_gsva$test %>% 
    map(filter, 
        regulation != 'ns') %>% 
    map(mutate, regulation = droplevels(regulation)) %>% 
    map(blast, regulation) %>% 
    map(map, ~.x$response) %>% 
    transpose
  
  ## common signatures: up- or down-regulated in at least four cohorts
  
  dge_gsva$cmm_sets <- names(dge_gsva$analysis_tbl) %>% 
    combn(m = 4, simplify = FALSE)
  
  for(i in names(dge_gsva$significant)) {
    
    dge_gsva$cmm_signatures[[i]] <- dge_gsva$cmm_sets %>% 
      map(~dge_gsva$significant[[i]][.x]) %>% 
      map(reduce, intersect) %>% 
      reduce(union)
    
  }
  
# Heat map plot ------
  
  insert_msg('Heat map of the means')
  
  ## variables and means
  
  dge_gsva$hm_plot$variables <- dge_gsva$cmm_signatures %>% 
    unlist %>% 
    unique
  
  dge_gsva$hm_plot$data <- dge_gsva$analysis_tbl %>% 
    map(select, clust_id, all_of(dge_gsva$hm_plot$variables)) %>% 
    map(blast, clust_id) %>% 
    map(map, select, -clust_id) %>% 
    map(map, map_dbl, mean) %>% 
    map(map, 
        compress, 
        names_to = 'variable', 
        values_to = 'mean_score') %>% 
    map(compress, names_to = 'clust_id') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(clust_id = stri_extract(as.character(clust_id), 
                                   regex = ('hi|int|low')), 
           clust_id = factor(clust_id, c('low', 'int', 'hi')))
  
  ## appending with the signature classification
  
  dge_gsva$hm_plot$sign_classes <- 
    read_tsv('./data/reactome_classification.tsv') %>% 
    mutate(class = stri_replace_all(class, fixed = '/', replacement = '\n'), 
           class = stri_replace_all(class, fixed = ' ', replacement = '\n'), 
           class = factor(class, 
                          c('ECM\ncollagen\ncytoskeleton', 
                            'adhesion\nmotility\nscavenging', 
                            'angiogenesis\ndevelopment\nhemostasis', 
                            'GF\nsignaling',
                            'immune\nsignaling', 
                            'GPCR\nWNT\nNOTCH', 
                            'metabolism\ntransport\ntranslation')))
  
  dge_gsva$hm_plot$data <- 
    left_join(dge_gsva$hm_plot$data, 
              dge_gsva$hm_plot$sign_classes, 
              by = 'variable') %>% 
    mutate(sign_label = exchange(variable, dge_gsva$lexicon))
  
  ## appending the plotting data with effect sizes
  ## and ANOVA p values
  
  dge_gsva$hm_plot$data <- dge_gsva$hm_plot$data %>% 
    left_join(dge_gsva$test %>% 
                map(~.x[c('response', 'p_adjusted', 'eff_size')]) %>% 
                compress(names_to = 'cohort') %>% 
                mutate(variable = response), 
              by = c('variable', 'cohort'))

  ## heat map
  
  dge_gsva$hm_plot$plot <- dge_gsva$hm_plot$data %>% 
    ggplot(aes(x = cohort, 
               y = reorder(sign_label, eff_size), 
               fill = mean_score)) + 
    geom_tile() + 
    facet_grid(class ~ clust_id, 
               scales = 'free', 
               space = 'free') +
    scale_x_discrete(labels = globals$study_labels, 
                     limits = names(dge_gsva$analysis_tbl)) + 
    scale_fill_gradient2(low = 'steelblue', 
                         mid = 'black', 
                         high = 'firebrick', 
                         midpoint = 0, 
                         limits = c(-0.5, 0.5), 
                         oob = scales::squish) + 
    scale_y_discrete(label = function(x) gsva_labeller(x)) + 
    globals$common_theme + 
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(hjust = 1, 
                                     angle = 90), 
          axis.title = element_blank()) + 
    labs(title = 'GSVA, Reactome pathway signatures', 
         subtitle = 'Regulated in at least four cohorts', 
         fill = 'mean ssGSEA')
  
# Top signatures per category -----
  
  insert_msg('Top signatures for category')

  ## signatures differentiating between the clusters from 
  ## the common regulated ones
  
  ## such signatures will be subsequently presented in detailed plots
  
  ## mean effect sizes
  
  dge_gsva$top_signatures$effect_sizes <- dge_gsva$hm_plot$data %>% 
    filter(clust_id == 'hi') %>% 
    select(variable, cohort, eff_size, class) %>% 
    summarise(eff_size = mean(eff_size), .by = c(variable, class))
  
  ## top 5 signatures per class
  
  dge_gsva$top_signatures$top_variables <- 
    dge_gsva$top_signatures$effect_sizes %>% 
    group_by(class) %>% 
    top_n(n = 5, eff_size) %>% 
    ungroup %>% 
    arrange(class)
  
# END ------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()