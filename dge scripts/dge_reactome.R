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
  
  ## one-way ANOVA: the ssGSVA scores have usually a nice normal distribution
  
  dge_gsva$test <- dge_gsva$analysis_tbl %>% 
    future_map(test_anova, 
               split_fct = 'clust_id', 
               variables = dge_gsva$lexicon$variable, 
               adj_method = 'BH', 
               .parallel = FALSE, 
               .options = furrr_options(seed = TRUE))
  
# Effect size computation ------
  
  insert_msg('Eta-squared')
  
  for(i in names(dge_gsva$analysis_tbl)) {
    
    dge_gsva$eff_size[[i]] <- dge_gsva$lexicon$variable %>% 
      future_map_dbl(get_etasq, 
                     split_factor = 'clust_id', 
                     data = dge_gsva$analysis_tbl[[i]], 
                     .options = furrr_options(seed = TRUE))
    
  }
  
  ## formatting the results
  ## significant effect is defined as eta^2 >= 0.06
  ## (moderate effect)
  
  dge_gsva$eff_size <- dge_gsva$eff_size %>% 
    map(set_names, dge_gsva$lexicon$variable) %>% 
    map(compress, 
        names_to = 'response', 
        values_to = 'eta_sq') %>% 
    map(mutate, 
        significant = ifelse(eta_sq >= 0.06, 'yes', 'no'))
  
# Significance in ANOVA ------
  
  insert_msg('Significance in ANOVA')
  
  ## pFDR < 0.05 and large effect size
  
  dge_gsva$anova_signif <- 
    list(anova = dge_gsva$test %>% 
           map(~.x$anova), 
         eff_size = dge_gsva$eff_size) %>% 
    map(map, 
        filter, 
        significant == 'yes') %>% 
    map(map, ~.x$response)

  dge_gsva$anova_signif <- 
    map2(dge_gsva$anova_signif$anova, 
         dge_gsva$anova_signif$eff_size, 
         intersect)
  
# Significant for the clusters as compared with the collagen low subset ------
  
  insert_msg('Comparison between the clusters')
  
  dge_gsva$lm_signif <- dge_gsva$test %>% 
    map(~.x$lm) %>% 
    map(filter, 
        level != '(Intercept)', 
        significant == 'yes') %>% 
    map2(., dge_gsva$anova_signif, 
         ~filter(.x, response %in% .y)) %>% 
    map(mutate, 
        regulation = ifelse(estimate > 0, 'upregulated', 'downregulated'), 
        regulation = factor(regulation, c('upregulated', 'downregulated')), 
        level = stri_extract(level, regex = 'int|hi'), 
        level = factor(level, c('int', 'hi'))) %>% 
    map(blast, level) %>% 
    transpose %>% 
    map(map, blast, regulation) %>% 
    map(transpose)
  
# Common significantly regulated signatures -----
  
  insert_msg('Commmon significantly regulated signatures')
  
  ## up- or downregulated in at least four cohorts
  
  dge_gsva$cmm_sets <- names(dge_gsva$analysis_tbl) %>% 
    combn(m = 4, simplify = FALSE)
  
  for(i in names(dge_gsva$lm_signif)) {
    
    dge_gsva$cmm_signatures[[i]] <- dge_gsva$lm_signif[[i]] %>% 
      map(function(reg) dge_gsva$cmm_sets %>% 
            map(~reg[.x]) %>% 
            map(map, ~.x$response) %>% 
            map(reduce, intersect) %>% 
            reduce(union))
    
  }
  
# Heat map plots ------
  
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
                          c('ECM\ncollagen', 
                            'adhesion\nmotility\nscavenging', 
                            'angiogenesis\ndevelopment\nhemostasis', 
                            'GF\nsignaling',
                            'immune\nsignaling', 
                            'other\nsignaling', 
                            'metabolism\ntransport\ntranslation')))
  
  dge_gsva$hm_plot$data <- 
    left_join(dge_gsva$hm_plot$data, 
              dge_gsva$hm_plot$sign_classes, 
              by = 'variable') %>% 
    mutate(sign_label = exchange(variable, dge_gsva$lexicon))
  
  ## heat map
  
  dge_gsva$hm_plot$plot <- dge_gsva$hm_plot$data %>% 
    ggplot(aes(x = cohort, 
               y = reorder(sign_label, mean_score), 
               fill = mean_score)) + 
    geom_tile() + 
    facet_grid(class ~ clust_id, 
               scales = 'free', 
               space = 'free') +
    scale_x_discrete(labels = globals$study_labels) + 
    scale_fill_gradient2(low = 'steelblue', 
                         mid = 'black', 
                         high = 'firebrick', 
                         midpoint = 0, 
                         limits = c(-0.65, 0.65)) + 
    scale_y_discrete(label = function(x) gsva_labeller(x)) + 
    globals$common_theme + 
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(hjust = 1, 
                                     angle = 90), 
          axis.title = element_blank()) + 
    labs(title = 'GSVA, Reactome pathway signatures', 
         subtitle = 'Regulated in at least four cohorts, \u03B7\u00B2 \u2265 0.06', 
         fill = 'mean ssGSEA')

# END ------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()