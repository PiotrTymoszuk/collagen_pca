# Expression of collagen genes and the Collagen Score as a function 
# of demographic and clinical features

  insert_head()
  
# Container ------
  
  coll_clinic <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')

  plan('multisession')
    
# Globals -----
  
  insert_msg('Globals')
  
  ## explanatory variables and responses

  coll_clinic$lexicon <- globals$clinical_lexicon
  
  coll_clinic$responses <- c(globals$genes_interest$gene_symbol, 
                             'collagen_score')
  
  coll_clinic$resp_lexicon <- 
    ifelse(coll_clinic$responses == 'collagen_score', 
                                  'Collagen Score', 
                                  coll_clinic$responses) %>% 
    set_names(coll_clinic$responses) %>% 
    compress(names_to = 'variable', 
             values_to = 'label') %>% 
    mutate(html_label = ifelse(variable == 'collagen_score', 
                               label, 
                               paste0('<em>', label, '</em>')))
  
  ## analysis data, appending with the collagen score values
  
  coll_clinic$analysis_tbl <- study_data %>% 
    map(~.x$expression) %>% 
    map(extract_tumor_samples) %>%
    map(select, 
        patient_id, 
        any_of(coll_clinic$lexicon$variable), 
        any_of(coll_clinic$resp_lexicon$variable)) %>% 
    map2(., 
         map(coll_score$score_tbl, 
             ~.x[c('patient_id', 'collagen_score')]), 
         left_join, by = 'patient_id')
  
  ## formatting of the stage and other clinical variables
  
  coll_clinic$analysis_tbl <- coll_clinic$analysis_tbl %>% 
    map(format_clinical, main_stage = TRUE)

  ## data set-specific variable lists, separate one with factor variables
  
  coll_clinic$var_list <- coll_clinic$analysis_tbl %>%  
    map(names)
  
  for(i in coll_clinic$responses) {
    
    coll_clinic$analysis_tbl <- coll_clinic$analysis_tbl %>% 
      map(mutate, !!i:= as.numeric(.data[[i]]))
    
  }

# Collagen genes and age: Pearson's correlation -------
  
  insert_msg('Collagen genes and age')
  
  ## correlation
  
  coll_clinic$age_corr <- 
    coll_clinic$analysis_tbl[c("GSE16560", "GSE70768", "GSE116918", "tcga")] %>% 
    future_map(function(cohort) map(coll_clinic$responses, 
                                    ~c('age', .x)) %>% 
                 map_dfr(~correlate_variables(cohort, 
                                              variables = .x, 
                                              what = 'correlation', 
                                              type = 'pearson', 
                                              ci = TRUE, 
                                              pub_styled = TRUE)), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(re_adjust)
  
  ## plot titles and captions
  
  coll_clinic$age_corr <- coll_clinic$age_corr %>% 
    map2(., names(.), 
         format_corr_results, 
         dict = coll_clinic$resp_lexicon) %>% 
    map(mutate, x_lab = 'Age, years')

  ## significant correlations
  
  coll_clinic$age_significant <- coll_clinic$age_corr %>% 
    map(filter, p_adjusted < 0.05)
  
  ## correlation plots
  
  coll_clinic$age_plots <- 
    list(data = coll_clinic$analysis_tbl[names(coll_clinic$age_corr)], 
         test_results = coll_clinic$age_corr, 
         point_color = globals$study_colors[names(coll_clinic$age_corr)]) %>% 
    future_pmap(plot_clinical_correlation, 
                .options = furrr_options(seed = TRUE))
  
# Collagen genes and PSA an diagnosis: Spearman correlation -----
  
  insert_msg('Collagen genes and PSA')
  
  ## correlation
  
  coll_clinic$psa_corr <- 
    coll_clinic$analysis_tbl[c("GSE70768", "GSE70769", "GSE116918", "tcga")] %>% 
    future_map(function(cohort) map(coll_clinic$responses, 
                                    ~c('psa_at_diagnosis', .x)) %>% 
                 map_dfr(~correlate_variables(cohort[.x] %>% 
                                                filter(complete.cases(.)), 
                                              variables = .x, 
                                              what = 'correlation', 
                                              type = 'spearman', 
                                              ci = TRUE, 
                                              pub_styled = TRUE)), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(re_adjust)
  
  ## plot title, subtitle and y axis labels
  
  coll_clinic$psa_corr <- coll_clinic$psa_corr %>% 
    map2(., names(.),
         format_corr_results, 
         dict = coll_clinic$resp_lexicon) %>% 
    map(mutate, x_lab = 'PSA at diagnosis')

  ## significant correlations
  
  coll_clinic$psa_significant <- coll_clinic$psa_corr %>% 
    map(filter, p_adjusted < 0.05)
  
  ## correlation plots
  
  coll_clinic$psa_plots <- 
    list(data = coll_clinic$analysis_tbl[names(coll_clinic$psa_corr)], 
         test_results = coll_clinic$psa_corr, 
         point_color = globals$study_colors[names(coll_clinic$psa_corr)]) %>% 
    future_pmap(plot_clinical_correlation, 
                .options = furrr_options(seed = TRUE))

# Collagen genes and numeric Gleason score: Spearman correlation -------  
  
  insert_msg('Collagen genes and Gleason score, correlation')
  
  ## correlation
  
  coll_clinic$gleason_corr <- coll_clinic$analysis_tbl %>% 
    future_map(function(cohort) map(coll_clinic$responses, ~c('gleason', .x)) %>% 
                 map_dfr(~correlate_variables(cohort[.x] %>% 
                                                filter(complete.cases(.)), 
                                              variables = .x, 
                                              what = 'correlation', 
                                              type = 'spearman', 
                                              ci = TRUE, 
                                              pub_styled = TRUE)), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(re_adjust)
  
  ## plot title, subtitle and Y axis
  
  coll_clinic$gleason_corr <- coll_clinic$gleason_corr %>% 
    map2(., names(.), 
         format_corr_results, 
         dict = coll_clinic$resp_lexicon) %>% 
    map(mutate, x_lab = 'Gleason score, sum')
  
  ## significant correlations
  
  coll_clinic$gleason_significant <- coll_clinic$gleason_corr %>% 
    map(filter, p_adjusted < 0.05)
  
  ## correlation plots
  
  coll_clinic$gleason_plots <- 
    list(data = coll_clinic$analysis_tbl, 
         test_results = coll_clinic$gleason_corr, 
         point_color = globals$study_colors[names(coll_clinic$analysis_tbl)]) %>% 
    future_pmap(plot_clinical_correlation, 
                .options = furrr_options(seed = TRUE))

# Pathological tumor stage -------
  
  insert_msg('Collagen genes and pathology stage')
  
  ## descriptive stats
  
  coll_clinic$path_tumor_stats <- 
    coll_clinic$analysis_tbl[c("GSE70768", "GSE70769", "GSE116918", "tcga")] %>% 
    future_map(~explore(filter(.x, !is.na(pathology_stage_tumor)), 
                        split_factor = 'pathology_stage_tumor', 
                        variables = coll_clinic$responses, 
                        what = 'table', 
                        pub_styled = TRUE) %>% 
                 reduce(left_join, by = 'variable'), 
               .options = furrr_options(seed = TRUE))
  
  ## testing
  
  coll_clinic$path_tumor_test <- 
    coll_clinic$analysis_tbl[c("GSE70768", "GSE70769", "GSE116918", "tcga")] %>% 
    future_map(~compare_variables(.x %>% 
                                    filter(!is.na(pathology_stage_tumor)), 
                                  split_factor = 'pathology_stage_tumor', 
                                  variables = coll_clinic$responses, 
                                  what = 'eff_size',
                                  types = 'etasq', 
                                  ci = FALSE, 
                                  exact = FALSE, 
                                  adj_method = 'BH', 
                                  pub_styled = TRUE), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(mutate, plot_cap = paste(eff_size, significance, sep = ', '))
  
  ## significant effects
  
  coll_clinic$path_tumor_significant <- coll_clinic$path_tumor_test %>% 
    map(filter, p_adjusted < 0.05) %>% 
    map(~.x$variable)
  
  ## plotting
  
  coll_clinic$path_tumor_plots <- 
    list(data = coll_clinic$analysis_tbl[c("GSE70768", "GSE70769", "GSE116918", "tcga")], 
         cohort_name = c("GSE70768", "GSE70769", "GSE116918", "tcga"), 
         test_results = coll_clinic$path_tumor_test) %>% 
    future_pmap(plot_clinical_violin, 
                split_factor = 'pathology_stage_tumor', 
                variables = coll_clinic$responses, 
                dict = coll_clinic$resp_lexicon, 
                palette = 'Blues', 
                x_lab = 'pathological stage', 
                .options = furrr_options(seed = TRUE))
  
# Pathological node and metastasis stages ------
  
  insert_msg('Pathological node and metastasis stages')
  
  ## descriptive stats
  
  coll_clinic$stages_stats <- 
    coll_clinic$analysis_tbl[c("GSE70768", "GSE70769", "tcga")] %>% 
    future_map(function(cohort) c(pathology_stage_node = 'pathology_stage_node', 
                                  pathology_stage_meta = 'pathology_stage_meta') %>% 
                 map(~explore(filter(cohort, !is.na(.data[[.x]])), 
                              split_factor = .x, 
                              variables = coll_clinic$responses, 
                              what = 'table', 
                              pub_styled = TRUE) %>% 
                       reduce(left_join, by = 'variable') %>% 
                       set_names(c('variable', levels(factor(cohort[[.x]]))))), 
               .options = furrr_options(seed = TRUE))
  
  ## testing for differences: one-way ANOVA
  
  coll_clinic$stages_test <- 
    coll_clinic$analysis_tbl[c("GSE70768", "GSE70769", "tcga")] %>% 
    future_map(function(cohort)  
      map2(c(pathology_stage_node = 'pathology_stage_node', 
             pathology_stage_meta = 'pathology_stage_meta'), 
           c('cohen_d', 'cohen_d'), 
           ~safely(compare_variables)(filter(cohort, !is.na(.data[[.x]])), 
                                      split_factor = .x, 
                                      variables = coll_clinic$responses, 
                                      what = 'eff_size',
                                      types = .y, 
                                      ci = FALSE, 
                                      pub_styled = TRUE, 
                                      adj_method = 'BH')) %>% 
        map(~.x$result) %>% 
        compact %>% 
        map(mutate, plot_cap = paste(eff_size, significance, sep = ', ')), 
      .options = furrr_options(seed = TRUE))

  ## significant effects
  
  coll_clinic$stages_significant <- coll_clinic$stages_test %>% 
    map(map, filter, p_adjusted < 0.05) %>% 
    map(map, ~.x$variable)

# Plotting of the pathological node stage -----
  
  insert_msg('Plotting of the pathological node stages')
  
  coll_clinic$node_stage_plots <- 
    list(data = coll_clinic$analysis_tbl[c("GSE70768", "tcga")],
         cohort_name = c("GSE70768", "tcga"), 
         test_results = coll_clinic$stages_test[c("GSE70768", "tcga")] %>% 
           map(~.x$pathology_stage_node)) %>% 
    future_pmap(plot_clinical_violin, 
                split_factor = 'pathology_stage_node', 
                variables = coll_clinic$responses, 
                dict = coll_clinic$resp_lexicon, 
                palette = 'Blues', 
                x_lab = 'Node pathology stage',
                .options = furrr_options(seed = TRUE))
  
# Plotting of the pathological metastasis stage -----
  
  insert_msg('Plotting of the metastasis node stages')
  
  coll_clinic$meta_stage_plots <- 
    list(data = coll_clinic$analysis_tbl[c("GSE70769", "tcga")],
         cohort_name = c("GSE70769", "tcga"), 
         test_results = coll_clinic$stages_test[c("GSE70769", "tcga")] %>% 
           map(~.x$pathology_stage_meta)) %>% 
    future_pmap(plot_clinical_violin, 
                split_factor = 'pathology_stage_meta', 
                variables = coll_clinic$responses, 
                dict = coll_clinic$resp_lexicon, 
                palette = 'Blues', 
                x_lab = 'Metastasis pathology stage', 
                .options = furrr_options(seed = TRUE))
  
# Extra-capsular extension -------
  
  insert_msg('Collagen genes and extracapsular extension')
  
  ## descriptive stats

  coll_clinic$ece_stats <- 
    coll_clinic$analysis_tbl[c('GSE70768', 'GSE70769')] %>% 
    map(~explore(filter(.x, !is.na(extra_capsular_extension)), 
                 split_factor = 'extra_capsular_extension', 
                 variables = coll_clinic$responses, 
                 what = 'table', 
                 pub_styled = TRUE)) %>% 
    map(reduce, left_join, by = 'variable') %>% 
    map(set_names, c('variable', 'ece_negative', 'ece_positive'))
  
  ## testing
  
  coll_clinic$ece_test <- 
    coll_clinic$analysis_tbl[c('GSE70768', 'GSE70769')] %>% 
    future_map(~compare_variables(.x %>% 
                                    filter(!is.na(extra_capsular_extension)) %>% 
                                    mutate(extra_capsular_extension = factor(extra_capsular_extension, 
                                                                             c('yes', 'no'))), 
                                  split_factor = 'extra_capsular_extension', 
                                  variables = coll_clinic$responses, 
                                  what = 'eff_size',
                                  types = 'cohen_d', 
                                  ci = FALSE, 
                                  adj_method = 'BH', 
                                  pub_styled = TRUE), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(mutate, plot_cap = paste(eff_size, significance, sep = ', '))
  
  ## significant effects
  
  coll_clinic$ece_significant <- coll_clinic$ece_test %>% 
    map(filter, p_adjusted < 0.05) %>% 
    map(~.x$variable)
  
  ## plotting
  
  coll_clinic$ece_plots <- 
    list(data = coll_clinic$analysis_tbl[c("GSE70768", "GSE70769")], 
         cohort_name = c("GSE70768", "GSE70769"), 
         test_results = coll_clinic$ece_test) %>% 
    future_pmap(plot_clinical_violin, 
                split_factor = 'extra_capsular_extension', 
                variables = coll_clinic$responses, 
                dict = coll_clinic$resp_lexicon, 
                palette = 'Blues', 
                x_lab = 'Extra-capsular extension', 
                .options = furrr_options(seed = TRUE))

# Positive surgical margins and collagen genes ------
  
  insert_msg('Positive surgical margins and collagen genes')
  
  ## descriptive stats
  
  coll_clinic$margin_stats <- 
    coll_clinic$analysis_tbl[c('GSE70768', 'GSE70769')] %>% 
    map(~explore(filter(.x, !is.na(positive_surgical_margins)), 
                 split_factor = 'positive_surgical_margins', 
                 variables = coll_clinic$responses, 
                 what = 'table', 
                 pub_styled = TRUE)) %>% 
    map(reduce, left_join, by = 'variable') %>% 
    map(set_names, c('variable', 'margin_negative', 'margin_positive'))
  
  ## testing
  
  coll_clinic$margin_test <- 
    coll_clinic$analysis_tbl[c('GSE70768', 'GSE70769')] %>% 
    future_map(~compare_variables(.x %>% 
                             filter(!is.na(positive_surgical_margins)) %>% 
                             mutate(positive_surgical_margins = factor(positive_surgical_margins, 
                                                                       c('yes', 'no'))), 
                           split_factor = 'positive_surgical_margins', 
                           variables = coll_clinic$responses, 
                           what = 'eff_size',
                           types = 'cohen_d', 
                           ci = FALSE, 
                           adj_method = 'BH', 
                           pub_styled = TRUE), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(mutate, plot_cap = paste(eff_size, significance, sep = ', '))
  
  ## plotting
  
  coll_clinic$margin_plots <- 
    list(data = coll_clinic$analysis_tbl[c("GSE70768", "GSE70769")], 
         cohort_name = c("GSE70768", "GSE70769"), 
         test_results = coll_clinic$margin_test) %>% 
    future_pmap(plot_clinical_violin, 
                split_factor = 'positive_surgical_margins', 
                variables = coll_clinic$responses, 
                dict = coll_clinic$resp_lexicon, 
                palette = 'Blues', 
                x_lab = 'Positive surgical margins', 
                .options = furrr_options(seed = TRUE))

# caching the results ------
  
  insert_msg('Caching the results')
  
  save(coll_clinic, file = './cache/coll_clinic.RData')  
  
# END -----
  
  rm(i)
  
  plan('sequential')