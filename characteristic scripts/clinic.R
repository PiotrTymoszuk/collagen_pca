# Clinical characteristic of the collagen clusters. 
# 
# Differences are explored by Chi-squared test or Mann-Whitney test with Cramer's 
# V or r effect size metrics. FDR correction for multiple testing within each 
# collective. 

  insert_head()
  
# container -----
  
  ana_clinic <- list()
  
# Parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis variables -------
  
  insert_msg('Analysis variables')
  
  ana_clinic$variables <- 
    c('age', 'psa_diagnosis', 'ct_stage', 'pt_stage', 'pn_stage', 
      'gleason_simple', 'surgical_margins', 'ece')
  
  ana_clinic$lexicon <- globals$clinical_lexicon %>% 
    filter(variable %in% ana_clinic$variables) %>% 
    mutate(test_type = ifelse(format == 'numeric', 'wilcoxon_r', 'cramer_v'), 
           plot_type = ifelse(format == 'numeric', 'box', 'stack'),
           label = ifelse(variable == 'gleason_simple', 
                          'ISUP', label))

# Analysis data --------
  
  insert_msg('Analysis data')
  
  ana_clinic$data <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(select, 
        sample_id, 
        any_of(ana_clinic$variables))
  
  ## cluster assignment
  ## coding GS as ISUP
  
  ana_clinic$data <- 
    map2(ana_globals$assignment, 
         ana_clinic$data[names(ana_globals$assignment)], 
         left_join, by = 'sample_id') %>% 
    map(filter, !is.na(clust_id))
  
  ana_clinic$data <- ana_clinic$data %>% 
    map(mutate, 
        gleason_simple = car::recode(gleason_simple, 
                                     "'5 - 6' = 'ISUP1'; 
                                     '7' = 'ISUP2'; 
                                     '8' = 'ISUP3+'"), 
        gleason_simple = factor(gleason_simple,  
                                c('ISUP1', 'ISUP2', 'ISUP3+'))) %>%  
    map(map_dfc, function(x) if(is.factor(x)) droplevels(x) else x)

  ## clinical variables present in the given collective
  
  ana_clinic$var_list <- ana_clinic$data %>% 
    map(names) %>% 
    map(~ana_clinic$lexicon$variable[ana_clinic$lexicon$variable %in% .x])
  
# Descriptive stats -------
  
  insert_msg('Descriptive stats')
  
  ana_clinic$stats <- 
    map2(ana_clinic$data, 
         ana_clinic$var_list, 
         explore, 
         split_factor = 'clust_id', 
         what = 'table') %>% 
    map(format_desc)
  
# Testing ---------
  
  insert_msg('Testing')
  
  ana_clinic$test <- 
    list(x = ana_clinic$data, 
         y = ana_clinic$var_list) %>% 
    future_pmap(function(x, y) x %>% 
                  compare_variables(variables = y, 
                                    split_factor = 'clust_id', 
                                    what = 'eff_size', 
                                    types = exchange(y, 
                                                     ana_clinic$lexicon, 
                                                     value = 'test_type'), 
                                    exact = FALSE, 
                                    ci = FALSE, 
                                    pub_styled = TRUE, 
                                    adj_method = 'BH'), 
                .options = furrr_options(seed = TRUE)) %>% 
    map(mutate, 
        plot_cap = paste(eff_size, significance, sep = ', '), 
        plot_cap = ifelse(is.na(significance), '', plot_cap))
  
# Significant factors ---------
  
  insert_msg('Significant variables')
  
  ana_clinic$significant <- ana_clinic$test %>% 
    map(filter, p_adjusted < 0.05) %>% 
    map(~.x$variable)
  
# A common result table -------
  
  insert_msg('A common result table')
  
  ana_clinic$result_tbl <- 
    map2(ana_clinic$stats,
         map(ana_clinic$test, 
             ~.x[c('variable', 'significance', 'eff_size')]), 
         left_join, by = 'variable') %>% 
    map(format_summ_tbl, 
        rm_n = FALSE, 
        rm_complete = TRUE) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(variable = exchange(variable, ana_clinic$lexicon), 
           cohort = globals$study_labels[cohort]) %>% 
    relocate(cohort)
  
# Plot panels for ISUP, pathology stages, PSA and age ------
  
  insert_msg('Plots of ISUP, tumor stages, PSA and age')
  
  ## stack plots of the Gleason scores and umor stages
  
  ana_clinic$plot_panels[c('pt_stage', 'gleason_simple')] <- 
    list(variable = c('pt_stage', 'gleason_simple'), 
         plot_title = c('Pathological tumor stage', 
                        'ISUP')) %>% 
    pmap(plot_stage_stack,
         data_lst = ana_clinic$data, 
         test_lst = ana_clinic$test, 
         show_clust_n = TRUE)
  
  ##PSA at diagnosis and age
  
  ana_clinic$plot_panels[c('age', 'psa_diagnosis')] <- 
    list(variable = c('age', 'psa_diagnosis'), 
         plot_title = c('Age at diagnosis', 'PSA at diagnosis'), 
         x_lab = c('years', 'ng/mL')) %>% 
    pmap(plot_clinic_box, 
         data_lst = ana_clinic$data, 
         test_lst = ana_clinic$test, 
         show_clust_n = TRUE, 
         normalize = FALSE, 
         point_size = 1, 
         point_alpha = 0.5) %>% 
    map2(c('identity', 'sqrt'), 
         ~.x + scale_x_continuous(trans = .y))

# END ------
  
  plan('sequential')
  
  insert_tail()