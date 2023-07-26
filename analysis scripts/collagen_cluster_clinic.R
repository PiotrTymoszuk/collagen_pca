# Clinical characteristic and values of the Collagen score 
# in the Collagen Clusters

  insert_head()
  
# container -----
  
  cs_cluster <- list()

# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals ------
  
  insert_msg('Globals')
  
  ## analysis variables
  
  cs_cluster$lexicon <- globals$clinical_lexicon %>% 
    filter(!variable %in% c('death', 'relapse', 'vitality_fup', 'relapse_fup')) %>% 
    rbind(tibble(variable = 'collagen_score', 
                 label = 'Collagen Score', 
                 format = 'numeric')) %>% 
    mutate(label_long = label, 
           label = stri_replace(label, regex = ',.*', replacement = ''), 
           unit = ifelse(format == 'factor', '% of cluster', 
                         car::recode(variable, 
                                     "'age' = 'years'; 
                                     'psa_at_diagnosis' = 'AU'; 
                                     'gleason' = 'score'; 
                                     'collagen_score' = 'Collagen Score'")), 
           test_type = ifelse(format == 'numeric', 'kruskal_etasq', 'cramer_v'), 
           plot_type = ifelse(format == 'numeric', 'violin', 'stack'))
  
  ## variables by thier format: the vestor is used later for plot styling
  
  cs_cluster$variables <- cs_cluster$lexicon %>% 
    blast(format) %>% 
    map(~.x$variable)
  
  ## analysis tables
  
  cs_cluster$analysis_tbl$assignment <- coll_clust$assignment
  
  cs_cluster$analysis_tbl$score <- coll_score$score_tbl %>% 
    map(select, patient_id, collagen_score)
  
  cs_cluster$analysis_tbl$clinic <- study_data %>% 
    map(~.x$expression) %>% 
    map(extract_tumor_samples) %>% 
    map(select, patient_id, any_of(cs_cluster$lexicon$variable)) %>% 
    map(format_clinical)
  
  cs_cluster$analysis_tbl <- 
    map2(cs_cluster$analysis_tbl$assignment, 
         cs_cluster$analysis_tbl$score, 
         left_join, by = 'patient_id') %>% 
    map2(cs_cluster$analysis_tbl$clinic, 
         left_join, by = 'patient_id') %>% 
    map(mutate, 
        clust_id = stri_replace(clust_id, fixed = 'Collagen ', replacement = ''), 
        clust_id = factor(clust_id, c('low', 'int', 'hi')), 
        gleason_factor = cut(gleason, 
                             c(-Inf, 7, Inf), 
                             c('6 - 7', '8+')))
  
  ## vectors with variables for single analysis tables
  
  cs_cluster$var_list <- cs_cluster$analysis_tbl %>% 
    map(select, - patient_id, - clust_id) %>% 
    map(names)

# Descriptive stats ------
  
  insert_msg('Descriptive stats')
  
  cs_cluster$stats <- 
    future_map2(cs_cluster$analysis_tbl, 
                cs_cluster$var_list, 
                explore, 
                split_factor = 'clust_id', 
                what = 'table', 
                pub_styled = TRUE, 
                .options = furrr_options(seed = TRUE)) %>% 
    map(reduce, left_join, by = 'variable') %>% 
    map(set_names, c('variable', levels(cs_cluster$analysis_tbl[[1]]$clust_id)))

# Serial testing ------
  
  insert_msg('Serial testing')
  
  cs_cluster$test <- 
    list(x = cs_cluster$analysis_tbl, 
         y = cs_cluster$var_list) %>% 
    future_pmap(function(x, y) x %>% 
                  compare_variables(variables = y, 
                                    split_factor = 'clust_id', 
                                    what = 'eff_size', 
                                    types = exchange(y, 
                                                     cs_cluster$lexicon, 
                                                     value = 'test_type'), 
                                    exact = FALSE, 
                                    ci = FALSE, 
                                    pub_styled = TRUE, 
                                    adj_method = 'BH'), 
                .options = furrr_options(seed = TRUE)) %>% 
    map(mutate, plot_cap = paste(eff_size, significance, sep = ', '))
  
  cs_cluster$signifcant <- cs_cluster$test %>% 
    map(filter, p_adjusted < 0.05) %>% 
    map(~.x$variable)
  
# Result tables ------
  
  insert_msg('Result tables')
  
  cs_cluster$result_tbl <- 
    map2(cs_cluster$stats,
         map(cs_cluster$test, 
             ~.x[c('variable', 'significance', 'eff_size')]), 
         left_join, by = 'variable') %>% 
    map(format_summ_tbl, 
        rm_n = FALSE, 
        rm_complete = FALSE) %>% 
    map(mutate, 
        variable = exchange(variable, 
                            cs_cluster$lexicon, 
                            value = 'label_long'))

# Plots -------
  
  insert_msg('Plots')
  
  for(i in names(cs_cluster$analysis_tbl)) {
    
    cs_cluster$plots[[i]] <- 
      list(variable = cs_cluster$var_list[[i]], 
           type = cs_cluster$var_list[[i]] %>% 
             exchange(cs_cluster$lexicon, 
                      value = 'plot_type'), 
           plot_title = cs_cluster$var_list[[i]] %>% 
             exchange(cs_cluster$lexicon, 
                      value = 'label') %>% 
             paste(globals$study_labels[[i]], 
                   sep = ', '), 
           plot_subtitle = cs_cluster$test[[i]]$plot_cap, 
           y_lab = cs_cluster$var_list[[i]] %>% 
             exchange(cs_cluster$lexicon, 
                      value = 'unit')) %>% 
      pmap(plot_variable, 
           cs_cluster$analysis_tbl[[i]], 
           split_factor = 'clust_id', 
           scale = 'percent', 
           x_lab = 'Collagen Cluster', 
           cust_theme = globals$common_theme, 
           x_n_labs = TRUE, 
           txt_size = 2.5) %>% 
      set_names(cs_cluster$var_list[[i]])
    
    ## styling of the numeric variable plots

    cs_cluster$plots[[i]][cs_cluster$variables$numeric] <- 
      cs_cluster$plots[[i]][cs_cluster$variables$numeric] %>% 
      map(~.x + 
            scale_fill_manual(values = set_names(globals$cluster_colors, 
                                                 c('hi', 'int', 'low'))))
    
    ## styling of the factor variable plots
    
    cs_cluster$plots[[i]][cs_cluster$variables$factor] <- 
      cs_cluster$plots[[i]][cs_cluster$variables$factor] %>% 
      map(~.x + 
            scale_fill_brewer(palette = 'Reds'))
    
    cs_cluster$plots[[i]] <- compact(cs_cluster$plots[[i]])
        
  }

# END -----
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()