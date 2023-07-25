# characteristic of the investigates cohorts

  insert_head()
  
# container ------
  
  cohorts <- list()
  
# globals ------
  
  insert_msg('Globals')
  
  ## variables
  
  cohorts$lexicon <- globals$clinical_lexicon

  ## data set-specific variable lists
  
  cohorts$analysis_tbl <- study_data %>% 
    map(~.x$expression) %>% 
    map(extract_tumor_samples) %>% 
    map(select, any_of(cohorts$lexicon$variable))

  ## analysis tables
  ## wrangling for a common variable format
  ## stages are reduced to the main ones
  
  cohorts$analysis_tbl <- cohorts$analysis_tbl %>% 
    map(format_clinical, main_stage = TRUE)
  
  cohorts$var_list <- cohorts$analysis_tbl %>% 
    map(names)

# Descriptive statistic ------
  
  insert_msg('Descriptive stats')
  
  cohorts$stats <- 
    map2(cohorts$analysis_tbl, 
         cohorts$var_list, 
         ~explore(data = .x, 
                  variables = .y, 
                  what = 'table', 
                  pub_styled = TRUE)) %>% 
    reduce(full_join, by = 'variable') %>% 
    set_names(c('variable', 
                globals$study_labels[names(cohorts$analysis_tbl)])) %>% 
    format_summ_tbl(rm_n = FALSE) %>% 
    mutate(variable = factor(variable, cohorts$lexicon$variable)) %>% 
    arrange(variable) %>% 
    mutate(variable = exchange(variable, cohorts$lexicon))
  
# END -----
  
  insert_tail()