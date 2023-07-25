# Normalization and batch effect removal with ComBat 
# in the collagen gene expression dataset
#
# The tables will be further used for survival modeling and 
# semi-supervised clustering

  insert_head()

# container ------
  
  combatch <- list()

# normalized datasets -------
  
  insert_msg('Input datasets')
  
  combatch$variables <- globals$genes_interest$gene_symbol
  
  ## as matrices with genes coded by the rows and samples coded by the columns
  
  combatch$input_mtx <- study_data %>% 
    map(~.x$expression) %>% 
    map(extract_tumor_samples) %>% 
    map(~.x[c('patient_id', combatch$variables)]) %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(column_to_rownames, 'patient_id') %>% 
    map(center_data, type = 'mean')
  
  combatch$input_mtx <- combatch$input_mtx %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = factor(cohort, names(combatch$input_mtx)))
  
  combatch$cohort_assignment <- combatch$input_mtx %>% 
    rownames_to_column('patient_id') %>% 
    select(patient_id, cohort)

  combatch$input_mtx <- combatch$input_mtx %>% 
    select(- cohort) %>% 
    t
  
# Batch effect adjustment -------
  
  insert_msg('Batch effect adjustment')
  
  combatch$adjusted_mtx <- combatch$input_mtx %>% 
    ComBat(batch = combatch$cohort_assignment$cohort)
  
# Retrieving the adjusted dataset --------
  
  insert_msg('Adjusted data frame')
  
  combatch$adjusted_data <- combatch$adjusted_mtx %>% 
    t %>% 
    as.data.frame %>% 
    rownames_to_column('patient_id') %>% 
    left_join(combatch$cohort_assignment, ., by = 'patient_id') %>% 
    as_tibble
  
# Diagnostic stats and tests -------
  
  insert_msg('Diagnostic stats')
  
  ## n numbers
  
  combatch$n_numbers <- combatch$cohort_assignment %>% 
    count(cohort) %>% 
    mutate(cohort = as.character(cohort))
  
  ## descriptive stats
  
  combatch$stats <- combatch$adjusted_data %>% 
    explore(split_factor = 'cohort', 
            variables = combatch$variables, 
            what = 'table', 
            pub_styled = TRUE) %>% 
    reduce(left_join, by = 'variable') %>% 
    set_names(c('variable', levels(combatch$cohort_assignment$cohort)))
  
  ## comparison of the cohorts with one-way ANOVA
  ## and eta-square effect size stat
  
  combatch$test <- combatch$adjusted_data %>% 
    compare_variables(variables = combatch$variables, 
                      split_factor = 'cohort',
                      what = 'eff_size', 
                      types = 'etasq', 
                      exact = FALSE, 
                      ci = FALSE, 
                      pub_styled = TRUE) %>% 
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))
  
  ## result table
  
  combatch$result_tbl <- 
    left_join(combatch$stats, 
              combatch$test[c('variable', 'significance', 'eff_size')], 
              by = 'variable') %>% 
    format_summ_tbl %>% 
    full_rbind(tibble(!!combatch$n_numbers$cohort[1] := combatch$n_numbers$n[1], 
                      !!combatch$n_numbers$cohort[2] := combatch$n_numbers$n[2], 
                      !!combatch$n_numbers$cohort[3] := combatch$n_numbers$n[3], 
                      !!combatch$n_numbers$cohort[4] := combatch$n_numbers$n[4], 
                      !!combatch$n_numbers$cohort[5] := combatch$n_numbers$n[5]), .)
  
# Diagnostic plots -------
  
  insert_msg('Diagnostic plots')
  
  combatch$plots <- 
    list(variable = combatch$variables, 
         plot_title = combatch$variables, 
         plot_subtitle = combatch$test$plot_cap) %>% 
    pmap(plot_variable, 
         combatch$adjusted_data %>% 
           mutate(cohort = car::recode(as.character(cohort), 
                                       "'tcga' = 'TCGA'")), 
         split_factor = 'cohort', 
         type = 'violin', 
         cust_theme = globals$common_theme,
         y_lab = expression('normalized log'[2] * ' expression'), 
         x_n_labs = TRUE) %>% 
    map(~.x + 
          scale_fill_manual(values = c(globals$study_colors, 
                                       TCGA = unname(globals$study_colors['tcga']))) + 
          theme(plot.title = element_text(face = 'bold.italic'))) %>% 
    set_names(combatch$variables)
  
# Splitting the adjusted dataset by the cohort ------
  
  insert_msg('Splitting the adjusted dataset by the cohort')
    
  combatch$adjusted_data <- combatch$adjusted_data %>% 
    blast(cohort) %>% 
    map(select, -cohort)
  
# END -----
  
  insert_tail()