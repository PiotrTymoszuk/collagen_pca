# Normalization and batch effect removal with ComBat 
# done for the collagen-related genes
#
# The tables may be further used for survival modeling and 
# semi-supervised clustering

  insert_head()

# container ------
  
  combat <- list()

# normalized data sets -------
  
  insert_msg('Input datasets')
  
  combat$variables <- globals$genes_interest$gene_symbol
  
  ## as matrices with genes coded by the rows and samples coded by the columns
  
  combat$input_mtx <- globals$study_exprs %>% 
    eval %>%
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(~.x[c('sample_id', combat$variables)]) %>% 
    map(column_to_rownames, 'sample_id')

  combat$input_mtx <- combat$input_mtx %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = factor(cohort, names(combat$input_mtx)))
  
  combat$cohort_assignment <- combat$input_mtx %>% 
    rownames_to_column('sample_id') %>% 
    select(sample_id, cohort)

  combat$input_mtx <- combat$input_mtx %>% 
    select(- cohort) %>% 
    t
  
# Batch effect adjustment -------
  
  insert_msg('Batch effect adjustment')
  
  combat$adjusted_mtx <- combat$input_mtx %>% 
    ComBat(batch = combat$cohort_assignment$cohort)
  
# Retrieving the adjusted dataset --------
  
  insert_msg('Adjusted data frame')
  
  combat$adjusted_data <- combat$adjusted_mtx %>% 
    t %>% 
    as.data.frame %>% 
    rownames_to_column('sample_id') %>% 
    left_join(combat$cohort_assignment, ., by = 'sample_id') %>% 
    as_tibble
  
# Diagnostic stats and tests -------
  
  insert_msg('Diagnostic stats')
  
  ## n numbers
  
  combat$n_numbers <- combat$cohort_assignment %>% 
    count(cohort) %>% 
    mutate(cohort = as.character(cohort))
  
  ## descriptive stats
  
  combat$stats <- combat$adjusted_data %>% 
    explore(split_factor = 'cohort', 
            variables = combat$variables, 
            what = 'table', 
            pub_styled = TRUE) %>% 
    reduce(left_join, by = 'variable') %>% 
    set_names(c('variable', levels(combat$cohort_assignment$cohort)))
  
  ## comparison of the cohorts with one-way ANOVA
  ## and eta-square effect size stat
  
  combat$test <- combat$adjusted_data %>% 
    compare_variables(variables = combat$variables, 
                      split_factor = 'cohort',
                      what = 'eff_size', 
                      types = 'etasq', 
                      exact = FALSE, 
                      ci = FALSE, 
                      pub_styled = TRUE, 
                      .parallel = TRUE, 
                      .paropts = furrr_options(seed = TRUE, 
                                               globals = c('combat'))) %>% 
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))
  
  ## result table
  
  combat$result_tbl <- 
    left_join(combat$stats, 
              combat$test[c('variable', 'significance', 'eff_size')], 
              by = 'variable') %>% 
    format_summ_tbl %>% 
    full_rbind(tibble(!!combat$n_numbers$cohort[1] := combat$n_numbers$n[1], 
                      !!combat$n_numbers$cohort[2] := combat$n_numbers$n[2], 
                      !!combat$n_numbers$cohort[3] := combat$n_numbers$n[3], 
                      !!combat$n_numbers$cohort[4] := combat$n_numbers$n[4], 
                      !!combat$n_numbers$cohort[5] := combat$n_numbers$n[5], 
                      !!combat$n_numbers$cohort[6] := combat$n_numbers$n[6], 
                      !!combat$n_numbers$cohort[7] := combat$n_numbers$n[7]), .)
  
# Diagnostic plots -------
  
  insert_msg('Diagnostic plots')
  
  plan('multisession')
  
  combat$plots <- 
    list(variable = combat$variables, 
         plot_title = combat$variables, 
         plot_subtitle = combat$test$plot_cap) %>% 
    future_pmap(plot_variable, 
                combat$adjusted_data %>% 
                  mutate(cohort = unname(globals$study_labels[cohort])), 
                split_factor = 'cohort', 
                type = 'violin', 
                cust_theme = globals$common_theme,
                y_lab = expression('normalized log'[2] * ' expression'), 
                x_n_labs = TRUE, 
                .options = furrr_options(seed = TRUE)) %>% 
    map(~.x + 
          scale_fill_manual(values = c(globals$study_colors, 
                                       TCGA = unname(globals$study_colors['tcga']))) + 
          theme(plot.title = element_text(face = 'bold.italic'))) %>% 
    set_names(combat$variables)
  
# Splitting the adjusted dataset by the cohort ------
  
  insert_msg('Splitting the adjusted dataset by the cohort')
    
  combat$adjusted_data <- combat$adjusted_data %>% 
    blast(cohort) %>% 
    map(select, -cohort)
  
# END -----
  
  plan('sequential')
  
  insert_tail()