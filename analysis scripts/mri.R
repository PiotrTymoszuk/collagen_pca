# Comparison of the ADC values (MRI) parameter between tumor and benign tissue
# Done for 10 patient pairs.
# At the moment, I'm working with simulated data as a placeholder 
# and waiting for the measuments from the study team

  insert_head()
  
# container ------
  
  mri <- list()
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## analysis table
  ## Warning: at the moment only simulated data with the mean and SD
  ## from the GraphPad plot
  
  mri$analysis_tbl <- 
    tibble(patient_id = paste0('pat_', rep(1:10, 2)), 
           adc = c(rnorm(10, mean = 700, sd = 150), 
                   rnorm(10, mean = 1240, sd = 250)), 
           tissue = c(rep('Tumor', 10), rep('Benign', 10))) %>% 
    mutate(tissue = factor(tissue, c('Tumor', 'Benign')))
  
# Descriptive stats ------
  
  insert_msg('Descriptive stats')
  
  mri$stats <- 
    explore(mri$analysis_tbl, 
            split_factor = 'tissue', 
            variables = 'adc', 
            what = 'table', 
            pub_styled = TRUE) %>% 
    reduce(left_join, by = 'variable') %>% 
    set_names(c('variable', levels(mri$analysis_tbl$tissue)))
  
# Testing: paired T test --------
  
  insert_msg('Paired T test')
  
  mri$test <- mri$analysis_tbl %>% 
    compare_variables(variables = 'adc', 
                      split_factor = 'tissue', 
                      what = 'eff_size', 
                      types = 'paired_cohen_d', 
                      exact = FALSE, 
                      ci = FALSE, 
                      pub_styled = TRUE) %>% 
    mutate(plot_cap = paste(eff_size, significance, sep = ', '), 
           plot_cap = paste0('n = ', n /2, ', ', plot_cap))
  
# Paired-plot ------
  
  insert_msg('Plot')
  
  mri$plot <- mri$analysis_tbl %>% 
    plot_variable(variable = 'adc', 
                  split_factor = 'tissue', 
                  type = 'paired', 
                  cust_theme = globals$common_theme, 
                  plot_title = 'Tissue diffusion capacity', 
                  plot_subtitle = mri$test$plot_cap, 
                  y_lab = expression('mean ADC, 10'^6 * ' mm'^2 * '/s')) + 
    scale_fill_manual(values = c(Tumor = 'coral3', 
                                 Benign = 'darkolivegreen4'), 
                      name = '') + 
    theme(axis.title.x = element_blank(), 
          plot.tag = element_blank())
  
# END -----
  
  insert_tail()