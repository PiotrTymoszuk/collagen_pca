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

  mri$analysis_tbl <- read_xlsx('./data/MRT Langer.xlsx')
  
  mri$analysis_tbl <- mri$analysis_tbl[1:10, ] %>% 
    transmute(patient_id = as.character(Patientennummer), 
              date = `Datum des MRT`, 
              tumor = `SI Avg Tumor`, 
              benign = `SI Avg kontralat gesunde Seite`) %>% 
    pivot_longer(cols = all_of(c('tumor', 'benign')), 
                 names_to = 'tissue', 
                 values_to = 'adc') %>%
    mutate(tissue = factor(tissue, c('benign', 'tumor')))

# Descriptive stats ------
  
  insert_msg('Descriptive stats')
  
  mri$stats <- 
    explore(mri$analysis_tbl, 
            split_factor = 'tissue', 
            variables = 'adc', 
            what = 'table', 
            pub_styled = TRUE) %>% 
    format_desc
  
# Testing: paired T test --------
  
  insert_msg('Paired T test')
  
  mri$test <- mri$analysis_tbl %>% 
    compare_variables(variables = 'adc', 
                      split_factor = 'tissue', 
                      what = 'eff_size', 
                      types = 'paired_cohen_d', 
                      exact = FALSE, 
                      ci = FALSE, 
                      pub_styled = FALSE) %>% 
    mutate(estimate = -estimate, 
           eff_size = paste('d =', signif(estimate, 2)), 
           plot_cap = paste(eff_size, significance, sep = ', '), 
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
    scale_fill_manual(values = c(tumor = 'coral3', 
                                 benign = 'darkolivegreen4'), 
                      name = '') + 
    theme(axis.title.x = element_blank(), 
          plot.tag = element_blank())
  
# END -----
  
  insert_tail()