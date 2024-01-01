# Differences in total mutation burden between the collagen clusters
#
# Statistical significance is assessed by Mann-Whitney test with r effect size 
# statistic

  insert_head()
  
# container ------
  
  ana_tmb <- list()
  
# analysis data ------
  
  insert_msg('Data')
  
  ana_tmb$data <- globals$study_exprs %>% 
    eval %>%
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(safely(select), sample_id, tmb) %>% 
    map(~.x$result) %>% 
    compact
  
  ## cluster assignment 
  
  ana_tmb$data <- 
    map2(ana_globals$assignment[names(ana_tmb$data)], 
         ana_tmb$data, 
         left_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x)))
  
# Numeric stats ------
  
  insert_msg('Numeric stats')
  
  ana_tmb$stats <- ana_tmb$data %>% 
    map(explore, 
        variables = 'tmb', 
        split_factor = 'clust_id',
        what = 'table') %>% 
    map(format_desc) %>% 
    compress(names_to = 'cohort')
  
# Testing --------
  
  insert_msg('Testing')
  
  ana_tmb$test <- ana_tmb$data %>% 
    map(compare_variables, 
        variables = 'tmb', 
        split_factor = 'clust_id', 
        what = 'eff_size', 
        types = 'wilcoxon_r', 
        pub_styled = TRUE, 
        ci = FALSE) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(plot_cap = paste(eff_size, significance))

# Plots --------
  
  insert_msg('Plots')
  
  ana_tmb$plots <- 
    list(x = ana_tmb$data, 
         y = paste('TMB,', globals$study_labels[names(ana_tmb$data)]), 
         z = ana_tmb$test$plot_cap) %>% 
    pmap(function(x, y, z) x %>% 
           plot_variable(variable = 'tmb', 
                         split_factor = 'clust_id', 
                         type = 'box', 
                         cust_theme = globals$common_theme, 
                         plot_title = y, 
                         plot_subtitle = z, 
                         y_lab = 'TMB, arbitrary units', 
                         x_n_labs = TRUE)) %>% 
    map2(., c('log10', 'identity'), 
         ~.x + 
           scale_fill_manual(values = globals$cluster_colors) + 
           scale_y_continuous(trans = .y))
  
# END ------
  
  ana_tmb$data <- NULL
  
  ana_tmb <- compact(ana_tmb)
  
  insert_tail()