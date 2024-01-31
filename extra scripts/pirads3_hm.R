# Confusion matrix for ISUP2+ detection by the SVM model in PIRADS3 patients

  insert_head()
  
# container -------
  
  extra_pirads <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# Data ------
  
  insert_msg('Data')
  
  ## a data frame in a format accepted by `multiClassSummary()`
  ## 'positive' denotes an ISUP2+ lesion. 
  ## factor levels set in a 'caret manner', i.e. the outcome-of-interest
  ## first
  
  extra_pirads$data$positive <- 
    data.frame(obs = rep('positive', 5), 
               pred = c(rep('positive', 4), 
                        rep('negative', 1)))
  
  extra_pirads$data$negative <- 
    data.frame(obs = rep('negative', 9), 
               pred = c(rep('negative', 8), 
                        rep('positive', 1)))
  
  extra_pirads$data <- extra_pirads$data %>% 
    map_dfr(map_dfc, factor, c('positive', 'negative')) %>% 
    as.data.frame
  
# Stats: general ROC stats, bootstap CI for kappa, sensitivity and specificity -------
  
  insert_msg('Stats')
  
  extra_pirads$stats <- roc_stats(data = extra_pirads$data)
  
  ## bootstrap stats for the sensitivity, specificity and kappa
  
  extra_pirads$ci_stats <- c(Kappa = 'Kappa', 
                             Sensitivity = 'Sensitivity', 
                             Specificity = 'Specificity') %>% 
    future_map(~bmap(extra_pirads$data, 
                     FUN = function(x) roc_stats(x)[[.x]], 
                     B = 1000, 
                     ci_method = 'bca'), 
               .options = furrr_options(seed = TRUE))
  
  extra_pirads$ci_stats <- extra_pirads$ci_stats %>% 
    map(~tibble(statistic_value = .x[['dataset']], 
                lower_ci = .x[['bootstrap']]$boot_lower_ci, 
                upper_ci = .x[['bootstrap']]$boot_upper_ci)) %>% 
    compress(names_to = 'statistic_name') %>% 
    relocate(statistic_name)
  
# Heat map plots -------
  
  insert_msg('Plotting the heat map of the confusion matrix')
  
  ## plotting data
  
  extra_pirads$plot_data <- 
    table(observed = extra_pirads$data$obs, 
          predicted = extra_pirads$data$pred) %>% 
    as.data.frame %>% 
    mutate(observed = factor(observed, c('negative', 'positive')),
           predicted = factor(predicted, c('negative', 'positive')), 
           n = Freq, 
           n_total = sum(Freq), 
           percent = n/n_total * 100, 
           n_lab = paste('n =', n), 
           percent_lab = paste0(signif(percent, 2), '%')) %>%
    blast(observed) %>% 
    map_dfr(mutate, 
            fraction_lab = paste0(n, ' of ', sum(n), '\nobserved')) %>% 
    as_tibble
  
  ## plots
  
  extra_pirads$confusion_hm <- 
    list(x = c('n_lab', 'percent_lab', 'fraction_lab'), 
         y = paste0(c('count', '% of total cases', 'fraction of ISUP strata'), 
                   ', total: n = ', extra_pirads$plot_data$n_total[[1]], 
                   ' PIRADS3 cases')) %>% 
    pmap(function(x, y) extra_pirads$plot_data %>% 
           ggplot(aes(x = observed, 
                      y = predicted, 
                      fill = n)) + 
           geom_tile(color = 'black') + 
           geom_text(aes(label = .data[[x]]), 
                     size = 2.75, 
                     color = 'black') + 
           scale_x_discrete(labels = c(negative = 'ISUP1', 
                                       positive = 'ISUP2+'), 
                            name = 'prostate biopsy') + 
           scale_y_discrete(labels = c(negative = 'ISUP1', 
                                       positive = 'ISUP2+'), 
                            name = 'SVM collagen urinome model') + 
           scale_fill_gradient2(low = 'steelblue', 
                                mid = 'white', 
                                high = 'firebrick') + 
           guides(fill = 'none') + 
           globals$common_theme +
           theme(plot.title.position = 'plot') + 
           labs(title = 'High-risk PCa, PIRADS3 lesions', 
                subtitle = y)) %>% 
    set_names(c('count', 'percentage', 'fraction'))
  
# Saving the results --------
  
  insert_msg('Saving the results')
  
  extra_pirads[c("stats", "ci_stats")] %>% 
    write_xlsx(path = './report/manuscript extras/pirads3_stats.xlsx')
  
  list(plot = extra_pirads$confusion_hm, 
       filename = paste0('./report/manuscript extras/', 
                         c('pirads3_hm_count.pdf', 
                           'pirads3_hm_percent.pdf', 
                           'pirads3_hm_fraction.pdf'))) %>% 
    pwalk(ggsave, 
          width = 8, 
          height = 7, 
          units = 'cm', 
          device = cairo_pdf)
  
# END ------
  
  plan('sequential')

  insert_tail()