# Spearman's correlation between the Collagen Score 
# and non-malignant cell content

  insert_head()
  
# container ------
  
  cs_infil <- list()
  
# Parallel backend ----
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals -------
  
  insert_msg('Globals')
  
  ## variables
  
  cs_infil$variables <- infil[c("xcell_types", "mcp_counter_types")] %>% 
    set_names(c('xcell', 'mcp'))
  
  ## analysis tables
  
  cs_infil$analysis_tbl[c('xcell', 'mcp')] <- 
    infil[c("xcell", "mcp_counter")] %>% 
    map(function(algo) map2(algo, 
                            map(coll_score$score_tbl, 
                                ~.x[c('patient_id', 'collagen_score')]), 
                            inner_join, 
                            by = 'patient_id'))
  
  cs_infil$analysis_tbl <- cs_infil$analysis_tbl %>% 
    map(map, ~filter(.x, complete.cases(.x)))

  ## n numbers
  
  cs_infil$n_numbers <- cs_infil$analysis_tbl %>% 
    map(map, nrow)

  ## cohort labels
  
  cs_infil$n_tags <- cs_infil$n_numbers %>% 
    map(~map2_chr(.x, names(.x), 
                  ~paste(globals$study_labels[.y], .x, sep = '\nn = ')))
  
  ## algorithm labels
  
  cs_infil$algo_labs <- c(xcell = 'xCell', 
                          mcp = 'MCP counter')

# Serial correlation analysis ------
  
  insert_msg('Serial correlation analysis')
  
  for(i in names(cs_infil$analysis_tbl)) {
    
    cs_infil$test[[i]] <- cs_infil$analysis_tbl[[i]] %>% 
      future_map(function(data)  map(cs_infil$variables[[i]], 
                                     ~c('collagen_score', .x)) %>% 
                   map(~safely(correlate_variables)(data, 
                                                    variables = .x, 
                                                    what = 'correlation', 
                                                    type = 'spearman', 
                                                    ci = TRUE, 
                                                    pub_styled = FALSE)) %>% 
                   map_dfr(~.x$result), 
                 .options = furrr_options(seed = TRUE)) %>% 
      map(re_adjust) %>% 
      map(mutate, 
          correlation = ifelse(p_adjusted >= 0.05, 
                               'ns', 
                               ifelse(estimate > 0, 'positive', 'negative')), 
          correlation = factor(correlation, c('positive', 'negative', 'ns')), 
          plot_cap = paste0('\u03C1 = ', signif(estimate, 2), 
                            ' [', signif(lower_ci, 2), 
                            ' - ', signif(upper_ci, 2), ']'), 
          plot_cap = paste(plot_cap, significance, sep = ', '), 
          plot_cap = paste(plot_cap, n, sep = ', n = '))
    
  }

# Bubble plots with the significant correlation coefficients ------
  
  insert_msg('Correlograms')
  
  ## plotting tables, only significant correlations
  
  cs_infil$bubble$data <- cs_infil$test %>% 
    map(compress, names_to = 'cohort') %>% 
    map(mutate,
        font_face = ifelse(correlation == 'ns', 'plain', 'bold'), 
        significance = ifelse(correlation == 'ns', 'no', 'yes'))

  ## plots
  
  cs_infil$bubble$plots <- 
    list(data = cs_infil$bubble$data, 
         plot_title = paste('Collagen score and infiltration,', 
                            c('xCell', 'MCP Counter')), 
         x_scale = cs_infil$n_tags) %>% 
    pmap(function(data, plot_title, x_scale) data %>% 
           ggplot(aes(x = cohort, 
                      y = variable2, 
                      size = abs(estimate), 
                      fill = estimate)) + 
           geom_point(shape = 21) + 
           geom_text(aes(label = signif(estimate, 2), 
                         alpha = significance, 
                         fontface = font_face), 
                     size = 2.5, 
                     hjust = 0.5, 
                     vjust = -1.6) + 
           scale_alpha_manual(values = c(no = 0.35, 
                                         yes = 1)) + 
           scale_fill_gradient2(low = 'steelblue', 
                                high = 'firebrick', 
                                mid = 'white', 
                                midpoint = 0, 
                                limits = c(-0.6, 0.6),
                                breaks = round(seq(-0.6, 0.6, by = 0.2), 2), 
                                name = expression(rho)) + 
           scale_radius(limits = c(0, 0.6), 
                        breaks = round(seq(0, 0.6, by = 0.1), 2), 
                        name = expression('abs(' * rho * ')')) +
           scale_x_discrete(labels = x_scale) + 
           guides(alpha = 'none') + 
           globals$common_theme + 
           theme(axis.title = element_blank())  +
           labs(title = plot_title, 
                subtitle = 'Spearman correlation'))
  
# Scatter plots for cancer associated fibroblasts -------
  
  insert_msg('Correlation plots for cancer assocaited fibroblasts')
  
  for(i in names(cs_infil$analysis_tbl)) {
    
    cs_infil$caf_plots[[i]] <- 
      list(data = cs_infil$analysis_tbl[[i]], 
           plot_title = globals$study_labels[names(cs_infil$analysis_tbl[[i]])], 
           plot_subtitle = cs_infil$test[[i]] %>% 
             map(filter, variable2 == 'Cancer associated fibroblast') %>% 
             map(~.x$plot_cap), 
           point_color = globals$study_colors[names(cs_infil$analysis_tbl[[i]])]) %>% 
      pmap(plot_correlation, 
           variables = c('collagen_score', 
                         'Cancer associated fibroblast'), 
           type = 'correlation', 
           cust_theme = globals$common_theme, 
           x_lab = 'Collagen Score', 
           y_lab = paste('Cancer-associated fibroblasts', 
                         cs_infil$algo_labs[[i]])) %>% 
      map(~.x + theme(plot.tag = element_blank()))
    
  }

# END ------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()