# Gleason score and expression of the collagen-related genes of interest. 
#
# Differential gene expression is assessed by one-way ANOVA and linear modeling 
# with the Gleason 5 - 6 group serving as a baseline. Effect size is measured 
# by eta-square and Cohen's d. 
#
# Differentially regulated genes are defined by the pFDR(ANOVA) < 0.05 with 
# at least weak effect size (eta-square >= 0.02). Shared significant effects 
# are differentially regulated genes in at least five cohorts excluding 
# GSE165060

  insert_head()
  
# container -----
  
  gs_uni <- list()
  
# parallel backend --------
  
  insert_msg('Parallel backend')
  
  plan('multisession')

# analysis globals -------
  
  insert_msg('Analysis globals')
  
  gs_uni$gleason_colors <- c('5 - 6' = 'steelblue', 
                             '7' = 'coral2', 
                             '8+' = 'coral4')
  
# gene expression data -------
  
  insert_msg('Gene expression data')
  
  ## variables
  
  gs_uni$genes <- globals$genes_interest$gene_symbol
  
  ## expression of the collagen-related genes
  
  gs_uni$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(select, 
        sample_id, all_of(gs_uni$genes))
  
  ## Gleason scores
  
  gs_uni$clinic <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(safely(select), 
        sample_id, gleason_simple) %>% 
    map(~.x$result) %>% 
    compact
  
  gs_uni$data <- 
    map2(gs_uni$clinic, 
         gs_uni$expression[names(gs_uni$clinic)], 
         left_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(~mutate(.x, 
                gleason_simple = factor(gleason_simple, 
                                        levels(.x$gleason_simple), 
                                        ordered = TRUE)))
  
# Z scores and mean Z scores --------
  
  insert_msg('Z scores')
  
  gs_uni$z_scores <- gs_uni$data
  
  for(i in names(gs_uni$z_scores)) {
    
    gs_uni$z_scores[[i]][gs_uni$genes] <- 
      gs_uni$z_scores[[i]][gs_uni$genes] %>% 
      center_data('mean')
    
  }
  
  ## strata- and cohort-wise means
  
  gs_uni$mean_z_scores <- gs_uni$z_scores %>% 
    map(select, - sample_id) %>% 
    map(blast, gleason_simple, .skip = TRUE) %>% 
    map(map, colMeans) %>% 
    map(map, 
        compress, 
        names_to = 'gene_symbol', 
        values_to = 'z_score') %>% 
    map(compress, names_to = 'gleason_simple') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(gleason_simple = factor(gleason_simple, 
                                   levels(gs_uni$data[[1]]$gleason_simple)), 
           cohort = factor(cohort, names(gs_uni$data)))
  
# N numbers -------
  
  insert_msg('N numbers')
  
  gs_uni$n_numbers <- gs_uni$data %>% 
    map(count, gleason_simple)
  
  gs_uni$n_tags <- gs_uni$n_number %>% 
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = '))
  
# Descriptive stats -------
  
  insert_msg('Descriptive stats')
  
  gs_uni$stats <- gs_uni$data %>% 
    future_map(~explore(.x, 
                        split_factor = 'gleason_simple', 
                        variables = gs_uni$genes, 
                        what = 'table', 
                        pub_styled = TRUE), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(format_desc)
  
# Testing ---------
  
  insert_msg('Testing')
  
  gs_uni$test <- gs_uni$data %>% 
    future_map(compare_variables, 
               variables = gs_uni$genes, 
               split_factor = 'gleason_simple', 
               what = 'eff_size', 
               types = 'etasq', 
               ci = FALSE, 
               pub_styled = FALSE, 
               .options = furrr_options(seed = TRUE)) %>% 
    map(mutate, 
        eff_size = paste('\u03B7\u00B2 =', signif(estimate, 2)), 
        plot_cap = paste(eff_size, significance, sep = ', '))
  
# Significant effects ------
  
  insert_msg('Significant effects')
  
  ## ANOVA significance
  
  gs_uni$significant <- gs_uni$test %>% 
    map(filter, 
        p_adjusted < 0.05, 
        estimate >= 0.02) %>% 
    map(~.x$variable)
  
  ## common significant effects: shared by at least five cohorts
  ## excluding the gse16560
  
  gs_uni$common_significant <- 
    gs_uni$significant[names(gs_uni$significant) != 'gs16560'] %>% 
    shared_features(m = 5)
  
# Box plots for single variables --------
  
  insert_msg('Box plots for single variables')
  
  for(i in names(gs_uni$test)) {
    
    gs_uni$plots[[i]] <- 
      list(variable = gs_uni$test[[i]]$variable, 
           plot_title = gs_uni$test[[i]]$variable %>% 
             html_italic %>% 
             paste(globals$study_labels[[i]], sep = ', ') %>% 
             html_bold, 
           plot_subtitle = gs_uni$test[[i]]$plot_cap) %>% 
      pmap(plot_variable, 
           gs_uni$data[[i]], 
           split_factor = 'gleason_simple', 
           type = 'box', 
           cust_theme = globals$common_theme, 
           x_n_labs = TRUE) %>% 
      map(~.x + 
            scale_fill_manual(values = gs_uni$gleason_colors) + 
            theme(plot.title = element_markdown())) %>% 
      set_names(gs_uni$test[[i]]$variable)
    
  }
  
  gs_uni$plots <- transpose(gs_uni$plots)
  
# Classification of the genes for specificity for Gleason score strata --------
  
  insert_msg('Classification')
  
  gs_uni$classification <- gs_uni$data %>% 
    map(classify, 
        variables = gs_uni$genes, 
        split_fct = 'gleason_simple')
  
  gs_uni$classification <- transpose(gs_uni$classification)
  
# Metadata for heat map plots --------
  
  insert_msg('Metadata for heat map plots')
  
  ## axis labels with effect sizes and p values
  
  gs_uni$ax_labs <- gs_uni$test %>% 
    map(mutate, 
        plot_cap = paste(html_italic(variable), plot_cap, sep = '<br>'), 
        plot_cap = ifelse(p_adjusted < 0.05 & estimate >= 0.02, 
                          html_bold(plot_cap), plot_cap)) %>% 
    map(~set_names(.x$plot_cap, .x$variable))
  
  gs_uni$short_ax_labs <- gs_uni$test %>% 
    map(mutate, 
        var_lab = html_italic(variable), 
        var_lab = ifelse(p_adjusted < 0.05 & estimate >= 0.02, 
                         html_bold(var_lab), var_lab)) %>% 
    map(~set_names(.x$var_lab, .x$variable))
  
  ## plotting order for heat map of the means
  
  gs_uni$mean_eff_sizes <- gs_uni$test %>% 
    map(select, variable, estimate) %>% 
    compress(names_to = 'cohort') %>% 
    summarise(estimate = mean(estimate), .by = variable)
  
  gs_uni$plot_order <- gs_uni$classification$classification %>% 
    map(left_join, 
        gs_uni$mean_eff_sizes, 
        by = 'variable') %>% 
    compress(names_to = 'cohort') %>% 
    blast(variable) %>% 
    map_dfr(~tibble(gene_symbol = .x$variable[1], 
                    gleason_auc = names(sort(table(.x$gleason_simple), decreasing = TRUE)[1]), 
                    delta_auc = mean(.x$delta_auc), 
                    estimate = .x$estimate[1])) %>% 
    mutate(gleason_auc = factor(gleason_auc, levels(gs_uni$data[[1]]$gleason_simple))) %>% 
    arrange(gleason_auc, estimate)
  
# Ribbon plots -----
  
  insert_msg('Ribbon plots')
  
  gs_uni$ribbon_plots <-
    list(data = gs_uni$z_scores, 
         plot_title = globals$study_labels[names(gs_uni$data)]) %>% 
    pmap(draw_stat_panel, 
         variables = gs_uni$genes, 
         split_factor = 'gleason_simple', 
         stat = 'mean',
         err_stat = '2se', 
         form = 'line', 
         x_lab = 'Z-score, mean \u00B1 2\u00D7SEM', 
         cust_theme = globals$common_theme)
  
  ## styling
  
  gs_uni$ribbon_plots <- 
    list(x = gs_uni$ribbon_plots, 
         y = map(gs_uni$classification$classification, ~.x$variable), 
         v = gs_uni$short_ax_labs, 
         z = gs_uni$n_tags) %>% 
    pmap(function(x, y, v, z) x + 
           scale_fill_manual(values = gs_uni$gleason_colors) + 
           scale_color_manual(values = gs_uni$gleason_colors) + 
           scale_y_discrete(limits = y, 
                            labels = v) + 
           theme(axis.text.y = element_markdown(), 
                 axis.title.y = element_blank()) + 
           labs(subtitle = z))
  
# Heat maps for particular cohorts ------
  
  insert_msg('Heat maps for particular cohorts')
  
  gs_uni$heat_maps <- 
    list(data = gs_uni$z_scores, 
         plot_title = globals$study_labels[names(gs_uni$data)]) %>% 
    pmap(heat_map, 
         variables = gs_uni$genes, 
         split_fct = 'gleason_simple', 
         normalize = FALSE, 
         cust_theme = globals$common_theme, 
         x_lab = 'Cancer sample') %>% 
    map2(., gs_uni$short_ax_labs, 
         ~.x + 
           scale_fill_gradient2(low = 'steelblue', 
                                mid = 'black', 
                                high = 'firebrick', 
                                midpoint = 0, 
                                limits = c(-3, 3), 
                                oob = scales::squish, 
                                name = 'Z-scores') + 
           scale_y_discrete(labels = .y) + 
           theme(axis.title.y = element_blank(), 
                 axis.text.y = element_markdown(), 
                 strip.background.y = element_blank(), 
                 strip.text.y = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.x = element_blank(), 
                 axis.line.x = element_blank()))
  
# Heat map of cohort-wise mean Z scores -------
  
  insert_msg('Heat map of the means')
  
  ## for the common significant genes
  
  gs_uni$mean_heat_map <- gs_uni$mean_z_scores %>% 
    left_join(gs_uni$plot_order, by = 'gene_symbol') %>% 
    filter(gene_symbol %in% gs_uni$common_significant) %>% 
    ggplot(aes(x = cohort, 
               y = reorder(gene_symbol, estimate), 
               fill = z_score)) + 
    geom_tile() + 
    facet_grid(gleason_auc ~ gleason_simple, 
               scales = 'free', 
               space = 'free') + 
    scale_fill_gradient2(low = 'steelblue', 
                         mid = 'black', 
                         high = 'firebrick', 
                         midpoint = 0, 
                         limits = c(-1, 1), 
                         oob = scales::squish, 
                         name = 'Mean Z score') + 
    scale_x_discrete(labels = globals$study_labels) + 
    globals$common_theme + 
    theme(axis.title = element_blank(), 
          axis.text.y = element_text(face = 'italic'), 
          strip.background.y = element_blank(), 
          strip.text.y = element_blank()) + 
    labs(title = 'Collagen genes and Gleason scores', 
         subtitle = 'average Z scores') + 
    guides(x = guide_axis(angle = 90))
  
# Result table for the manuscript --------
  
  insert_msg('Result table for the manuscript')
  
  gs_uni$result_tbl <- 
    map2(gs_uni$stats, 
         map(gs_uni$test, ~.x[c('variable', 'significance', 'eff_size')]), 
         left_join, by = 'variable') %>% 
    map(format_summ_tbl) %>% 
    map2(., gs_uni$n_numbers, 
         ~full_rbind(tibble(variable = 'Samples, n', 
                            !!as.character(.y[[1]][1]) := .y[[2]][1], 
                            !!as.character(.y[[1]][2]) := .y[[2]][2],
                            !!as.character(.y[[1]][3]) := .y[[2]][3]), 
                     .x))
  
  gs_uni$result_tbl <- gs_uni$result_tbl %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    relocate(cohort) %>% 
    set_names(c('Cohort', 'Variable', 
                levels(gs_uni$data[[1]]$gleason_simple), 
                'Significance', 'Effect size'))
  
# END -------
  
  rm(i)

  gs_uni <- gs_uni[c("genes", "n_numbers", "stats", "test", 
                     "significant", "common_significant", 
                     "plots", "classification", "ribbon_plots", 
                     "heat_maps", "mean_heat_map", "result_tbl")]
  
  plan('sequential')
  
  insert_tail()