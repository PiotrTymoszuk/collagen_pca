# This script provides project specific tools ----

  library(plyr)
  library(tidyverse)
  library(ggvenn)
  library(cowplot)
  library(UpSetR)
  library(furrr)
  library(caret)
  library(doParallel)
  library(survminer)
  library(exda)
  library(furrr)
  library(stringi)
  library(rlang)
  library(rstatix)

# Clinical data formatting -------

  extract_pfield <- function(x, 
                             format = c('numeric', 'factor', 'character'), 
                             levels = NULL) {
    
    format <- match.arg(format[1], 
                        c('numeric', 'factor', 'character'))

    if(is.null(levels)) {
      
      levels <- sort(unique(as.character(x)))
      
    }
        
    x <- stri_replace(x, regex = '.*:\\s{1}', replacement = '')
    
    switch(format, 
           numeric = as.numeric(x), 
           factor = factor(x, levels = levels), 
           character = as.character(x))
    
    
  }

  regex_replacer <- function(x, regex, replacement = '') {
    
    if(is.character(x)) {
      
      return(stri_replace(x, regex = regex, replacement = replacement))
      
    } else {
      
      return(x)
      
    }
    
  }

  format_stage <- function(data, main_stage = TRUE) {
    
    ## formats stage information
    
    data <- data %>% 
      map_dfc(regex_replacer, regex = '(T|N|M)x', replacement = NA) %>% 
      map_dfc(regex_replacer, regex = '(T|N|M)X', replacement = NA) %>% 
      map_dfc(regex_replacer, regex = '(p|pT|N|M)NA', replacement = NA)
    
    for(i in names(data)) {
      
      if(is.numeric(data[[i]]) | is.factor(data[[i]])) next
      
      data <- data %>% 
        mutate(!!i := ifelse(stri_detect(.data[[i]], 
                                         regex = '(T|N|M)\\d{1}'), 
                             regex_replacer(.data[[i]], 
                                            regex = '^p',
                                            replacement = ''), 
                             .data[[i]]))
      
    }

    if(main_stage) {
      
      for(i in names(data)) {
        
        if(is.numeric(data[[i]]) | is.factor(data[[i]])) next
        
        data <- data %>% 
          mutate(!!i := ifelse(stri_detect(.data[[i]], 
                                           regex = '(T|N|M)\\d{1}[a-d]'), 
                               regex_replacer(.data[[i]], 
                                              regex = '[a-d]$', 
                                              replacement = ''), 
                               .data[[i]]))
        
      }

    }
    
    return(data)
    
  }

  format_clinical <- function(data, main_stage = TRUE) {
    
    ## takes care for the proper format of clinical variables
    
    data <- format_stage(data, main_stage = main_stage)
    
    if('psa_at_diagnosis' %in% names(data)) {
      
      data <- mutate(data, psa_at_diagnosis = as.numeric(psa_at_diagnosis))
      
    }

    if('age' %in% names(data)) {
      
      data <- mutate(data, age = as.numeric(age))
      
    }
    
    if('gleason' %in% names(data)) {
      
      data <- data %>% 
        mutate(gleason = as.numeric(gleason), 
               gleason_factor = factor(gleason))
      
    }
    
    if('vitality_fup' %in% names(data)) {
      
      data <- mutate(data, vitality_fup = as.numeric(vitality_fup))
      
    }
    
    if('relapse_fup' %in% names(data)) {
      
      data <- mutate(data, relapse_fup = as.numeric(relapse_fup))
      
    }
    
    if('death' %in% names(data)) {
      
      data <- data %>% 
        mutate(death = as.character(death), 
               death = car::recode(death, "'1' = 'yes'; '0' = 'no'"))
      
    }
    
    if('relapse' %in% names(data)) {
      
      data <- data %>% 
        mutate(relapse = as.character(relapse), 
               relapse = car::recode(relapse, "'1' = 'yes'; '0' = 'no'"))
      
    }
    
    data %>% 
      map_dfc(function(x) if(is.character(x)) factor(x) else x)
    
  }
  
  extract_tumor_samples <- function(data) {
    
    if('tissue' %in% names(data)) {
      
      return(filter(data, tissue == 'tumor'))
      
    } else {
      
      return(data)
    }
    
    
  }
  
# Venn and Upset plotting -----

  plot_venn <- function(plotting_lst,  
                        text = NULL, 
                        show_text = TRUE, 
                        colors = c('blue', 'yellow', 'green', 'red'), 
                        plot_title = NULL, 
                        plot_subtitle = NULL, 
                        short = TRUE, 
                        fct_per_line = 4, 
                        rel_widths = c(0.4, 0.6), ...) {
    
    ## generates a Venn plot to show an overlap between the significant features
    ## identified for different risk responses
    
    ## plot
    
    venn_plot <- plotting_lst %>% 
      ggvenn(show_percentage = FALSE, 
             fill_color = unname(colors), 
             set_name_size = 2.75, 
             text_size = 2.75) + 
      theme(plot.title = element_text(size = 8, face = 'bold'), 
            plot.subtitle = globals$common_text) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle)
    
    if(!show_text) return(venn_plot)
    
    ## feature listing
    
    if(is.null(text)) {
      
      feat_txt <- plotting_lst %>% 
        reduce(intersect) %>% 
        sort
      
    } else {
      
      feat_txt <- text
      
    }
    
    feat_txt <- feat_txt %>% 
      wrap_vector(line_length = fct_per_line)
    
    
    ## plot panel
    
    venn_panel <- plot_grid(venn_plot, 
                            ggdraw() + 
                              draw_text(feat_txt, 
                                        hjust = 0, 
                                        x = 0.05, 
                                        size = 8, ...) + 
                              theme(plot.margin = ggplot2::margin(t = 0, b = 0, r = 2, l = 3)), 
                            ncol = 2, 
                            rel_widths = rel_widths) + 
      theme(plot.margin = globals$common_margin)
    
    return(venn_panel)
    
  }
  
  plot_upset <- function(plotting_lst, 
                         plot_title = '', 
                         plot_subtitle = '', 
                         fct_per_line = 4, 
                         label_common = TRUE, 
                         query = NULL, 
                         feat_txt = NULL, 
                         y_lab = '# common variables', 
                         rel_widths = c(2, 1), ...) {
    
    ## makes an upset plot to visualize intersections
    
    if(label_common) {
      
      query <- list(list(query = intersects, 
                         params = list(names(plotting_lst)), 
                         color = 'coral3', 
                         active = T))
      
    }
    
    ## upset
    
    upset_plot <- upset(fromList(plotting_lst), 
                        order.by = 'freq', 
                        nsets = 5, 
                        set_size.show = FALSE, 
                        queries = query, 
                        main.bar.color = 'gray40', 
                        mainbar.y.label = y_lab, 
                        text.scale = 1.2)
    
    ## listing of the common factors
    
    if(is.null(feat_txt)) {
      
      feat_txt <- plotting_lst %>% 
        reduce(intersect) %>% 
        stri_sort %>% 
        wrap_vector(line_length = fct_per_line)
      
    }
    
    ## panel
    
    upset_panel <- plot_grid(upset_plot$Main_bar, 
                             upset_plot$Matrix, 
                             nrow = 2, 
                             align = 'hv', 
                             rel_heights = c(1.5, 0.5)) %>% 
      plot_grid(., 
                ggdraw() + 
                  draw_text(feat_txt, 
                            hjust = 0, 
                            x = 0.05, 
                            size = 8, ...) + 
                  theme(plot.margin = ggplot2::margin(t = 0, b = 0, r = 2, l = 3)), 
                ncol = 2, 
                rel_widths = rel_widths) %>% 
      plot_grid(ggdraw() + 
                  draw_text(plot_subtitle, 
                            size = 9, 
                            fontface = 'plain', 
                            hjust = 0, 
                            x = 0.05), 
                ., 
                nrow = 2, 
                rel_heights = c(0.05, 0.95)) %>% 
      plot_grid(ggdraw() + 
                  draw_text(plot_title, 
                            size = 9, 
                            fontface = 'bold', 
                            hjust = 0, 
                            x = 0.05), 
                ., 
                nrow = 2, 
                rel_heights = c(0.05, 0.95))
    
    upset_panel <- upset_panel +
      theme(plot.margin = globals$common_margin)
    
    return(upset_panel)
    
  }
  
# Text functions -----
  
  split_vec <- function(inp_vector, chunk_size) {
    
    return(split(inp_vector, ceiling(seq_along(inp_vector)/chunk_size)))
    
  }
  
  wrap_vector <- function(txt_vec, line_length = 5) {
    
    split_txt <- split_vec(txt_vec, 
                           line_length) %>% 
      map(paste, 
          collapse = ', ') %>%
      paste(collapse = ',\n')
    
    return(split_txt)
    
  }
  
  capitalize_first_char <- function(x) {
    
    ## capitalizes the first character of a string
    
    stri_trans_totitle(stri_sub(x, 1, 1)) %s+% stri_sub(x, 2)
    
  }
  
  space2break <- function(x, n) {
    
    ## n defines the n-th space to be converted to a line break
    
    if(length(x) > 1) {
      
      stop("'x' must be a single string, not a vector of strings.", 
           call. = FALSE)
      
    }
    
    split_string <- stri_split_fixed(x, " ")[[1]]
    
    if(length(split_string) < n + 1) return(x)

    paste(paste(split_string[1:n], 
                collapse = ' '), 
          paste(split_string[(n + 1):length(split_string)], 
                collapse = ' '), 
          sep = '\n')
    
  }
  
  calc_string_dist <- function(x) {
    
    ## calculates pairwise normalized Levenshtein distances
    ## between pairs of the list elements
    
    els <- names(x)
    
    res <- matrix(data = NA, nrow = length(x), ncol = length(x))
    
    res <- set_colnames.matrix(res, els)
    res <- set_rownames.matrix(res, els)
    
    for(i in els) {
      
      res[i, ] <- els %>% 
        future_map_dbl(~StrDist(as.character(as.numeric(factor(x[[i]]))), 
                                as.character(as.numeric(factor(x[[.x]]))), 
                                method = 'normlevenshtein'))
      
    }
    
    return(1 - res)
    
  }
  
  gsva_labeller <- function(x) {
    
    repl <- 
      c('EXTRACELLULAR MATRIX' = 'ECM', 
        'INTERLEUKIN ' = 'IL', 
        'INSULIN LIKE GROWTH FACTOR BINDING PROTEINS ' = '', 
        'INSULIN LIKE GROWTH FACTOR ' = '')
    
    for(i in names(repl)) {
      
      x <- stri_replace(x, fixed = i, replacement = repl[i])
      
    }
    
    x
    
  }
  
  biol_labeller <- function(x) {
    
    x <- stri_replace_all(x, fixed = '\n', replacement = ' ')
    
    repl <- c('ECM-receptor interaction' = 'ECM', 
              'Small cell lung cancer' = 'SCLC', 
              'Staphylococcus aureus' = 'S. aureus', 
              'Systemic lupus erythematosus' = 'SLE', 
              'Fatty acid oxidation' = 'FAOX', 
              'isoleucine' = 'Ile', 
              'Valine' = 'Val', 
              'leucine' = 'Leu', 
              'and ' = '', 
              ' infection' = '\ninfection', 
              ' metabolism' = '\nmetabolism', 
              ' synthesis' = '\nsynthesis', 
              'Phosphatidylinositol phosphate' = 'PIP', 
              'interconversion' = '\ninterconversion', 
              'Transport, extracellular' = 'Extracellular\ntransport', 
              'adhesion' = '\nadhesion', 
              'Regulation of actin cytoskeleton' = 'ACT\ncytoskeleton', 
              'TGF-beta signaling pathway' = 'TGF-\u03B2', 
              ' in' = '\nin')
    
    for(i in names(repl)) {
      
      x <- stri_replace_all(x, fixed = i, replacement = repl[[i]])
      
    }
    
    x
    
  }
  
# Result formatting --------
  
  format_summ_tbl <- function(data, 
                              rm_n = TRUE, 
                              rm_mean = TRUE, 
                              rm_complete = TRUE) {
    
    ## formats a summary table with descriptive stats
    
    data <- data %>% 
      map_dfc(stri_replace, regex = 'no:.*\\nyes:\\s{1}', replacement = '') %>% 
      map_dfc(stri_replace, regex = '\\nno:.*$', replacement = '\n') %>% 
      map_dfc(stri_replace_all, fixed = '% (', replacement = '% (n = ') %>% 
      map_dfc(stri_replace, fixed = 'Median =', replacement = 'median:') %>% 
      map_dfc(stri_replace, fixed = 'Mean =', replacement = 'mean:') %>% 
      map_dfc(stri_replace, fixed = 'Range', replacement = 'range') %>% 
      map_dfc(stri_replace, fixed = 'Complete', replacement = 'complete') %>% 
      map_dfc(stri_replace, fixed = ' [', replacement = '\n[')
    
    if(rm_n) {
      
      data <- data %>% 
        map_dfc(stri_replace, regex = '\\ncomplete.*$', replacement = '')
      
    }
    
    if(rm_mean) {
      
      data <- data %>% 
        map_dfc(stri_replace, 
                regex = 'mean.*\\nmedian:\\s{1}', 
                replacement = '') %>% 
        map_dfc(stri_replace, fixed = '\n[', replacement = ' [')
      
    }
    
    if(rm_complete) {
      
      data <- data %>% 
        map_dfc(stri_replace, fixed = 'complete: ', replacement = '')
      
    }
    
    data
    
  }
  
  re_adjust <- function(data, p_variable = 'p_value', adj_method = 'BH') {
    
    data %>% 
      mutate(p_adjusted = p.adjust(.data[[p_variable]], method = adj_method), 
             significance = ifelse(p_adjusted < 0.05, 
                                   paste('p =', signif(p_adjusted, 2)), 
                                   paste0('ns (p = ', signif(p_adjusted, 2), ')')))
    
  }
  
# Correlation result formatting --------
  
  format_corr_results <- function(data, cohort_name, dict, ...) {
    
    ## formats results of correlation of gene expression and a clinical factor
    
    data %>% 
      mutate(cohort = cohort_name, 
             eff_size = stri_replace(eff_size, 
                                     fixed = 'rho', 
                                     replacement = '\u03C1'), 
             plot_cap = paste(eff_size, significance, sep = ', '), 
             plot_cap = paste(plot_cap, n, sep = ', n = '), 
             plot_title = ifelse(variable2 != 'collagen_score', 
                                 paste0('<b><em>', 
                                        exchange(variable2, dict, ...), 
                                        '</em>, ', 
                                        globals$study_labels[cohort_name], 
                                        '</b>'), 
                                 paste0('<b>', 
                                        exchange(variable2, dict, ...), 
                                        ', ', 
                                        globals$study_labels[cohort_name], 
                                        '</b>')), 
             y_lab = ifelse(variable2 != 'collagen_score', 
                            paste0('<em>', 
                                   exchange(variable2, dict), 
                                   '</em>, log<sub>2</sub> expression'), 
                            exchange(variable2, dict)))
    
  }
  
# Clinical testing result plotting, collagen genes and score -------
  
  plot_clinical_correlation <- function(data, 
                                        test_results, 
                                        point_color = 'steelblue') {
    
    ## draws a series of scatter plots visualizing correlation between
    ## variables
    
    list(variables = map2(test_results$variable1, 
                          test_results$variable2, c), 
         plot_title = test_results$plot_title, 
         plot_subtitle = test_results$plot_cap, 
         y_lab = test_results$y_lab, 
         x_lab = test_results$x_lab) %>% 
      pmap(function(variables, 
                    plot_title,
                    plot_subtitle, 
                    y_lab, 
                    x_lab) data[variables] %>% 
             filter(complete.cases(.)) %>% 
             plot_correlation(variables = variables, 
                              type = 'correlation', 
                              point_color = point_color, 
                              plot_title = plot_title, 
                              plot_subtitle = plot_subtitle, 
                              y_lab = y_lab, 
                              x_lab = x_lab, 
                              cust_theme = globals$common_theme)) %>% 
      map(~.x + 
            theme(plot.tag = element_blank(), 
                  axis.title.y = element_markdown(), 
                  plot.title = element_markdown())) %>% 
      set_names(map2_chr(test_results$variable1, 
                         test_results$variable2, 
                         paste, sep = '_'))
    
  }
  
  plot_clinical_violin <- function(data, 
                                   split_factor, 
                                   variables, 
                                   cohort_name, 
                                   test_results, 
                                   dict = coll_clinic$resp_labs, 
                                   palette = 'Blues', 
                                   x_lab = NULL) {
    
    ## plots a series of violin plots
    
    data <- data %>% 
      filter(!is.na(.data[[split_factor]]))
    
    list(variable = variables, 
         plot_title = ifelse(variables == 'collagen_score', 
                             paste('<b>Collagen Score, ', 
                                   globals$study_labels[[cohort_name]], 
                                   '</b>'), 
                             paste('<b><em>', 
                                   exchange(variables, dict), 
                                   '</em>, ', 
                                   globals$study_labels[[cohort_name]], 
                                   '</b>')), 
         plot_subtitle = test_results$plot_cap, 
         y_lab = ifelse(variables == 'collage_score', 
                        'Collagen Score', 
                        paste('<em>', 
                              exchange(variables, dict), 
                              '</em>, log<sub>2</sub> expression'))) %>% 
      future_pmap(plot_variable, 
                  data, 
                  split_factor = split_factor, 
                  type = 'violin', 
                  cust_theme = globals$common_theme, 
                  x_lab = x_lab, 
                  x_n_labs = TRUE, 
                  .options = furrr_options(seed = TRUE)) %>% 
      map(~.x + 
            scale_fill_brewer(palette = palette) + 
            theme(plot.title = element_markdown(), 
                  axis.title.y = element_markdown())) %>% 
      set_names(variables)
    
  }
  
# Semi-supervised clustering -----
  
  adapt_knn <- function(object, newdata, max_kNN = 15, ...) {
    
    ## finds the optimal number of nearest neighbors for semi-supervised
    ## clustering based on the maximal clustering variance
    
    tst_k <- seq(3, max_kNN, by = 2)
    
    tst_k <- set_names(tst_k, paste0('knn_', tst_k))
    
    tst_objects <- list(kNN = tst_k) %>% 
      future_pmap(object = object, 
                  newdata = newdata, 
                  predict.clust_analysis, ..., 
                  .options = furrr_options(seed = TRUE))
    
    variances <- tst_objects %>%
      map(clustTools::var) %>% 
      map_dbl(~.x$frac_var) %>% 
      compress(names_to = 'model', 
               values_to = 'clust_variance') %>% 
      mutate(kNN = unname(tst_k), 
             optimum = ifelse(clust_variance == max(clust_variance), 'yes', 'no'))
    
    best_object <- variances %>% 
      filter(optimum == 'yes')
    
    best_object <- best_object[1, ]
    
    variance_plot <- variances %>% 
      ggplot(aes(x = kNN, 
                 y = clust_variance)) + 
      geom_vline(xintercept = best_object[['kNN']], 
                 linetype = 'dashed') + 
      geom_path(color = 'steelblue') + 
      globals$common_theme + 
      labs(title = 'kNN classifier tuning', 
           subtitle = 'Semi-supervised slustering', 
           y = 'Clusetring variance', 
           x = 'kNN')
    
    list(object = tst_objects[[best_object[['model']]]], 
         tuning_stats = variances, 
         kNN = best_object[['kNN']], 
         plot = variance_plot)
    
  }
  
  
# GO similarity -----
  
  go_similarity <- function(GOs, semData, dist = TRUE, ...) {
    
    ## computes a similarity or distance matrix for a vector of GO terms
    ## it will be used for representation of related GOs in plots
    
    mtx <- matrix(data = NA, nrow = length(GOs), ncol = length(GOs))
    
    rownames(mtx) <- GOs
    colnames(mtx) <- GOs
    
    for(i in GOs) {
      
      for(j in GOs) {
        
        mtx[i, j] <- goSim(i, j, semData = semData, ...)
        
      }
      
    }
    
    if(dist) mtx <- 1 - mtx
    
    mtx
    
  }
  
  hclust_gos <- function(GOs, semData, measure = 'Wang', method = 'ward.D2', k = 2, ...) {
    
    ## hierarchical clustering of GOs
    
    dist_mtx <- go_similarity(GOs = GOs, 
                              semData = semData, 
                              dist = TRUE, 
                              measure = measure, ...)
    
    hcl_object <- hclust(d = as.dist(dist_mtx), 
                         method = method)
    
    assignment <- cutree(hcl_object, 
                         k = k) %>% 
      compress(names_to = 'observation', 
               values_to = 'clust_id')
    
    list(data = quo(dist_mtx), 
         dist_mtx = dist_mtx, 
         dist_method = measure, 
         clust_obj = hcl_object, 
         clust_fun = 'prediction', 
         clust_assignment = assignment, 
         dots = NULL) %>% 
      clust_analysis
    
  }
  
# Metabolism ------
  
  extract_cmm_subs <- function(estimate_lst, 
                               m = 4) {
    
    ## extracts common up- and downregulated genes
    ## for a given metabolic subsystem shared by m cohorts
    
    cmm_sets <- names(estimate_lst) %>% 
      combn(m = m, 
            simplify = TRUE)
    
    reg <- list()
    
    for(i in c('upregulated', 'downregulated')) {
      
      reg[[i]] <- estimate_lst %>% 
        map(filter, regulation == i) %>% 
        map(~.x$gene_symbol)
      
      reg[[i]] <- cmm_sets %>% 
        map(~reg[[i]][.x]) %>% 
        map(reduce, intersect) %>% 
        reduce(union)
      
    }
    
    reg
    
  }
  
  plot_est_heat_map <- function(estimate_lst, 
                                gene_symbols = NULL, 
                                plot_title = NULL, 
                                plot_subtitle = NULL, 
                                plot_tag = NULL, ...) {
    
    plot_data <- estimate_lst %>% 
      compress(names_to = 'cohort')
    
    if(!is.null(gene_symbols)) {
      
      plot_data <- plot_data %>% 
        filter(gene_symbol %in% gene_symbols)
      
    }
    
    plot_data %>% 
      ggplot(aes(x = cohort, 
                 y = reorder(gene_symbol, estimate), 
                 fill = estimate)) + 
      geom_tile() + 
      scale_fill_gradient2(low = 'steelblue', 
                           mid = 'black', 
                           high = 'firebrick', 
                           midpoint = 0, ...) + 
      scale_x_discrete(labels = globals$study_labels) + 
      globals$common_theme + 
      theme(axis.title = element_blank(), 
            legend.title = element_markdown(), 
            axis.text.y = element_text(face = 'italic')) + 
      #facet_grid(. ~ level, 
       #          space = 'free', 
        #         scales = 'free') + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           fill = 'log<sub>2</sub> fold-regulation<br>vs Collagen low')
    
    
  }
  
# Markdown functions ------
  
  set_widths <- function(flextable_object,
                         widths,
                         unit = 'cm') {
    
    ## sets widths of flextable columns
    
    stopifnot(is.numeric(widths))
    
    if(length(widths) != length(flextable_object$col_keys)) {
      
      stop("Improper length of the 'widths' argument.", call. = FALSE)
      
    }
    
    for(i in seq_along(widths)) {
      
      flextable_object <- flextable_object %>%
        width(j = i, width = widths[i], unit = 'cm')
      
    }
    
    flextable_object
    
  }
  
  my_word <- function(...) {
    
    form <- word_document2(number_sections = FALSE, 
                           reference_docx = 'ms_template.docx')
    
    form$pandoc$lua_filters <- c(form$pandoc$lua_filters, 
                                 'scholarly-metadata.lua', 
                                 'author-info-blocks.lua')
    
    form
    
  }
  
# ANOVA ------
  
  get_etasq <- function(response, split_factor, data) {
    
    form <- paste(response, split_factor, sep = '~') %>% 
      as.formula
    
    aov_model <- aov(formula = form, 
                     data = data)
    
    return(rstatix::eta_squared(aov_model))
    
  }
  
# END -----