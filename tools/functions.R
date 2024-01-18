# This script provides project specific tools ----

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
  library(ggrepel)

# data import ------

  annotate_raw_symbol <- function(data) {
    
    data %>%
      filter(!is.na(entrez_id),
             entrez_id != '') %>%
      mutate(gene_symbol = mapIds(org.Hs.eg.db,
                                  column = 'SYMBOL',
                                  keys = entrez_id,
                                  keytype = 'ENTREZID')) %>%
      mutate(gene_symbol = unlist(gene_symbol)) %>%
      filter(gene_symbol != '') %>%
      filter(complete.cases(.))
    
  }

  integrate_expression <- function(data, annotation, fun = colMeans, ...) {
    
    ## integrates expression over multiple probes via arithmetic mean
    
    dupl_genes <- annotation %>%
      filter(duplicated(gene_symbol)) %>%
      .$gene_symbol %>%
      unique
    
    single_annot <- annotation %>%
      filter(!gene_symbol %in% dupl_genes)
    
    dupl_annot <- annotation %>%
      filter(gene_symbol %in% dupl_genes)
    
    if(nrow(dupl_annot) == 0) {
      
      data <- data[annotation$probe_id, ] %>%
        set_rownames(annotation$gene_symbol)
      
    } else {
      
      agg_data <- aggRows(x = data[dupl_annot$probe_id, ],
                          f = dupl_annot$gene_symbol,
                          fun = fun,
                          na.rm = TRUE, ...)
      
      if(nrow(single_annot) == 0) {
        
        data <- agg_data
        
      } else {
        
        single_data <- data[single_annot$probe_id, ] %>%
          set_rownames(single_annot$gene_symbol)
        
        data <- rbind(single_data, agg_data)
        
      }
      
    }
    
    data <- data %>%
      t %>%
      as.data.frame %>%
      rownames_to_column('sample_id')
    
    as_tibble(data)
    
  }

  extract_mutations <- function(data,
                                id_variable = 'Tumor_Sample_Barcode',
                                gene_id = 'Hugo_Symbol') {
    
    ## counts mutations per gene and sample
    
    data[[gene_id]] <- factor(data[[gene_id]])
    
    data_splits <- split(data[[gene_id]], data[[id_variable]])
    
    mut_counts <- data_splits %>%
      map(table)
    
    sample_ids <- names(mut_counts)
    
    mut_counts <- do.call('rbind', mut_counts)
    
    mut_counts <- ifelse(mut_counts > 0, 1, mut_counts)
    
    mut_counts <- as_tibble(mut_counts) %>%
      future_map_dfc(~car::recode(as.character(.x),
                                  "'0' = 'WT'; '1' = 'mutated'"),
                     .options = furrr_options(seed = TRUE)) %>%
      future_map_dfc(factor, c('WT', 'mutated'),
                     .options = furrr_options(seed = TRUE))
    
    mut_counts %>%
      mutate(sample_id = sample_ids) %>%
      relocate(sample_id)
    
  }

  box_from_stats <- function(data,
                             width = 0.8,
                             plot_title = NULL,
                             plot_subtitle = NULL,
                             y_lab = 'log<sub>2</sub> expression',
                             labeller = c('expression' = 'raw',
                                          'corr_expression' = 'COMBAT-corrected'),
                             scales = 'free') {
    
    ## plots a box plot given a data frame with the columns
    ## `median`, `perc25`, `perc75`, `perc025` and `perc975`
    ## storing the median expression, interquartile and 95 percentile ranges
    
    data <- data %>%
      mutate(center_pos = as.numeric(cohort))
    
    data %>%
      ggplot(aes(x = cohort,
                 y = median,
                 fill = cohort)) +
      facet_grid(type ~ .,
                 labeller = as_labeller(labeller),
                 scales = scales) +
      geom_errorbar(aes(ymin = perc025,
                        ymax = perc975),
                    width = 0) +
      geom_rect(aes(xmin = center_pos - 0.5 * width,
                    xmax = center_pos + 0.5 * width,
                    ymin = perc25,
                    ymax = perc75),
                color = 'black') +
      geom_segment(aes(x = center_pos - 0.5 * width,
                       xend = center_pos + 0.5 * width,
                       y = median,
                       yend = median),
                   color = 'black') +
      scale_fill_manual(values = globals$cohort_colors,
                        labels = globals$cohort_labs,
                        name = '') +
      scale_x_discrete(labels = globals$cohort_labs) +
      guides(fill = 'none') +
      globals$common_theme +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_markdown(),
            axis.text.x = element_text(hjust = 1,
                                       vjust = 0.5,
                                       angle = 90)) +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           y = y_lab)
    
  }

# data transformation -------

  safely_mutate <- function(.data, ...) {
    
    ## enables for error-resistant change of data frame columns
    ## in case of any error, the genuine data frame is returned
    
    new_data <- try(mutate(.data, ...),
                    silent = TRUE)
    
    if(inherits(new_data, 'try-error')) return(.data)
    
    return(new_data)
    
  }
  
  safely_select <- function(.data, cols, fill = NULL) {
    
    ## selects columns with a character vector of column names
    ## in safe mode, i.e. without errors if the column does not exist
    ##
    ## if `fill` is specified, the non-existing columns are added and filled
    ## with the defined value
    
    stopifnot(is.character(cols))
    
    new_data <- select(.data, any_of(cols))
    
    if(is.null(fill)) return(new_data)
    
    absent_cols <- cols[!cols %in% names(.data)]
    
    for(i in absent_cols) {
      
      new_data[[i]] <- fill
      
    }
    
    new_data
    
  }
  
  minimum_shift <- function(df) {
    
    ## shifts values of variables of a numeric data frame by their minima
    
    min_values <- map_dbl(df, min, ns.rm = TRUE)
    
    map2_dfc(df, abs(min_values), `+`)
    
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
  
# Cluster characteristic ---------
  
  plot_stage_stack <- function(data_lst, 
                               test_lst = NULL, 
                               variable = 'gleason_simple', 
                               palette = c('cornsilk', 'cornsilk3', 
                                           'coral2', 'coral4'), 
                               plot_title = NULL, 
                               plot_subtitle = NULL, 
                               fill_lab = variable, 
                               x_lab = '% of cluster', 
                               show_clust_n = FALSE) {
    
    ## plotting, setting the common levels
    
    plot_data <- data_lst %>% 
      map(safely(select), 
          sample_id, 
          clust_id, 
          all_of(variable)) %>% 
      map(~.x$result) %>% 
      compact %>% 
      map(~filter(.x, complete.cases(.x)))
    
    plot_data <- plot_data %>% 
      map(blast, clust_id) %>% 
      map(map, count, .data[[variable]]) %>% 
      map(map, arrange, desc(.data[[variable]])) %>% 
      map(map, 
          mutate, 
          n_total = sum(n), 
          percent = n/n_total  * 100, 
          x_pos = cumsum(percent))
    
    plot_data <- plot_data %>% 
      map(compress, names_to = 'clust_id') %>% 
      compress(names_to = 'cohort') %>% 
      mutate(clust_id = factor(clust_id, levels(data_lst[[1]]$clust_id)), 
             cohort = factor(cohort, rev(names(data_lst))))
    
    plot_data <- plot_data %>% 
      mutate(!!variable := factor(.data[[variable]], 
                                  rev(levels(plot_data[[variable]]))))
    
    plot_data <- plot_data %>% 
      mutate(ax_lab = stri_extract(clust_id, regex = 'low|hi'))
    
    if(show_clust_n) {
      
      plot_data <- plot_data %>% 
        mutate(ax_lab = paste(ax_lab, n_total, sep = ': n = '))
      
    } 
    
    ## labeller ---------
    
    if(is.null(test_lst)) {
      
      labeller <- globals$study_labels
      
    } else {
      
      test_lst <- test_lst %>% 
        map(filter, variable == .env[['variable']]) %>% 
        compress(names_to = 'cohort') %>% 
        mutate(plot_cap = paste(eff_size, significance, sep = '\n'), 
               cohort_lab = globals$study_labels[cohort], 
               label = paste(cohort_lab, plot_cap, sep = '\n'))
      
      labeller <- set_names(test_lst$label, test_lst$cohort)
      
    }
    
    ## plotting --------
    
    plot_data %>% 
      ggplot(aes(x = percent, 
                 y = ax_lab, 
                 fill = .data[[variable]]))  +
      geom_bar(stat = 'identity', 
               color = 'black') + 
      facet_grid(cohort ~ ., 
                 labeller = as_labeller(labeller), 
                 space = 'free', 
                 scales = 'free') + 
      scale_fill_manual(values = rev(palette), 
                        name = fill_lab) + 
      globals$common_theme + 
      theme(axis.title.y = element_blank()) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           x = x_lab)
    
  }
  
  plot_clinic_box <- function(data_lst, 
                              test_lst = NULL, 
                              variable = 'psa_diagnosis', 
                              normalize = TRUE, 
                              norm_center = 'mean', 
                              palette = globals$cluster_colors, 
                              plot_title = variable, 
                              plot_subtitle = NULL, 
                              x_lab = 'PSA at diagnosis', 
                              fill_lab = '', 
                              show_clust_n = FALSE, 
                              point_size = 2, 
                              point_alpha = 0.75, 
                              point_hjitter = 0.1, 
                              point_wjitter = 0, 
                              shape_alpha = 0.25, 
                              dodge_w = 0.9) {
    
    ## plotting data -------
    
    n_numbers <- data_lst %>% 
      map(count, clust_id) %>% 
      map(mutate, n_total = sum(n)) %>% 
      compress(names_to = 'cohort')
    
    plot_data <- data_lst %>% 
      map(safely(select), sample_id, clust_id, all_of(variable)) %>% 
      map(~.x$result) %>% 
      compact
    
    if(normalize) {
      
      norm_fun <- switch(norm_center, 
                         mean = function(x) mean(x, na.rm = TRUE), 
                         median = function(x) median(x, na.rm = TRUE))
      
      plot_data <- plot_data %>% 
        map(mutate, 
            !!variable := scale(.data[[variable]], norm_fun(.data[[variable]])))
      
    }
    
    plot_data <- plot_data %>% 
      compress(names_to = 'cohort') %>% 
      left_join(n_numbers, by = c('cohort', 'clust_id')) %>% 
      mutate(ax_lab = stri_extract(clust_id, regex = 'low|hi'), 
             cohort = factor(cohort, rev(names(data_lst))))
    
    if(show_clust_n) {
      
      plot_data <- plot_data %>% 
        mutate(ax_lab = paste(ax_lab, n_total, sep = ': n = '))
      
    } 
    
    ## labeller ---------
    
    if(is.null(test_lst)) {
      
      labeller <- globals$study_labels
      
    } else {
      
      test_lst <- test_lst %>% 
        map(filter, variable == .env[['variable']]) %>% 
        compress(names_to = 'cohort') %>% 
        mutate(plot_cap = paste(eff_size, significance, sep = '\n'), 
               cohort_lab = globals$study_labels[cohort], 
               label = paste(cohort_lab, plot_cap, sep = '\n'))
      
      labeller <- set_names(test_lst$label, test_lst$cohort)
      
    }
    
    plot_data %>% 
      ggplot(aes(x = .data[[variable]], 
                 y = ax_lab, 
                 fill = clust_id, 
                 color = clust_id))  +
      geom_boxplot(outlier.color = NA, 
                   color = 'black', 
                   alpha = shape_alpha, 
                   position = position_dodge(dodge_w)) +
      geom_point(shape = 16, 
                 size = point_size, 
                 position = position_jitterdodge(dodge.width = dodge_w, 
                                                 jitter.width = point_hjitter, 
                                                 jitter.height = point_wjitter)) + 
      facet_grid(cohort ~ ., 
                 labeller = as_labeller(labeller), 
                 space = 'free', 
                 scales = 'free') + 
      scale_fill_manual(values = palette, 
                        name = fill_lab) + 
      scale_color_manual(values = palette, 
                         name = fill_lab) + 
      globals$common_theme + 
      theme(axis.title.y = element_blank()) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           x = x_lab)
    
  }
  
  extract_n_numbers <- function(data, 
                                split_factor = 'clust_id', 
                                variable_txt = 'Samples, n') {
    
    n_numbers <- data[[split_factor]] %>% 
      table
    
    n_numbers %>% 
      reduce(cbind) %>% 
      set_colnames(names(n_numbers)) %>% 
      as_tibble %>% 
      set_names(names(n_numbers)) %>% 
      mutate(variable = 'Samples, n') %>% 
      relocate(variable)
    
  }
  
  format_desc <- function(desc_data) {
    
    levs <- names(desc_data)
    
    desc_data %>% 
      reduce(left_join, by = 'variable') %>% 
      set_names(c('variable', levs))
    
  }
  
  format_two_test <- function(test_res) {
    
    test_res %>% 
      re_adjust(p_variable = 'p_value', 
                method = 'none') %>% 
      mutate(variable = response, 
             eff_lab = paste(effect_size_name, 
                             signif(effect_size, 2), 
                             sep = ' = '), 
             plot_cap = paste(eff_lab, significance, sep = ', '))
    
  }
  
  plot_box_panel <- function(data_lst, 
                             test_lst, 
                             variable, 
                             split_factor = 'clust_id', 
                             normalize = FALSE, 
                             norm_center = 'mean', 
                             p_cutoff = 0.05, 
                             effect_cutoff = 0.1, 
                             show_stats = TRUE, 
                             point_size = 2, 
                             point_alpha = 0.75, 
                             point_hjitter = 0.1, 
                             point_wjitter = 0, 
                             shape_alpha = 0.25, 
                             dodge_w = 0.9, 
                             palette = globals$cluster_colors, 
                             plot_title = variable, 
                             plot_subtitle = NULL, 
                             x_lab = 'variable level', 
                             fill_lab = NULL) {
    
    ## draws a panel of box plots for a given variable:
    ## cohorts are shown in the Y axis, variable levels in the X axis.
    ## significant effects are labeled in bold (defined by the p and effect 
    ## size cutoffs)
    
    ## axis labels --------
    
    plot_stats <- test_lst %>% 
      map(filter, response == .env[['variable']]) %>% 
      compress(names_to = 'cohort') %>% 
      format_two_test %>% 
      mutate(cohort_lab = globals$study_labels[cohort], 
             significant = ifelse(p_adjusted < p_cutoff & abs(effect_size >= effect_cutoff), 
                                  'yes', 'no'), 
             plot_cap = paste(cohort_lab, plot_cap, sep = '<br>'), 
             plot_cap = ifelse(significant == 'yes', 
                               html_bold(plot_cap), plot_cap), 
             cohort_lab = ifelse(significant == 'yes', 
                                 html_bold(cohort_lab), cohort_lab))
    
    
    if(show_stats) {
      
      ax_labs <- set_names(plot_stats$plot_cap, plot_stats$cohort)
      
    } else {
      
      ax_labs <- set_names(plot_stats$cohort_lab, plot_stats$cohort)
      
    }
    
    ## plotting data -------
    
    plot_data <- data_lst %>% 
      map(select, all_of(c(split_factor, variable)))
    
    if(normalize) {
      
      fun <- switch(norm_center, 
                    mean = function(x) mean(x, na.rm = TRUE), 
                    median = function(x) median(x, na.rm = TRUE))
      
      plot_data <- plot_data %>% 
        map(mutate, 
            !!variable := scale(.data[[variable]], 
                                center = fun(.data[[variable]]))[, 1])
      
    }
    
    plot_data <- plot_data %>% 
      compress(names_to = 'cohort') %>% 
      mutate(cohort = factor(cohort, names(data_lst)))
    
    ## box plots -------
    
    plot_data %>% 
      ggplot(aes(x = .data[[variable]], 
                 y = cohort, 
                 fill = .data[[split_factor]], 
                 color = .data[[split_factor]])) +
      geom_boxplot(alpha = shape_alpha, 
                   color = 'black', 
                   outlier.color = NA, 
                   position = position_dodge(dodge_w)) + 
      geom_point(shape = 16, 
                 size = point_size, 
                 alpha = point_alpha, 
                 position = position_jitterdodge(jitter.width = point_hjitter, 
                                                 jitter.height = point_wjitter, 
                                                 dodge.width = dodge_w)) + 
      scale_fill_manual(values = palette, 
                        name = fill_lab) + 
      scale_color_manual(values = palette, 
                         name = fill_lab) + 
      scale_y_discrete(labels = ax_labs) + 
      globals$common_theme + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_markdown()) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           x = x_lab)
    
  }
  
  cohort_heat_map <- function(data_lst, 
                              variables, 
                              split_fct = 'clust_id', 
                              normalize = FALSE, ...) {
    
    ## average feature levels in the cohorts
    
    ## plotting data: cohort averages
    
    mean_data <-  data_lst %>% 
      map(select, 
          all_of(c(split_fct, variables))) %>% 
      map(blast, all_of(split_fct), .skip = TRUE) %>% 
      map(map, map_dfc, mean) %>% 
      map(compress, names_to = split_fct) %>% 
      compress(names_to = 'cohort') %>% 
      mutate(!!split_fct:= factor(.data[[split_fct]], 
                                  levels(data_lst[[1]][[split_fct]])))
    
    cohort_lexicon <- 
      tibble(observation = paste0('obs_', 1:nrow(mean_data)), 
             cohort = mean_data$cohort)
    
    ## plotting 
    
    hm <-  heat_map(data = mean_data, 
                    variables = variables, 
                    split_fct = split_fct, 
                    normalize = normalize, 
                    cust_theme = globals$common_theme, ...)
    
    ## appending with the cohort names and their order
    
    hm$data <- hm$data %>% 
      mutate(observation = exchange(observation, 
                                    dict = cohort_lexicon, 
                                    key = 'observation', 
                                    value = 'cohort'),
             observation = factor(observation, names(data_lst)))
    
    hm + 
      scale_x_discrete(labels = globals$study_labels)
    
  }
  
  regulation_bubble <- function(data_lst, 
                                variables, 
                                label_variable = 'Name', 
                                regulation_variable = 'tA', 
                                status_variable = 'regulation', 
                                plot_title = NULL, 
                                plot_subtitle = NULL, 
                                fill_lab = NULL, 
                                size_lab = 'Pathway regulation\ntA', 
                                palette = c(activated = 'firebrick', 
                                            inhibited = 'steelblue', 
                                            ns = 'gray60'), 
                                abs_estimate = FALSE, 
                                show_estimates = TRUE, ...) {
    
    ## draws a bubble plot with point size coding for absolute regulation
    ## magnitude and point color coding for the regulation sign
    
    ## plotting data --------
    
    plot_tbl <- data_lst %>% 
      compress(names_to = 'cohort') %>% 
      filter(.data[[label_variable]] %in% variables) %>% 
      mutate(cohort = factor(cohort, names(data_lst)), 
             significant = ifelse(.data[[status_variable]] == 'ns', 
                                  'no', 'yes'), 
             fontface = ifelse(.data[[status_variable]] != 'ns', 
                               'bold', 'plain'))
    
    ## plotting -------
    
    bubble <- plot_tbl %>% 
      ggplot(aes(x = cohort, 
                 y = reorder(.data[[label_variable]], 
                             .data[[regulation_variable]]), 
                 fill = .data[[status_variable]], 
                 color = .data[[status_variable]], 
                 size = abs(.data[[regulation_variable]]))) + 
      geom_point(shape = 21, 
                 color = 'black')
    
    if(show_estimates) {
      
      if(abs_estimate) {
        
        fun <- abs
        
      } else {
        
        fun <- identity
        
      }
      
      bubble <- bubble + 
        geom_text(aes(label = signif(fun(.data[[regulation_variable]]), 2),
                      alpha = significant, 
                      x = as.numeric(cohort) + 0.25), 
                  size = 2.4, 
                  hjust = 0, 
                  vjust = 0.5, 
                  show.legend = FALSE)
      
    }
    
    bubble + 
      scale_color_manual(values = palette, 
                         drop = FALSE, 
                         name = fill_lab) + 
      scale_fill_manual(values = palette, 
                        drop = FALSE, 
                        name = fill_lab) + 
      scale_size_area(name = size_lab, ...) + 
      scale_x_discrete(labels = globals$study_labels) + 
      scale_alpha_manual(values = c(no = 0.5, 
                                    yes = 1)) + 
      guides(size = 'legend', 
             fill = 'legend', 
             alpha = 'none') + 
      globals$common_theme + 
      theme(axis.title = element_blank()) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle)
    
  }
  
  
# Text functions -----
  
  html_italic <- function(x) paste0('<em>', x, '</em>')
  
  html_bold <- function(x) paste0('<b>', x, '</b>')

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
  
  clust_labeller <- function(x) {
    
    x <- x %>% 
      stri_replace(fixed = 'squared euclidean', 
                   replacement = 'Euclidean\u00B2') %>% 
      stri_replace(fixed = 'euclidean', 
                   replacement = 'Euclidean') %>% 
      stri_replace(fixed = 'manhattan',
                   replacement = 'Manhattan')
    
    ifelse(x == 'PAM, cosine', html_bold('PAM, cosine'), x)
    
  }
  
  progeny_labeller <- function(x) {
    
    x %>% 
      stri_replace(fixed = 'TGFb', 
                   replacement = 'TGF-\u03B2') %>% 
      stri_replace(fixed = 'TNFa', 
                   replacement = 'TNF-\u03B1')
    
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

# Survival analysis --------
  
  cut_tertiles <- function(x) {
    
    terts <- quantile(x, c(1/3, 2/3), na.rm = TRUE)
    
    cuts <- c(-Inf, terts, Inf)
    
    cut(x, cuts, c('low', 'int', 'high'))
    
  }
  
  plot_reg_cox_tuning <- function(data, 
                                  best_tune, 
                                  plot_title = NULL, 
                                  x_lab = expression(lambda), 
                                  y_lab = 'model deviance, CV') {
    
    ## plots for the process of lambda tuning
    
    data <- data %>% 
      mutate(best = ifelse(cvm == min(cvm), 
                           'yes', 'no'))
    
    data %>% 
      ggplot(aes(x = lambda, 
                 y = cvm, 
                 fill = best)) + 
      geom_point(shape = 21, 
                 size = 2) + 
      scale_fill_manual(values = c(no = 'steelblue', 
                                   yes = 'coral3'), 
                        name = 'best tune') + 
      globals$common_theme + 
      labs(title = plot_title, 
           subtitle = paste('Best tune: \u03BB =', 
                            signif(best_tune$lambda[1], 3)), 
           x = x_lab,
           y = y_lab)
    
  }
  
  plot_surv_stats <- function(stats, 
                              plot_title = 'Transcriptional Collagen Score, RFS prediction', 
                              plot_subtitle = 'Ridge Cox algorithm', 
                              x_lab = 'C-index', 
                              y_lab = '1 - IBS', 
                              palette = surv_globals$model_colors, 
                              labels = surv_globals$study_labels, 
                              label_variable = 'cohort', 
                              txt_size = 2.75, 
                              color_variable = 'dataset') {
    
    ## a scatter plot of 1 - IBS and C-index
    ## cohort type is color-coded
    
    stats %>% 
      mutate(!!label_variable := labels[.data[[label_variable]]]) %>% 
      ggplot(aes(x = c_index, 
                 y = 1 - ibs_model, 
                 fill = .data[[color_variable]], 
                 color = .data[[color_variable]])) + 
      geom_vline(xintercept = 0.5, 
                 linetype = 'dashed') +
      geom_hline(yintercept = 0.75, 
                 linetype = 'dashed') + 
      geom_point(size = 2, 
                 shape = 21, 
                 color = 'black') + 
      geom_text_repel(aes(label = .data[[label_variable]]), 
                      size = txt_size, 
                      show.legend = FALSE) + 
      scale_fill_manual(values = palette, 
                        name = '') + 
      scale_color_manual(values = palette, 
                         name = '') + 
      globals$common_theme + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           x = x_lab, 
           y = y_lab)
    
  }
  
  plot_tertile_km <- function(fit, 
                              n_numbers, 
                              p_value, 
                              palette = surv_globals$tertile_colors, 
                              plot_title = NULL, 
                              x_lab = 'Relapse-free survival, months',
                              legend_title = 'Score tertile', 
                              txt_size = 2.75) {
    
    ## Kaplan-Meier plots for tertiles of a transcriptional score
    
    plot_subtitle <- paste0('total: n =', sum(n_numbers$n_total), 
                            ', events: n = ', sum(n_numbers$n_events))
    
    scale_labs <- n_numbers %>% 
      mutate(scale_lab = paste0(score_cuts, '\ntotal: n = ', 
                                n_total, '\nevents: n = ', 
                                n_events, '\n'))
    
    scale_labs <- set_names(scale_labs$scale_lab, 
                            as.character(scale_labs$score_cuts))
    
    
    km_plot <- ggsurvplot(fit = fit, 
                          pval = p_value, 
                          pval.size = txt_size, 
                          pval.coord = c(0.01, 0.05))
    
    km_plot$plot + 
      scale_color_manual(values = unname(palette),
                         labels = unname(scale_labs), 
                         name = legend_title) +
      globals$common_theme + 
      theme(legend.position = 'bottom') + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           x = x_lab)
    
  }
  
  plot_surv_importance <- function(data, 
                                   imp_stat = 'coef', 
                                   n_top = NULL, 
                                   form = c('bar', 'bubble'), 
                                   labeller = c(positive = 'unfavorable', 
                                                negative = 'favorable', 
                                                ns = 'ns'), 
                                   plot_title = 'Transcriptional Collagen Score, variable importance', 
                                   plot_subtitle = 'Ridge Cox regression', 
                                   size_title = 'abs(log HR)', 
                                   max_size = 5,
                                   x_lab = expression('log HR'[Ridge]), 
                                   line_color = 'black', 
                                   palette = c(positive = 'firebrick', 
                                               negative = 'steelblue', 
                                               ns = 'gray70')) {
    
    ## a bar or bubble plot of variable importance stats
    
    form <- match.arg(form, c('bar', 'bubble'))
    
    data <- data %>% 
      mutate(regulation = ifelse(.data[[imp_stat]] > 0, 'positive', 
                                 ifelse(.data[[imp_stat]]< 0, 'negative', 'ns')),
             regulation = factor(regulation, c('positive', 'negative', 'ns')), 
             variable = stri_replace(variable, 
                                     fixed = '_sq', 
                                     replacement = '\u00B2'))
    
    if(!is.null(n_top)) {
      
      data <- data %>% 
        blast(regulation) %>% 
        map_dfr(top_n, n = n_top, abs(.data[[imp_stat]]))
      
    }
    
    if(form == 'bar') {
      
      imp_plot <- data %>% 
        ggplot(aes(x = .data[[imp_stat]], 
                   y = reorder(variable, .data[[imp_stat]]), 
                   fill = regulation)) + 
        geom_bar(stat = 'identity', 
                 color = line_color) 
      
    } else {
      
      imp_plot <- data %>% 
        ggplot(aes(x = .data[[imp_stat]], 
                   y = reorder(variable, .data[[imp_stat]]), 
                   fill = regulation, 
                   size = abs(.data[[imp_stat]]))) + 
        geom_point(shape = 21) + 
        scale_size_area(max_size = max_size, 
                        name = size_title)
      
    }
    
    imp_plot + 
      facet_grid(regulation ~ ., 
                 scales = 'free', 
                 space = 'free', 
                 labeller = as_labeller(labeller)) + 
      scale_fill_manual(values = palette) + 
      globals$common_theme + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_text(face = 'italic')) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle,
           x = x_lab)
    
  }
  
  hr_from_cut <- function(survcut_obj) {
    
    ## fits a Cox PH model to expression data stratified by the optimal cutoff
    ## used to identify favorable and unfavorable markers.
    ## For this reason we're not carrying much about assumptions
    
    cox_data <- survcut_obj %>% 
      model.frame %>% 
      select(rfs_months, relapse, ends_with('_strata'))
    
    cox_model <- coxph(Surv(rfs_months, relapse) ~ ., 
                       data = cox_data)
    
    log_hr <- coef(cox_model)
    
    log_hr_ci <- confint(cox_model)
    
    tibble(beta = log_hr, 
           beta_lower_ci = log_hr_ci[1], 
           beta_upper_ci = log_hr_ci[2], 
           hr = exp(log_hr), 
           hr_lower_ci = exp(log_hr_ci)[1], 
           hr_upper_ci = exp(log_hr_ci)[2]) %>% 
      mutate(significant = ifelse(all(log_hr_ci > 0) | all(log_hr_ci < 0), 
                                  'yes', 'no'), 
             marker = ifelse(significant == 'no', 'ns', 
                             ifelse(beta > 0, 'unfavorable', 
                                    ifelse(beta < 0, 'favorable', 'ns'))), 
             marker = factor(marker, c('unfavorable', 'favorable', 'ns'))) %>% 
      select(-significant)
    
  } 
  
  
# gene, GO and signature similarity -----
  
  clust_gos <- function(GOs, 
                        semData, 
                        measure = 'Wang', 
                        mds_dim = 2, 
                        fun = kcluster, ...) {
    
    ## clustering of GOs by semantic (functional) similarity
    
    ## distance matrix
    
    dist_mtx <- go_sem(GOs, 
                       semData = semData, 
                       as_dist = TRUE, 
                       measure = measure)
    
    ## MDS: a dimensionality reduction step is necessary for large GO vectors
    
    mds_scores <- cmdscale(dist_mtx, k = mds_dim) %>% 
      as.data.frame %>% 
      set_names(c('comp_1', 'comp_2')) %>% 
      rownames_to_column('observation')
    
    model_frame <- enexpr(mds_scores)
    
    mds_obj <- list(red_obj = NULL, 
                    red_fun = 'mds', 
                    dist_method = 'custom', 
                    component_tbl = mds_scores, 
                    loadings = NULL, 
                    data = quo(model_frame)) %>% 
      red_analysis
    
    fun(mds_obj, ...)
    
  }
  
  clust_signature <- function(data, 
                              variables, 
                              mds_distance = 'cosine', 
                              mds_dim = 2, 
                              fun = hcluster, ...) {
    
    ## clustering by similarity 
    
    clust_data <- data %>% 
      select(sample_id, all_of(variables)) %>% 
      column_to_rownames('sample_id') %>% 
      select(all_of(variables)) %>% 
      center_data('mean') %>% 
      t %>% 
      reduce_data(distance_method = mds_distance, 
                  kdim = mds_dim, 
                  red_fun = 'mds')
    
    fun(clust_data, ...)

  }
  
# Metabolism ------

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

  my_word <- function(...) {
    
    form <- word_document2(number_sections = FALSE, 
                           reference_docx = 'ms_template.docx')
    
    form$pandoc$lua_filters <- c(form$pandoc$lua_filters, 
                                 'scholarly-metadata.lua', 
                                 'author-info-blocks.lua')
    
    form
    
  }
  
  format_strata_n <- function(count_data) {
    
    map2_chr(count_data[[1]], count_data[[2]], 
             paste, sep = ': n = ') %>% 
      paste(collapse = ', ')
    
  }
  
# END -----