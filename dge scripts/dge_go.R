# GO enrichment analysis
# for the transcripts up- and downregulated in the Collagen Score strata

  insert_head()
  
# container ------
  
  dge_go <- list()
  
# globals -------
  
  insert_msg('Analysis globals')
  
  ## plot title prefixes
  
  dge_go$prefixes <-
    c(int = 'Collagen int vs low', 
      high = 'Collagen high vs low')
  
  ## ancestor annotation database
  
  dge_go$go_db <- godata('org.Hs.eg.db', ont = "BP", computeIC = FALSE)

# serial GO enrichment analysis --------
  
  insert_msg('GO enrichment analysis')
  
  plan('multisession')
  
  dge_go$test$int <- dge$entrez_collagen_int %>% 
    future_map(goana, 
               .options = furrr_options(seed = TRUE))
  
  dge_go$test$high <- dge$entrez_collagen_hi %>% 
    future_map(goana, 
               .options = furrr_options(seed = TRUE))
  
  plan('sequential')
  
# formatting of the results: BD GO and FDR correction -----
  
  insert_msg('Formatting of the results')

  ## appending the result tables with the total number of genes
  ## and number of regulated genes
  
  dge_go$test$int <- 
    list(x = dge_go$test$int, 
         y = dge$entrez_collagen_int, 
         z = dge$annotation) %>% 
    pmap(function(x, y, z) x %>% 
           mutate(n_de = length(y), 
                  n_total = nrow(z)))
  
  dge_go$test$high <- 
    list(x = dge_go$test$high, 
         y = dge$entrez_collagen_hi, 
         z = dge$annotation) %>% 
    pmap(function(x, y, z) x %>% 
           mutate(n_de = length(y), 
                  n_total = nrow(z)))
  
  ## effect sizes of the enrichment: odds ratio
  ## this is computed with the following formula
  ## (n(DE, GO)/n(GO))/(n(DE)/n(total))
  
  dge_go$test <- dge_go$test %>% 
    map(map, 
        mutate, 
        OR = (DE/N)/(n_de/n_total))
  
  ## multiple testing adjustment
  ## significant GOs: pFDR < 0.05
  ## effect size (Cohen 1988): at least small
  
  dge_go$test <- dge_go$test %>% 
    map(~map(.x, rownames_to_column, 'go_id') %>% 
          map(mutate, 
              p_adjusted = p.adjust(P.DE), 
              eff_size = interpret_oddsratio(OR, rules = 'cohen1988')) %>% 
          map(filter, 
              Ont == 'BP') %>% 
          map(as_tibble))
  
# GO lexicon -----
  
  insert_msg('GO lexicon')
  
  dge_go$lexicon <- dge_go$test[[1]][[1]] %>% 
    select(go_id, Term)
  
# identification of significantly enriched GOs ------
  
  insert_msg('Significant and near-significant GOs')
  
  ## significant GOs: pFDR < 0.05, at least small effect size for enrichment
  ## OR > 1 (I'm looking for the positive enrichment)
  
  dge_go$significant <- dge_go$test %>% 
    map(map, 
        filter, 
        p_adjusted < 0.05, 
        OR > 1, 
        eff_size %in% c('small', 'medium', 'large'))
  
# Plotting top enriched GOs per cohort ------
  
  insert_msg('Plots of the top enriched GOs')
  
  ## top 10 GOs
  
  dge_go$top_go$data <- dge_go$significant %>% 
    map(map, 
        filter, 
        stri_count_words(Term) < 5) %>% 
    map(map, top_n, n = 10, OR) %>% 
    map(map, 
        mutate, 
        Term = stri_replace(Term, 
                            fixed = 'extracellular matrix', 
                            replacement = 'ECM'))
  
  ## plots
  
  for(i in names(dge_go$top_go$data)) {
    
    dge_go$top_go$plots[[i]] <- 
      list(x = dge_go$top_go$data[[i]], 
           y = paste(dge_go$prefixes[[i]], 
                     globals$study_labels[names(dge_go$top_go$data[[i]])], 
                     sep = ', ')) %>% 
      pmap(function(x, y) x %>% 
             ggplot(aes(x = OR, 
                        y = reorder(Term, OR))) + 
             geom_bar(stat = 'identity', 
                      color = 'black', 
                      fill = set_names(globals$cluster_colors[c("Collagen int", 
                                                                "Collagen hi")], 
                                       c('int', 'high'))[[i]]) + 
             scale_y_discrete(labels = function(x) map(x, space2break, n = 3)) + 
             globals$common_theme + 
             theme(axis.title.y = element_blank()) + 
             labs(title = y, 
                  x = 'OR over genomic frequency'))
    
  }
  
# Common enriched GOs ------
  
  insert_msg('Common enriched GOs')
  
  ## in at least 4 out of 5 investigated cohorts
  
  dge_go$cmm_sets <- names(dge_go$test$int) %>% 
    combn(m = 4, simplify = FALSE)
  
  for(i in names(dge_go$significant)) {
    
    dge_go$common[[i]] <- dge_go$cmm_sets %>% 
      map(function(set) dge_go$significant[[i]][set] %>% 
            map(~.x$Term) %>% 
            reduce(intersect)) %>% 
      reduce(union)
    
    dge_go$common_ids[[i]] <- dge_go$cmm_sets %>% 
      map(function(set) dge_go$significant[[i]][set] %>% 
            map(~.x$go_id) %>% 
            reduce(intersect)) %>% 
      reduce(union)
    
  }
  
# similarity between the common enriched GOs -------
  
  insert_msg('Wang similatrity between common GOs')

  ## hierarchical clustering of the common regulated GOS
  
  set.seed(12345)
  
  dge_go$go_clust <- 
    list(GOs = dge_go$common_ids, 
         k = c(3, 3)) %>% 
    pmap(hclust_gos, 
         semData = dge_go$go_db)
  
  ## GO cluster assignment, annotation with human-friendly names
  
  dge_go$go_assignment <- dge_go$go_clust %>% 
    map(~.x$clust_assignment) %>% 
    map(set_names, c('go_id', 'clust_id')) %>% 
    map(mutate, 
        clust_id = paste0('clust_', clust_id), 
        clust_id = factor(clust_id)) %>% 
    map(left_join, 
        dge_go$lexicon, 
        by = 'go_id')
  
  ## setting descriptive names of the clusters
  
  dge_go$go_clust_desc$int <- 
    c(clust_1 = 'ECM\ncollagen', 
      clust_2 = 'tissue\ndevelopment\nangiogenesis', 
      clust_3 = 'adhesion\nsignaling\ngrowth factors')
  
  dge_go$go_clust_desc$high <- 
    c(clust_1 = 'ECM\ncollagen', 
      clust_2 = 'tissue\ndevelopment\nangiogenesis', 
      clust_3 = 'motility')

# Plotting the common significant GOs -----
  
  insert_msg('Heat bubble plots with the significances of common enriched GOs')
  
  ## plotting tables
  
  dge_go$cmm_plot_tbl <- map2(dge_go$test, 
                              dge_go$common, 
                              function(data, term) data %>% 
                                map(filter, Term %in% term) %>%
                                compress(names_to = 'cohort')) %>% 
    map2(dge_go$go_assignment, 
         left_join, by = c('go_id', 'Term')) %>% 
    map(mutate,
        cohort = factor(cohort, names(dge_go$test[[1]])), 
        significant = ifelse(p_adjusted < 0.05 & 
                               OR > 1 & 
                               eff_size %in% c('small', 'medium', 'large'), 
                             'yes', 'no'), 
        fontface = ifelse(significant == 'yes', 'bold', 'plain'))
    
  
  ## heat bubble plots
  
  dge_go$cmm_plots <- 
    map2(dge_go$cmm_plot_tbl, 
         paste('Common GO terms', 
               dge_go$prefixes, 
               sep = ', '), 
         ~ggplot(.x, 
                 aes(x = cohort, 
                     y = reorder(go_id, OR), 
                     fill = OR, 
                     size = OR)) + 
           geom_point(shape = 21) + 
           geom_text(aes(label = signif(OR, 2), 
                         alpha = significant, 
                         fontface = fontface, 
                         x = as.numeric(cohort) + 0.25), 
                     size = 2.5, 
                     hjust = 0, 
                     vjust = 0.5) + 
           scale_x_discrete(labels = globals$study_labels) + 
           scale_fill_gradient(low = 'white', 
                               high = 'firebrick') + 
           scale_radius(range = c(1, 4)) + 
           scale_alpha_manual(values = c(no = 0.5, 
                                         yes = 1)) + 
           guides(fill = 'legend', 
                  size = 'legend', 
                  alpha = 'none') + 
           globals$common_theme + 
           theme(axis.title = element_blank(), 
                 axis.text.y = element_text(size = 8), 
                 strip.text.y = element_text(size = 8, 
                                             angle = 0, 
                                             hjust = 0)) + 
           labs(title = .y, 
                fill = 'enrichment\nOR', 
                size = 'enrichment\nOR'))
  
  ## faceting by the GO cluster assignment
  
  dge_go$cmm_plots <- 
    map2(dge_go$cmm_plots, 
         dge_go$go_clust_desc, 
         ~.x + 
           facet_grid(clust_id ~ ., 
                      space = 'free', 
                      scales = 'free', 
                      labeller = as_labeller(.y)))
  
# END -----
  
  rm(i)
  
  insert_tail()