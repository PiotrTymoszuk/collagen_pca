# Identification of differentially regulated genes 
# between the Collagen clusters
# Differentially regulated genes are identified by one-way ANOVA
# corrected for multiple testing with the Benjamini-Hochberg method
# Regulation cutoff: 1.25-fold

  insert_head()
  
# container ------
  
  dge <- list()
  
# Gene variables -----
  
  insert_msg('Gene variables')

  ## gene variables per study
  ## only probes with annotated Entrez ID are included in the analysis

  dge$variables <- study_data %>% 
    map(~.x$annotation$symbol)
  
  dge$variables <- 
    map2(dge$variables, 
         map(study_data, ~.x$expression), 
         ~.x[.x %in% names(.y)]) %>% 
    map(~.x[!duplicated(.x)])
  
  ## annotation lexicons
  
  dge$annotation_symbol <- dge$variables %>% 
    map(~mapIds(org.Hs.eg.db, 
                keys = .x, 
                column = 'ENTREZID', 
                keytype = 'SYMBOL')) %>% 
    map(~tibble(gene_symbol = names(.x), 
                entrez_id = .x))
  
  dge$annotation_alias <- dge$variables %>% 
    map(~mapIds(org.Hs.eg.db, 
                keys = .x, 
                column = 'ENTREZID', 
                keytype = 'ALIAS')) %>% 
    map(~tibble(gene_symbol = names(.x), 
                entrez_id = .x))
  
  dge$annotation <- map2(dge$annotation_symbol, 
                         dge$annotation_symbol, 
                         rbind) %>% 
    map(filter, 
        !is.na(gene_symbol), 
        !is.na(entrez_id)) %>% 
    map(filter, 
        !duplicated(gene_symbol))

  dge$annotation_symbol <- NULL
  dge$annotation_alias <- NULL
  
  dge <- compact(dge)
  
  ## restricting the variable list to the genes with known Entrez IDs
  
  dge$variables <- map2(dge$variables, 
                        dge$annotation, 
                        ~.x[.x %in% .y$gene_symbol])
  
# Analysis tables -------
  
  insert_msg('Analysis tables')

  ## analysis tables
  
  dge$assignment <- coll_clust$assignment
  
  dge$analysis_tbl <- study_data %>% 
    map(~.x$expression) %>% 
    map(extract_tumor_samples) %>% 
    map2(., dge$variables, 
         ~.x[c('patient_id', .y)]) %>% 
    map2(dge$assignment, ., 
         left_join, by = 'patient_id') %>% 
    map(mutate, 
        clust_id = factor(clust_id, 
                          c('Collagen low', 
                            'Collagen int', 
                            'Collagen hi')))
  
# Serial ANOVA -------
  
  insert_msg('Serial testing')
  
  dge$test_results <- list(data = dge$analysis_tbl, 
                           variables = dge$variables) %>% 
    pmap(test_anova, 
         split_fct = 'clust_id', 
         type = 't', 
         adj_method = 'BH', 
         .parallel = TRUE)
  
# Annotation of the results with Entrez IDs -------
  
  insert_msg('Annotation of the results')
  
  for(i in names(dge$analysis_tbl)) {
    
    dge$test_results[[i]]$anova <- dge$test_results[[i]]$anova %>% 
      mutate(gene_symbol = response, 
             entrez_id = exchange(gene_symbol, 
                                  dge$annotation[[i]], 
                                  key = 'gene_symbol', 
                                  value = 'entrez_id'))
    
    dge$test_results[[i]]$lm <- dge$test_results[[i]]$lm %>% 
      mutate(gene_symbol = response, 
             entrez_id = exchange(gene_symbol, 
                                  dge$annotation[[i]], 
                                  key = 'gene_symbol', 
                                  value = 'entrez_id'))
    
    
  }

# Formatting the results and identifying differentially regulated genes ------
  
  insert_msg('Identification of differentially regulated genes')
  
  dge$anova_signif_results <- dge$test_results %>% 
    map(~.$anova) %>% 
    map(filter, significant == 'yes')
  
  dge$lm_signif_results <- 
    map2(map(dge$test_results, ~.x$lm), 
         map(dge$anova_signif_results, ~.x$response), 
         ~filter(.x, response %in% .y))
  
  ## indicating gene regulation sign
  
  dge$lm_signif_results <- dge$lm_signif_results %>% 
    map(mutate, 
        level = ifelse(level == '(Intercept)', 
                       'Collagen low', 
                       stri_extract(level, regex = 'Collagen.*$')), 
        level = factor(level, 
                       c('Collagen low', 'Collagen int', 'Collagen hi')), 
        regulation = ifelse(significant == 'no', 
                            'ns', 
                            ifelse(estimate > log2(1.25), 
                                   'upregulated', 
                                   ifelse(estimate < -log2(1.25), 
                                          'downregulated', 'ns'))), 
        regulation = factor(regulation, c('upregulated', 'downregulated', 'ns')))
  
  ## identifying genes regulated between the clusters
  
  dge$dge_collagen_int <- dge$lm_signif_results %>% 
    map(filter, 
        level == 'Collagen int', 
        regulation != 'ns')
  
  dge$dge_collagen_hi <- dge$lm_signif_results %>% 
    map(filter, 
        level == 'Collagen hi', 
        regulation != 'ns')

# Extraction of the Entrez IDs and regulation vectors ------
  
  insert_msg('Extracting the significxantly regulated genes')

  ## Entrez IDs for the regulated genes
  
  dge[c('entrez_collagen_int', 
        'entrez_collagen_hi')] <- dge[c('dge_collagen_int', 
                                        'dge_collagen_hi')] %>% 
    map(map,~.x$entrez_id)

  ## regulation vectors
  
  dge[c('regulation_int', 
        'regulation_hi')] <- dge[c('dge_collagen_int', 
                                   'dge_collagen_hi')] %>% 
    map(map, ~set_names(.x$estimate, .x$entrez_id))
  
  ## universe vectors
  
  dge$all_vectors <- dge$test_results %>% 
    map(~.x$anova$entrez_id)
  
# Identification of the common regulated genes ------
  
  insert_msg('Common regulated genes')
  
  dge$entrez_common_int <- dge$entrez_collagen_int %>% 
    reduce(intersect)
  
  dge$symbol_common_int <- dge$anova_signif_results$tcga %>% 
    filter(entrez_id %in% dge$entrez_common_int) %>% 
    .$gene_symbol
  
  dge$entrez_common_hi <- dge$entrez_collagen_hi %>% 
    reduce(intersect)
  
  dge$symbol_common_hi <- dge$anova_signif_results$tcga %>% 
    filter(entrez_id %in% dge$entrez_common_hi) %>% 
    .$gene_symbol
  
# caching the results --------
  
  insert_msg('Caching the results')
  
  save(dge, file = './cache/dge.RData')

# END -----
  
  rm(i)
  
  insert_tail()