# Identification of differentially regulated genes 
# between the Collagen clusters
# Differentially regulated genes are identified by two-tailed T test
# corrected for multiple testing with the Benjamini-Hochberg method
# Regulation cutoff: 1.2-fold

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
                            'Collagen hi')), 
        clust_id = droplevels(clust_id))
  
# Serial t test -------
  
  insert_msg('Serial testing')
  
  dge$test_results <- list(data = dge$analysis_tbl, 
                           variables = dge$variables) %>% 
    pmap(test_two_groups, 
         split_fct = 'clust_id', 
         type = 't', 
         adj_method = 'BH', 
         .parallel = TRUE)
  
# Annotation of the results with Entrez IDs -------
  
  insert_msg('Annotation of the results')
  
  dge$test_results <- 
    map2(dge$test_results, 
         dge$annotation, 
         ~mutate(.x, 
                 gene_symbol = response, 
                 entrez_id = exchange(gene_symbol, 
                                      .y, 
                                      key = 'gene_symbol', 
                                      value = 'entrez_id')))


# Formatting the results and identifying differentially regulated genes ------
  
  insert_msg('Identification of differentially regulated genes')
  
  # pFDR < 0.05 and at least weak regulation as assessed by Cohen's d
  
  dge$test_results <- dge$test_results %>% 
    map2(., dge$analysis_tbl, 
         ~mutate(.x, 
                 n = nrow(.y), 
                 se = abs(estimate/stat), 
                 sd = se * sqrt(n), 
                 eff_size = abs(estimate/sd), 
                 regulation = ifelse(significant == 'yes' & 
                                       eff_size > 0.2, 
                                     ifelse(estimate > 0, 
                                            'upregulated', 'downregulated'), 
                                     'ns'), 
                 regulation = factor(regulation, 
                                     c('upregulated', 'downregulated', 'ns')))) %>% 
    map(filter, 
        !is.na(entrez_id), 
        !is.na(regulation))
  
  ## significant results
  
  dge$signif_results <- dge$test_results %>% 
    map(filter, 
        regulation != 'ns')

# Extraction of the Entrez IDs and regulation vectors ------
  
  insert_msg('Extracting the significantly regulated genes')

  ## Entrez IDs for the regulated genes
  
  dge$dge_entrez <- dge$signif_results %>% 
    map(~.x$entrez_id)

  ## regulation vectors
  
  dge$regulation_vectors <- dge$signif_results %>% 
    map(~set_names(.x$estimate, .x$entrez_id))

  ## universe vectors
  
  dge$all_vectors <- dge$test_results %>% 
    map(~.x$entrez_id)
  
# caching the results --------
  
  insert_msg('Caching the results')
  
  save(dge, file = './cache/dge.RData')

# END -----
  
  insert_tail()