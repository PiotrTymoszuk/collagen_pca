# Differential gene expression in the collagen clusters.
#
# Investigated by FDR-corrected two-tailed T test and Cohen's d effect size 
# statistic. Differentially regulated genes are defined by pFDR < 0.05 and 
# d >= 0.2, which indicates a significant, at last weak effect size.
# Common differentially regulated genes are shared by at least five cohorts 
# except of  GSE16560.

  insert_head()
  
# container -------
  
  ana_dge <- list()
  
# testing -------
  
  insert_msg('Testing')
  
  ana_dge$test <- list(data = ana_globals$expression, 
                       variables = ana_globals$genes) %>% 
    pmap(test_two_groups, 
         split_fct = 'clust_id', 
         type = 't', 
         adj_method = 'BH', 
         .parallel = TRUE)
  
# Annotation of the results -------
  
  insert_msg('Annotation')
  
  ana_dge$test <- 
    map2(ana_dge$test, 
         ana_globals$annotation, 
         ~mutate(.x, 
                 gene_symbol = response, 
                 entrez_id = exchange(gene_symbol, 
                                      dict = .y, 
                                      key = 'gene_symbol', 
                                      value = 'entrez_id')))
  
# Significant effects ---------
  
  insert_msg('Differentially regulated genes')
  
  ana_dge$test <- ana_dge$test %>% 
    map(mutate, 
        regulation = ifelse(p_adjusted >= 0.05, 
                            'ns', 
                            ifelse(effect_size >= 0.2, 
                                   'upregulated', 
                                   ifelse(effect_size <= -0.2, 
                                          'downregulated', 'ns'))), 
        regulation = factor(regulation, 
                            c('upregulated', 'downregulated', 'ns')))
  
  ## significantly up- and downregulated genes
  
  ana_dge$significant <- ana_dge$test %>% 
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>% 
    map(blast, regulation) %>% 
    transpose %>% 
    map(map, ~.x$gene_symbol)
  
  ## common significant genes
  
  ana_dge$common_significant <- ana_dge$significant %>% 
    map(~.x[names(.x) != 'gse16560']) %>% 
    map(shared_features, m = 5)
  
# Caching the results ---------
  
  insert_msg('Caching the results')
  
  save(ana_dge, file = './cache/ana_dge.RData')
  
# END ---------
  
  insert_tail()
  