# GO enrichment for genes differentially regulated between the clusters
# Significant effects defined by pFDR < 0.05 and enrichment odds ratio (OR) of 
# at least 1.44, i.e. at least weak effect size. Common significant effects 
# are shared by at least 5 cohorts except GSE16560.

  insert_head()
  
# container -------
  
  ana_go <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# GO annotation data base -------
  
  insert_msg('GO database')
  
  ## which will be later used for classification
  ## of significantly enriched terms
  
  ana_go$go_db <- godata('org.Hs.eg.db', ont = "BP", computeIC = FALSE)
  
# Vectors of Entrez ID of the significantly regulated genes -------
  
  insert_msg('Gene ID vectors')
  
  ## all genes
  
  ana_go$all <- ana_globals$annotation %>% 
    map(filter, 
        !is.na(entrez_id), 
        !duplicated(entrez_id), 
        entrez_id != '') %>% 
    map(~.x$entrez_id) %>% 
    map(unname)
  
  ## up- and downregulated genes
  
  ana_go$dge <- ana_dge$significant %>% 
    map(map, 
        mapIds, 
        x = org.Hs.eg.db, 
        keytype = 'SYMBOL', 
        column = 'ENTREZID')
  
  ana_go$dge <- ana_go$dge %>% 
    map(map, ~.x[!is.na(.x)]) %>% 
    map(map, ~.x[!duplicated(.x)]) %>% 
    map(map, ~.x[.x != '']) %>% 
    map(map, unname)
  
# Enrichment analysis ---------
  
  insert_msg('Enrichment testing')
  
  for(i in names(ana_go$dge)) {
    
    ana_go$test[[i]] <- 
      list(de = ana_go$dge[[i]], 
           universe = ana_go$all) %>% 
      future_pmap(GOana, 
                  ontology = 'BP', 
                  adj_method = 'BH', 
                  .options = furrr_options(seed = TRUE))
    
  }
  
# Significantly enriched terms -------
  
  insert_msg('Significantly enriched terms')
  
  ## I'm extracting GO IDs as well, which will be later utilized
  ## in clustering
  
  ana_go$significant <- ana_go$test %>% 
    map(map, 
        filter, 
        p_adjusted < 0.05, 
        or >= 1.44) %>% 
    map(map, ~.x$term)
  
  ana_go$significant_id <- ana_go$test %>% 
    map(map, 
        filter, 
        p_adjusted < 0.05, 
        or >= 1.44) %>% 
    map(map, ~.x$go_id)
  
  ## common significant terms
  
  ana_go$common_significant <- ana_go$significant %>% 
    map(~.x[names(.x) != 'gse165060']) %>% 
    map(shared_features, m = 5)
  
  ana_go$common_significant_id <- ana_go$significant_id %>% 
    map(~.x[names(.x) != 'gse165060']) %>% 
    map(shared_features, m = 5)
  
# similarity between the common enriched GOs -------
  
  insert_msg('Wang similatrity between common GOs')
  
  ## lexicon of GO terms and their IDs used later for annotation
  ## of the clustering results
  
  ana_go$lexicon <- ana_go$test %>% 
    map(map, select, go_id, term) %>% 
    map(reduce, rbind) %>% 
    reduce(rbind)
  
  ana_go$lexicon <- ana_go$lexicon %>% 
    filter(!duplicated(go_id), 
           !duplicated(term))

  ## hierarchical clustering of the common regulated GOS
  
  set.seed(12345)

  ana_go$go_clust <- 
    list(GOs = ana_go$common_significant_id, 
         k = c(3, 3)) %>% 
    pmap(clust_gos, 
         semData = ana_go$go_db, 
         fun = kcluster, 
         clust_fun = 'pam')

  ## GO cluster assignment, annotation with human-friendly names
  
  ana_go$go_assignment <- ana_go$go_clust %>% 
    map(extract, 'assignment') %>% 
    map(set_names, c('go_id', 'clust_id')) %>% 
    map(left_join, 
        ana_go$lexicon, 
        by = 'go_id')
  
  ## setting descriptive names of the clusters
  
  ana_go$go_clust_desc$upregulated <- 
    c('1' = 'ECM\ncytoskeleton\nadhesion\nmotilty\nAg and T cells', 
      '2' = 'BMP and TGF-\u03B2\nMAPK and WNT\nTLR, IFN and T cells', 
      '3' = 'development\nangiogenesis')
  
  ana_go$go_clust_desc$downregulated <- 
    c('1' = 'RNA processing', 
      '2' = 'mitochondrion\nribosome', 
      '3' = 'oxidative\nmetabolism')
  
# Visualization of classification of the common enriched GOs -------
  
  insert_msg('Visualization of the common enriched GOs')
  
  ## MDS layouts
  
  ana_go$plots <- ana_go$go_clust %>% 
    map(plot,
        type = 'data', 
        cust_theme = globals$common_theme)
  
  ## titles and cluster descriptions
  
  ana_go$plots <-
    list(x = ana_go$plots, 
         y = c('GO enrichment, upregulated in Collagen hi', 
               'GO enrichment, downregulated in Collagen hi'), 
         z = ana_go$go_clust_desc) %>% 
    pmap(function(x, y, z) x + 
           labs(title = y, 
                subtitle = 'GO terms significant in at least 5 cohorts, Wang similarity') + 
           scale_fill_manual(values = c('1' = 'cornsilk3', 
                                        '2' = 'plum4', 
                                        '3' = 'orangered3'), 
                             labels = z))
  

# Caching the results ---------
  
  insert_msg('Caching the results')
  
  ana_go <- ana_go[c("test", "significant", 
                     "common_significant", "go_clust", 
                     "go_assignment", "go_clust_desc", 
                     "plots")]
  
  save(ana_go, file = './cache/ana_go.RData')
  
# END ------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()
  