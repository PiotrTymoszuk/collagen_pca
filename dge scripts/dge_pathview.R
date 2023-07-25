# Visualizes the regulated components of common pathways activated 
# in the collagen int and high tumors. Working with the 'meta-analysis'
# estimates of regulation calculated by the inverted variance method

  insert_head()
  
# container ------
  
  dge_path <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals -----
  
  insert_msg('Globals')
  
  ## KEGG IDs of the regulated samples
  
  dge_path$kegg_ID <- 
    map2(dge_spia$test, 
         dge_spia$common %>% 
           reduce(union), 
         function(data, name) data[[1]] %>% 
           filter(Name %in% name) %>% 
           .$ID)
  
  ## vectors of differentially regulated genes in at least one cohort
  
  dge_path$collagen_genes[c('int', 'hi')] <- 
    dge[c("dge_collagen_int", "dge_collagen_hi")] %>% 
    map(map, ~.x$entrez_id) %>% 
    map(reduce, union)
  
  ## analysis tables
  
  dge_path$collagen_tbl <- dge$test_results %>% 
    map(~.x$lm) %>% 
    map(filter, 
        level %in% c('Collagen hi', 'Collagen int')) %>% 
    map(blast, level) %>% 
    map(function(data) map2(data, 
                            dge_path$collagen_genes, 
                            ~filter(.x, 
                                    entrez_id %in% .y))) %>% 
    transpose %>% 
    set_names(c('int', 'hi')) %>% 
    map(compress, names_to = 'cohort') %>% 
    map(select, cohort, entrez_id, estimate, se)

# calculation of the meta-regulation estimates -----
  
  insert_msg('Computing meta estimates')
  
  for(i in names(dge_path$collagen_tbl)) {
    
    dge_path$meta_tbl[[i]] <- dge_path$collagen_tbl[[i]] %>% 
      blast(entrez_id) %>% 
      future_map(~safely(metagen)(TE = .x$estimate, 
                                  seTE = .x$se)) %>% 
      map(~.$result) %>% 
      map2_dfr(., names(.), ~tibble(entrez_id = .y, 
                                    estimate = .x$TE.common, 
                                    se = .x$seTE.common))
    
  }
  
  ## regulation vectors, transformation from log2 to identity
  
  dge_path$regulation <- dge_path$meta_tbl %>% 
    map(~set_names(2^.x$estimate, 
                   .x$entrez_id))

# Generating the pathway images -------
  
  insert_msg('Generating the pathway images')
  
  for(i in names(dge_path$kegg_ID)) {
    
    enter_directory(paste('./report/kegg pathviews', i))
    
    dge_path$path_images[[i]] <- dge_path$kegg_ID[[i]] %>% 
      map(~pathview(gene.data = dge_path$regulation[[i]], 
                    pathway.id = .x, 
                    low = list(gene = 'steelblue', 
                               cpd = 'steelblue'), 
                    mid = list(gene = 'white', 
                               cpd = 'white'), 
                    high = list(gene = 'firebrick', 
                                cpd = 'firebrick'), 
                    limit = list(gene = 5, cpd = 5)))
    
    go_proj_directory()
    
  }

# END -----
  
  plan('sequential')
  
  insert_tail()