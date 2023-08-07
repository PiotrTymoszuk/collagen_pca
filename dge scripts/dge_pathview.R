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
  
  dge_path$kegg_ID <- dge_spia$test %>% 
    map(filter, Name %in% unique(unlist(dge_spia$common))) %>% 
    map(~.x$ID) %>% 
    reduce(union)

  ## vectors of differentially regulated genes in at least one cohort
  
  dge_path$collagen_genes <- dge$dge_entrez %>% 
    reduce(union)

  ## analysis tables
  
  dge_path$collagen_tbl <- dge$test_results %>% 
    map(filter, entrez_id %in% dge_path$collagen_genes) %>% 
    compress(names_to = 'cohort') %>% 
    select(cohort, entrez_id, estimate, se)

# calculation of the meta-regulation estimates -----
  
  insert_msg('Computing meta estimates')

  dge_path$meta_tbl <- dge_path$collagen_tbl %>% 
    blast(entrez_id) %>% 
    future_map(~safely(metagen)(TE = .x$estimate, 
                                seTE = .x$se)) %>% 
    map(~.$result) %>% 
    map2_dfr(., names(.), 
             ~tibble(entrez_id = .y, 
                     estimate = .x$TE.common, 
                     se = .x$seTE.common))
  
  ## regulation vectors, transformation from log2 to identity
  
  dge_path$regulation <- 
    set_names(2^dge_path$meta_tbl$estimate, 
              dge_path$meta_tbl$entrez_id)
  
# Generating the pathway images -------
  
  insert_msg('Generating the pathway images')

  enter_directory('./report/kegg pathviews')
  
  dge_path$path_images <- dge_path$kegg_ID %>% 
    map(~pathview(gene.data = dge_path$regulation, 
                  pathway.id = .x, 
                  low = list(gene = 'steelblue', 
                             cpd = 'steelblue'), 
                  mid = list(gene = 'white', 
                             cpd = 'white'), 
                  high = list(gene = 'firebrick', 
                              cpd = 'firebrick'), 
                  limit = list(gene = 5, cpd = 5)))
  
  go_proj_directory()

# END -----
  
  plan('sequential')
  
  insert_tail()