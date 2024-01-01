# Visualizes the regulated components of common pathways activated 
# in the collagen int and high tumors. Working with the average 
# estimates of regulation estimates

  insert_head()
  
  go_proj_directory()
  
# container ------
  
  ana_path <- list()

# globals -----
  
  insert_msg('Globals')
  
  ## KEGG IDs of the regulated samples
  
  ana_path$kegg_ID <- ana_spia$test %>% 
    map(filter, Name %in% unique(unlist(ana_spia$common_significant))) %>% 
    map(~.x$ID) %>% 
    reduce(union)

  ## vectors of differentially regulated genes in at least one cohort
  
  ana_path$genes <- ana_dge$significant %>%
    unlist %>% 
    unique

  ## analysis tables
  
  ana_path$estimates <- ana_dge$test %>% 
    map(filter, gene_symbol %in% ana_path$genes) %>% 
    compress(names_to = 'cohort') %>% 
    select(entrez_id, gene_symbol, estimate) %>% 
    group_by(gene_symbol, entrez_id) %>% 
    summarise(estimate = mean(estimate)) %>% 
    ungroup
  
  ## vectors with identity regulation estimates
  
  ana_path$regulation <- 
    set_names(2^ana_path$estimates$estimate, 
              ana_path$estimates$entrez_id)

# Generating the pathway images -------
  
  insert_msg('Generating the pathway images')

  enter_directory('./report/kegg pathviews')
  
  ana_path$path_images <- ana_path$kegg_ID %>% 
    map(~pathview(gene.data = ana_path$regulation, 
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
  
  ana_path$estimates <- NULL
  
  ana_path <- compact(ana_path)
  
  insert_tail()