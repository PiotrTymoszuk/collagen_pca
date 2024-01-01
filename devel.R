# Some experiments with metabolite fluxes

  library(BiGGR)

  data('Recon2')

  gene.info <- extractGeneAssociations(Recon2)
  
  reaction_ids <- names(gene.info)
  
  plan('multisession')
  
  res <- dge$significant_entrez$upregulated$tcga[1:1000] %>% 
    unname %>% 
    future_map(function(gene) map_lgl(gene.info, 
                                      ~stri_detect(.x, fixed = paste0('(', gene))))
  
  plan('sequential')
  
  res %>% 
    map_dbl(sum)

  res %>% 
    map(function(x) if(sum(x) > 0) x else NULL) %>% 
    compact
  
  library(tidyverse)
  
  tcga_data <- dge$test_results$tcga$lm %>% 
    filter(level == 'clust_idCollagen hi') %>% 
    mutate(gene_symbol = response, 
           entrez_id = dge$annotation_vec$tcga[gene_symbol], 
           regulation = ifelse(significant == 'no', 
                               'ns', 
                               ifelse(estimate > 0, 
                                      'upregulated', 'downregulated'))) %>% 
    select(gene_symbol, entrez_id, 
           estimate, 
           se, 
           lower_ci, upper_ci, 
           p_adjusted, 
           regulation)
  
  save(tcga_data, file = 'tcga_data.RData')
  
  
  