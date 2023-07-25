# Identification of differential regulated biochemical reactions 
# in the Recon2 BiGG model. 
# SEM brought back to the linear scale (from log2) by the first-differential 
# equation log(2) * SEM(log2-Expression) * 2^log2-Expression
# No multiple-testing correction!!!

  insert_head()
  
# container -----
  
  meta <- list()
  
# globals ------
  
  insert_msg('Globals')
  
  ## Tables with regulation estimates and their errors
  ## annotation with Entrez ID
  
  meta$analysis_tbl <- dge$test_results %>% 
    map(~.x$lm) %>% 
    map(filter, level != '(Intercept)') %>% 
    map(mutate, 
        level = stri_extract(level, regex = 'int|hi'), 
        level = factor(level, c('int', 'hi')), 
        se = log(2) * se * 2^estimate, 
        estimate = 2^estimate) %>% 
    map(~.x[c('level', 'gene_symbol', 'entrez_id', 'estimate', 'se')])

  meta$analysis_tbl <- meta$analysis_tbl %>% 
    map(blast, level) %>% 
    transpose

# construction of the models ------
  
  insert_msg('Construction of the SBML models')
  
  set.seed(1234)
  
  meta$models$int <- meta$analysis_tbl$int %>% 
    map(~build_geneSBML(x = set_names(.x$estimate, .x$entrez_id), 
                        err = set_names(.x$se, .x$entrez_id), 
                        database = Recon2D, 
                        or_fun = 'mean', 
                        and_fun = 'min', 
                        x_default = 1, 
                        err_method = 'mc', 
                        n_iter = 1010, 
                        ci_method = 'bca', 
                        .parallel = TRUE))
  
  meta$models$hi <- meta$analysis_tbl$hi %>% 
    map(~build_geneSBML(x = set_names(.x$estimate, .x$entrez_id), 
                        err = set_names(.x$se, .x$entrez_id), 
                        database = Recon2D, 
                        or_fun = 'mean', 
                        and_fun = 'min', 
                        x_default = 1, 
                        err_method = 'mc', 
                        n_iter = 1010, 
                        ci_method = 'bca', 
                        .parallel = TRUE))
  
# Identification of significantly regulated reactions -------
  
  insert_msg('Significantly regulated reactions')
  
  meta$regulation <- meta$models %>% 
    map(map, components, 'regulation') %>% 
    map(map, 
        mutate, 
        significant = ifelse(p_adjusted < 0.05, 'yes', 'no'), 
        regulation = ifelse(significant == 'no', 
                            'ns', 
                            ifelse(fold_reg > 1, 'activated', 
                                   ifelse(fold_reg < 1, 'inhibited', 'ns'))), 
        regulation = factor(regulation, c('activated', 'inhibited', 'ns')))
  
  meta$signif_reactions <- meta$regulation %>% 
    map(map, filter, regulation != 'ns') %>% 
    map(map, blast, regulation) %>% 
    map(map, map, ~.x$react_id)

# Identification of significant reactions common for all studies -----
  
  insert_msg('Common significant reactions')
  
  ## all cohorts except of GSE40272, where no significant regulation
  ## was found
  
  for(i in names(meta$signif_reactions)) {
    
    meta$common_reactions[[i]] <- 
      meta$signif_reactions[[i]][c('GSE16560', 'GSE70768', 'GSE70769', 'tcga')] %>% 
      transpose %>% 
      map(reduce, intersect)
    
  }
  
# caching the results --------
  
  insert_msg('Caching the results')
  
  save(meta, file = './cache/meta.RData')
  
# END -----
  
  rm(i)
  
  insert_tail()