# Identification of differential regulated biochemical reactions 
# in the Recon2 BiGG model. 
# SEM brought back to the linear scale (from log2) by the first-differential 
# equation log(2) * SEM(log2-Expression) * 2^log2-Expression
# FDR multiple testing correction

  insert_head()
  
# container -----
  
  meta <- list()
  
# globals ------
  
  insert_msg('Globals')
  
  ## Tables with regulation estimates and their errors
  ## annotation with Entrez ID
  
  meta$analysis_tbl <- dge$test_results %>% 
    map(mutate, 
        se = log(2) * se * 2^estimate, 
        estimate = 2^estimate) %>% 
    map(~.x[c('gene_symbol', 'entrez_id', 'estimate', 'se')])

# construction of the models ------
  
  insert_msg('Construction of the SBML models')
  
  set.seed(1234)
  
  meta$models <- meta$analysis_tbl %>% 
    map(~build_geneSBML(x = set_names(.x$estimate, .x$entrez_id), 
                        err = set_names(.x$se, .x$entrez_id), 
                        database = Recon2D, 
                        or_fun = 'mean', 
                        and_fun = 'min', 
                        x_default = 1, 
                        err_method = 'mc', 
                        n_iter = 2010, 
                        ci_method = 'bca', 
                        .parallel = TRUE))
  
# Identification of significantly regulated reactions -------
  
  insert_msg('Significantly regulated reactions')
  
  ## regulation: log-fold as well!
  
  meta$regulation <- meta$models %>% 
    map(components, 'regulation') %>% 
    map(mutate, 
        significant = ifelse(p_adjusted < 0.05, 'yes', 'no'),
        regulation = ifelse(significant == 'no', 
                            'ns', 
                            ifelse(fold_reg > 1, 'activated', 
                                   ifelse(fold_reg < 1, 'inhibited', 'ns'))), 
        regulation = factor(regulation, c('activated', 'inhibited', 'ns'))) %>% 
    map(mutate, 
        log_fold_reg = log2(fold_reg), 
        log_lower_ci = log2(lower_ci), 
        log_upper_ci = log2(upper_ci))

  ## significant reactions
  
  meta$signif_reactions <- meta$regulation %>% 
    map(filter, regulation != 'ns') %>% 
    map(blast, regulation) %>% 
    map(map, ~.x$react_id) %>% 
    transpose

# caching the results --------
  
  insert_msg('Caching the results')
  
  save(meta, file = './cache/meta.RData')
  
# END -----
  
  insert_tail()