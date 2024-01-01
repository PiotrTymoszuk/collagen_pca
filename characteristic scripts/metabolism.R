# Monte Carlo modeling of activity of metabolic reactions with biggrExtra. 
# Whole-transcriptome regulation estimates and their SEM are used to predict 
# significantly up- and downregulated metabolic reactions in the collagen high 
# versus collagen low cluster. Multiple testing adjustment with FDR separately 
# in each cohort.
# Enrichment of metabolic Recon subsystems in significantly activated and 
# inhibited reactions is investigated by comparison with 100000 random draws from
# the entire reaction pool. Significant enrichment is defined by pFDR < 0.05 
# and odds ratio (OR) >= 1.44, which indicates a significant enrichment with at 
# least weak effect size. Common enriched metabolic subsystems are defined as
# subsystems significantly enriched in activated or inhibited reactions in at 
# least five cohorts except GSE16560

  insert_head()
  
# container --------
  
  ana_meta <- list()
  
# gene regulation estimates -------
  
  insert_msg('Gene expression regulation estimates')
  
  ana_meta$estimates <- ana_dge$test %>% 
    map(filter, 
        !is.na(entrez_id), 
        !duplicated(entrez_id), 
        entrez_id != '') %>% 
    map(~set_names(.x$estimate, .x$entrez_id))
  
  ana_meta$errors <- ana_dge$test %>% 
    map(filter, 
        !is.na(entrez_id), 
        !duplicated(entrez_id), 
        entrez_id != '') %>% 
    map(~set_names(.x$error, .x$entrez_id))
  
# modeling -------
  
  insert_msg('Modeling')
  
  ana_meta$models <- 
    list(x = ana_meta$estimates, 
         err = ana_meta$errors) %>% 
    pmap(build_geneSBML, 
         database = Recon2D, 
         scale = 'log2', 
         or_fun = 'mean', 
         and_fun = 'min', 
         x_default = 1, 
         err_method = 'mc', 
         n_iter = 3000,
         mc_estimate = 'mean', 
         ci_method = 'bca',
         save_memory = TRUE, 
         .parallel = TRUE)
  
# Subsystem enrichment ------
  
  insert_msg('Subsystem enrichment')
  
  ana_meta$subsystems <- ana_meta$models %>% 
    map(suba, 
        signif_type = 'fdr', 
        method = 'simulation', 
        n_iter = 100000, 
        .parallel = TRUE)
  
  ## re-adjustment, I'm not interested at generally enrichment subsystems
  
  ana_meta$subsystems <- ana_meta$subsystems %>% 
    map(filter, status %in% c('activated', 'inhibited')) %>% 
    map(re_adjust, 'p_value')
  
# Significantly enriched metabolic subsystems -------
  
  insert_msg('Significantly enriched metabolic subsystems')
  
  # subsystems significantly enriched in activated and inhibited reactions
  
  ana_meta$significant_subs <- ana_meta$subsystems %>% 
    map(filter, 
        p_adjusted < 0.05, 
        OR >= 1.44) %>% 
    map(blast, status) %>% 
    transpose %>% 
    map(map, ~.x$subsystem) %>% 
    map(map, as.character)
  
  ## common enriched subsystems
  
  ana_meta$common_significant_subs <- ana_meta$significant_subs %>% 
    map(~.x[names(.x) != 'gse16560']) %>% 
    map(shared_features, m = 5)
  
# Caching the results --------
  
  insert_msg('Caching the results')
  
  ana_meta <- 
    ana_meta[c("models", "subsystems", 
               "significant_subs", "common_significant_subs")]
  
  save(ana_meta, file = './cache/ana_meta.RData')
  
# END -------
  
  insert_tail()