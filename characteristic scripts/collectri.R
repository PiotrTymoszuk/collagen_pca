# Analysis of differential activity of regulons (genes activated or inhibited 
# by a common transcription factor).
# The analysis is carried out essentially as described in the decoupleR's
# vignette (https://saezlab.github.io/decoupleR/articles/pw_bk.html).
#
# Estimates of whole-transcriptome gene regulation are derived from the DGE 
# analysis results with two-tailed T test. 
# The pathway activity statistic is the 'univariable linear model score'
# introduced by the package Authors.
#
# Significantly up- and downregulated pathways are defined by the LM score
# significance corrected for multiple testing separately for each cohort.
# Common significantly activated or inhibited pathways are defined as pathways
# found significantly up- or downregulated in at least five cohorts excluding 
# GSE16560.

  insert_head()
  
# container --------
  
  ana_collectri <- list()
  
# parallel backend --------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# regulon database ------
  
  insert_msg('Regulon database')
  
  if(file.exists('./data/collectri.RData')) {
    
    load('./data/collectri.RData')
    
  } else {
    
    collectri <- get_collectri(split_complexes = FALSE)
    
    save(collectri, file = './data/collectri.RData')
    
  }
  
  ana_collectri$collectri <- collectri
  
# GO term database for regulon classification -------
  
  insert_msg('GO database')
  
  ana_collectri$go_db <- godata('org.Hs.eg.db', 
                                ont= 'BP', 
                                computeIC = FALSE)
  
# Gene regulation estimates --------
  
  insert_msg('Gene expression regulation estimates')
  
  ## in a matrix with genes in rows and samples in columns 
  ## (in out case just one!)
  
  ana_collectri$estimates <- ana_dge$test %>% 
    map(column_to_rownames, 'gene_symbol') %>% 
    map(select, effect_size) %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(as.matrix)
  
# modeling ---------
  
  insert_msg('Modeling')
  
  ana_collectri$test <- ana_collectri$estimates %>%
    future_map(run_ulm,
               network = ana_collectri$collectri,
               .source = 'source',
               .target = 'target',
               .mor = 'mor',
               minsize = 5,
               .options = furrr_options(seed = TRUE))
  
# FDR correction and significant effects -------
  
  insert_msg('Formatting and significant effects')
  
  ana_collectri$test <- ana_collectri$test %>% 
    map(re_adjust, 'p_value') %>% 
    map(mutate, 
        regulation = ifelse(p_adjusted >= 0.05, 'ns', 
                            ifelse(score > 0, 'activated', 
                                   ifelse(score < 0, 'inhibited', 'ns'))), 
        regulation = factor(regulation, c('activated', 'inhibited', 'ns')))
  
  ## significance in single cohorts

  ana_collectri$significant <- ana_collectri$test %>% 
    map(filter, regulation %in% c('activated', 'inhibited')) %>% 
    map(blast, regulation) %>% 
    transpose %>% 
    map(map, ~.x$source)
  
  ## regulons shared by at least five cohorts

  ana_collectri$common_significant <- ana_collectri$significant %>% 
    map(~.x[names(.x) != 'gse16560']) %>% 
    map(shared_features, m = 5)
  
# Top 20 regulons -------
  
  insert_msg('Top regulons')
  
  ## out of the common modulated ones. The criterion is the average score
  
  ana_collectri$mean_scores <- ana_collectri$test %>% 
    map(filter, 
        source %in% reduce(ana_collectri$common_significant, union)) %>% 
    compress(names_to = 'cohort') %>% 
    group_by(source) %>% 
    summarise(score = mean(score)) %>% 
    ungroup %>% 
    mutate(regulation = ifelse(score > 0, 'activated', 'inhibited'))
  
  ana_collectri$top_regulons <- ana_collectri$mean_scores %>% 
    group_by(regulation) %>% 
    top_n(n = 20, abs(score)) %>% 
    ungroup
  
# Bubble plot for the top regulons -------
  
  insert_msg('Bubble plot')
  
  ana_collectri$cmm_plot <- 
    regulation_bubble(data_lst = ana_collectri$test,
                      variables = ana_collectri$top_regulons$source, 
                      label_variable = 'source', 
                      regulation_variable = 'score', 
                      status_variable = 'regulation', 
                      plot_title = 'collecTRI regulons, collagen hi vs low', 
                      plot_subtitle = paste('Top 20 regulons differentially', 
                                            'modulated in at least 5 cohorts'), 
                      size_lab = 'LM score')
  
# END ------
  
  ana_collectri$collectri <- NULL
  ana_collectri$mean_scores <- NULL
  ana_collectri$go_db <- NULL
  ana_collectri$estimates <- NULL
  
  ana_collectri <- compact(ana_collectri)
  
  plan('sequential')
  
  insert_tail()