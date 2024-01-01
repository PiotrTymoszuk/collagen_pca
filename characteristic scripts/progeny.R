# Analysis of differential modulation of signaling pathway by Progeny.
# The analysis is carried out essentially as described in the decoupleR's
# vignette (https://saezlab.github.io/decoupleR/articles/pw_bk.html).
#
# Effect sizes (Cohen's d) of whole-transcriptome gene regulation are derived 
# from the DGE analysis results with two-tailed T test. 
# The pathway activity statistic is the 'multivariable linear model score'
# introduced by the package Authors (sum of t scores of the linear model
# estimates that define hyperplane of gene expression?).
#
# Significantly up- and downregulated pathways are defined by the MLM score
# significance corrected for multiple testing separately for each cohort.
# Common significantly activated or inhibited pathways are defined as pathways
# found significantly up- or downregulated in at least five cohorts excluding 
# GSE16560.

  insert_head()
  
# container -------
  
  ana_progeny <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# The data base: resorting to the cache, it takes a while! ---------
  
  insert_msg('Pathway - gene database')
  
  if(file.exists('./data/progeny.RData')) {
    
    load('./data/progeny.RData')

  } else {
    
    progeny <- get_progeny(top = 500)
    
    save(progeny, file = './data/progeny.RData')
    
  }
  
  ana_progeny$progeny <- progeny
  
# Gene regulation estimates --------
  
  insert_msg('Gene expression regulation estimates')
  
  ## in a matrix with genes in rows and samples in columns 
  ## (in out case just one!)
  
  ana_progeny$estimates <- ana_dge$test %>% 
    map(column_to_rownames, 'gene_symbol') %>% 
    map(select, effect_size) %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(as.matrix)
  
# modeling --------
  
  insert_msg('Modeling')

  ana_progeny$test <- ana_progeny$estimates %>%
    future_map(run_mlm,
               network = ana_progeny$progeny,
               .source = 'source',
               .target = 'target',
               .mor = 'weight',
               minsize = 5,
               .options = furrr_options(seed = TRUE))
  
# FDR correction and identification of significantly regulated pathways -------
  
  insert_msg('Formatting and significant effects')
  
  ana_progeny$test <- ana_progeny$test %>% 
    map(re_adjust, 'p_value', method = 'BH') %>% 
    map(mutate, 
        regulation = ifelse(p_adjusted >= 0.05, 
                            'ns', 
                            ifelse(score > 0, 'activated', 
                                   ifelse(score < 0, 'inhibited', 'ns'))), 
        regulation = factor(regulation, 
                            c('activated', 'inhibited', 'ns')))
  
  ## significant pathways in single cohorts
  
  ana_progeny$significant <- ana_progeny$test %>% 
    map(filter, regulation %in% c('activated', 'inhibited')) %>%
    map(blast, regulation) %>% 
    transpose %>% 
    map(map, ~.x$source)
  
  ## common regulated pathways
  
  ana_progeny$common_significant <- ana_progeny$significant %>% 
    map(~.x[names(.x) != 'gse16560']) %>% 
    map(shared_features, m = 5)

# Bubble heat plot with the common regulated pathways ------
  
  insert_msg('Bubble plot with the common regulated pathways')
  
  ana_progeny$cmm_plot <-
    regulation_bubble(data_lst = ana_progeny$test,
                      variables = reduce(ana_progeny$common_significant, union), 
                      label_variable = 'source', 
                      regulation_variable = 'score', 
                      status_variable = 'regulation', 
                      plot_title = 'PROGENy signaling, collagen hi vs low', 
                      plot_subtitle = paste('Pathways differentially modulated', 
                                            'in at least 5 cohorts'), 
                      size_lab = 'MLM score', 
                      limits = c(0, 25), 
                      breaks = seq(0, 25, 
                                   by = 5), 
                      max_size = 5) + 
    scale_y_discrete(labels = progeny_labeller)
  
# END --------
  
  ana_progeny$progeny <- NULL
  ana_progeny$estimates <- NULL
  
  ana_progeny <- compact(ana_progeny)
  
  plan('sequential')
  
  insert_tail()