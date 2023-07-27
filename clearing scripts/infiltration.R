# Calculates infiltration of non-malignant cell in the tumor samples
# Two methods are applied: xCell and MCP Counter

  insert_head()
  
# container ------
  
  infil <- list()
  
# globals: analysis tables -------
  
  insert_msg('Analysis tables')
  
  infil$analysis_tbl <- study_data %>% 
    map(~.x$expression) %>% 
    map(function(x) if('tissue' %in% names(x)) filter(x, tissue == 'tumor') else x)
  
  
  infil$variables <- study_data %>% 
    map(~.x$annotation) %>% 
    map(filter, 
        !duplicated(symbol)) %>% 
    map(~.x$symbol)
  
  infil$analysis_tbl <- 
    map2(infil$analysis_tbl, 
         infil$variables, 
         ~.x[c('patient_id', .y)]) %>% 
    map(filter, !is.na(patient_id)) %>% 
    map(column_to_rownames, 'patient_id') %>%
    map(t)
  
  infil$analysis_tbl[c('GSE16560', 
                       'GSE40272', 
                       'GSE70768', 
                       'GSE70769')] <- infil$analysis_tbl[c('GSE16560', 
                                                            'GSE40272', 
                                                            'GSE70768', 
                                                            'GSE70769')] %>% 
    map(~2^.x)
  
  infil$analysis_tbl$tcga <- 2^infil$analysis_tbl$tcga - 1
  
  ## parallel backend
  
  plan('multisession')
  
# X cell deconvolution ------
  
  insert_msg('Xcell deconvolution')
  
  infil$xcell <- list(gene_expression = infil$analysis_tbl,
                      arrays = c(TRUE, TRUE, TRUE, TRUE, FALSE)) %>%
    future_pmap(safely(deconvolute),
                method = 'xcell',
                expected_cell_types = c('B cell', 
                                        'T cell CD8+', 
                                        'T cell CD4+ (non-regulatory)', 
                                        'T cell regulatory (Tregs)', 
                                        'Myeloid dendritic cell', 
                                        'Macrophage', 
                                        'Monocyte', 
                                        'Neutrophil', 
                                        'NK cell', 
                                        'Cancer associated fibroblast', 
                                        'Endothelial cell')) %>% 
    map(~.x$result) %>% 
    compact
  
# MCP counter deconvolution ------
  
  insert_msg('MCP counter deconvolution')
  
  infil$mcp_counter <- list(gene_expression = infil$analysis_tbl,
                            arrays = c(TRUE, TRUE, TRUE, TRUE, FALSE)) %>%
    future_pmap(safely(deconvolute),
                method = 'mcp_counter') %>% 
    map(~.x$result) %>% 
    compact
  
# saving the results ------
  
  insert_msg('saving the results')
  
  infil <- infil[c('xcell', 'mcp_counter')]
  
  save(infil, file = './data/infiltration.RData')
  
# END -----
  
  plan('sequential')
  
  insert_tail()