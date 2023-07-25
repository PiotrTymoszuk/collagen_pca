# This script imports prostate cancer data from GEO (GSE16560, GSE40272, GSE70768 and GSE70769)

  insert_head()
  
# data containers -----
  
  geo_globals <- list()
  geo_data <- list()
    
# definition of folders with data to import and strings to clear from the design data sets ----
  
  geo_globals$study_list <- c('GSE16560', 
                              'GSE40272-GPL15971',
                              'GSE40272-GPL15972', 
                              'GSE40272-GPL15973', 
                              'GSE40272-GPL9497', 
                              'GSE70768', 
                              'GSE70769') %>% 
    paste('./input data', ., sep = '/')
  
  geo_globals$clearing_strings <- c('gleason', ': ', 'major.', 'fup.month', 'diag.yr', 'status.all', 'extreme', 'cancer.percent', 
                                    'age', 'batch', 'minor.', 'ethnicity', 'individual', 'disease state', 'pre-operation treatment', 
                                    'tissue', 'disease free survival', 'sample type', 'tumour', '%', 'iclusterplus group', 
                                    'extra-capsular extension (ece)', 'positive surgical margins (psm)', 'biochemical relapse (bcr)',
                                    'time to bcr (months)', 'tmprss2ERG gene fusion status:', 'at diag', 'psa ', 'clinical ', 'pathology ', 
                                    'total follow up (months)', 'st', 'extra capsular extension (ece)', 'derived data ()')
  
# definition of objects exported to the processor cluster for parallel computing ---- 
  
  geo_globals$cl_export <- c('transpose_tibble', 'read_expression', 'read_design', 'read_annotation', 
                             'untangle_gleason_sum', 'read_study', 'read_excel', 'parse_guess', 'as_tibble', 'set_names', 
                             'geo_globals', 'left_join', 'filter', 'select')

# reading whole genome data ----
  
  cl <- makeCluster(7)
  
  registerDoParallel(cl)
  
  clusterExport(cl, geo_globals$cl_export)
  
  geo_data <- llply(geo_globals$study_list, 
                    function(x) read_study(x, 
                                           merge_clinical = T, 
                                           geo_globals$clearing_strings), 
                    .parallel = T) %>% 
    set_names(c('GSE16560', 
                'GSE40272-GPL15971',
                'GSE40272-GPL15972', 
                'GSE40272-GPL15973', 
                'GSE40272-GPL9497', 
                'GSE70768', 
                'GSE70769'))
  
  stopCluster(cl)
  
# END ----
  
  insert_tail()

  
  
