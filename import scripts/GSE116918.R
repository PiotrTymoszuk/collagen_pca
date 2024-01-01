# Import of the GSE116918 data set

  insert_head()
  
# container ------
  
  gse116918 <- list()
  
# reading from GEO ------
  
  insert_msg('Reading from GEO')
  
  gse116918$raw <- getGEO(GEO = 'GSE116918', destdir = './data/GSE116918')
  
# clinical information -------
  
  insert_msg('Clinical information')
  
  gse116918$clinic <- pData(gse116918$raw[[1]]) %>% 
    as_tibble
  
  gse116918$clinic <- gse116918$clinic %>% 
    transmute(sample_id = geo_accession, 
              patient_id = title, 
              tissue_type = factor('tumor', c('normal', 'tumor')), 
              relapse = as.numeric(`bcr event (1=yes, 0=no):ch1`), 
              rfs_months = as.numeric(`follow-up time (bcr, months):ch1`),
              metastasis = as.numeric(`met event (1=yes, 0=no):ch1`), 
              mfs_months = as.numeric(`follow-up time (met, months):ch1`), 
              gleason_sum = as.numeric(`gleason grade:ch1`), 
              gleason_simple = cut(gleason_sum, 
                                   c(-Inf, 6, 7, Inf), 
                                   c('5 - 6', '7', '8+')), 
              gleason_sum = factor(gleason_sum), 
              age = as.numeric(`patient age (years):ch1`),
              psa_diagnosis = as.numeric(`psa (ng/ml):ch1`), 
              pt_stage = stri_extract(`t-stage:ch1`, regex = 'T\\d{1}'), 
              pt_stage = factor(pt_stage, c('T1', 'T2', 'T3', 'T4')), 
              pt_stage = droplevels(pt_stage))
  
# Annotation ------
  
  insert_msg('Annotation')
  
  gse116918$annotation <- fData(gse116918$raw[[1]]) %>% 
    as_tibble
  
  gse116918$annotation <- gse116918$annotation %>% 
    transmute(probe_id = ID, 
              entrez_id = stri_extract(`Entrez Gene`, regex = '\\d+')) %>% 
    annotate_raw_symbol
  
# expression --------
  
  insert_msg('Expression')
  
  ## duplicated probes are aggregated by arithmetic means of log2 expression
  ## signals
  
  gse116918$expression <- exprs(gse116918$raw[[1]]) %>% 
    integrate_expression(gse116918$annotation) %>% 
    left_join(gse116918$clinic[c('sample_id', 'patient_id', 'tissue_type')], ., 
              by = 'sample_id')
  
  gse116918$annotation <- gse116918$annotation %>% 
    filter(!duplicated(gene_symbol), 
           !duplicated(entrez_id))
  
# Caching the results ------
  
  insert_msg('Caching the results')
  
  gse116918 <- gse116918[c("clinic", "annotation", "expression")]
  
  save(gse116918, file = './data/gse116918.RData')
  
# END --------
  
  insert_tail()
  