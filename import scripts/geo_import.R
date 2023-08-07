# Import of the GEO datasets 


insert_head()

# tools ------

  library(tidyverse)
  library(GEOquery)

# data containers -----

  geo_globals <- list()
  geo_data <- list()
  
# fetching the GEO objects ------
  
  insert_msg('Fetching the GEO objects')
  
  geo_globals$ids <- 
    c('GSE16560', 
      'GSE70768', 
      'GSE70769', 
      'GSE116918')
  
  geo_globals$ids <- 
    set_names(geo_globals$ids, geo_globals$ids)
  
  geo_data$objects <- geo_globals$ids %>% 
    map(~getGEO(.x, destdir = './data'))
  
# Unpacking the annotation data ------
  
  insert_msg('Unpacking the annotation data')
  
  ## in case of multiple Entrez IDs per probe, the first Entrez ID is taken
  
  geo_data$annotation <- geo_data$objects %>% 
    map(~.x[[1]]) %>% 
    map(fData) %>% 
    map(as_tibble)
  
  geo_data$annotation$GSE16560 <- geo_data$annotation$GSE16560 %>% 
    transmute(probe_ID = ID, 
              gene_ID = stri_extract(GeneID, regex = '\\d+'), 
              transcript_ID = GB_ACC, 
              symbol = Symbol, 
              description = Description)
  
  geo_data$annotation[c("GSE70768", "GSE70769")] <- 
    geo_data$annotation[c("GSE70768", "GSE70769")] %>% 
    map(transmute, 
        probe_ID = ID, 
        gene_ID = stri_extract(Entrez_Gene_ID, regex = '\\d+'), 
        transcript_ID = Source_Reference_ID, 
        symbol = Symbol, 
        description = Definition)
  
  geo_data$annotation$GSE116918 <- geo_data$annotation$GSE116918 %>% 
    transmute(probe_ID = ID, 
              gene_ID = stri_extract(`Entrez Gene`, regex = '\\d+'), 
              symbol = stri_extract(`Gene Symbol`, regex = '^\\w+'), 
              description = `Gene Description`)
  
  ## restricting the annotation data only to the genes Entrez ID
  
  geo_data$annotation <- geo_data$annotation %>% 
    map(filter, !is.na(gene_ID))
  
# Expression data -------
  
  insert_msg('Expression data')
  
  geo_data$expression <- geo_data$objects %>% 
    map(~.x[[1]]) %>% 
    map(exprs) %>% 
    map(t) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_ID') %>% 
    map(as_tibble)
  
# Phenotype data ------
  
  insert_msg('Phenotype data')
  
  geo_data$clinic <- geo_data$objects %>% 
    map(~.x[[1]]) %>% 
    map(pData) %>% 
    map(as_tibble)

# Clearing the container list -----
  
  insert_msg('Cleaning the container list')
  
  geo_data <- geo_data[c("clinic", "expression", "annotation")] %>% 
    transpose
  
# END -----
  
  insert_tail()