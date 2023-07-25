# Tools for reading GEO microarray experiments from series matrix files

# libraries ----

  library(plyr)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(readr)
  library(stringi)

  source('./tools/sys_tools.R')

# functions for reading experiment design data (series matrix) ----

  get_sample_info <- function(ser_matrix_path, output_file = 'sample_info.txt', filter_features = NULL) {
    
    ## extracts sample information form a given series_matrix file and saves it as
    ## a new Tab-delimited file. If provided, filter features enable selecting specific sample
    ## characteristics. The '!Series_geo_accession' feature is always selected
    
    mtx_file <- file(ser_matrix_path, 'r')
    
    mtx <- readLines(mtx_file)
    
    close(mtx_file)
    
    
    if(!is.null(filter_features)) {
      
      sample_info <- ''
      
      for(feature in unique(c('!Sample_geo_accession', filter_features))) {
        
        sample_info <- c(sample_info, mtx[startsWith(mtx, feature)])
        
      }
      
    } else {
      
      sample_info <- mtx[startsWith(mtx, '!Sample')]
      
    }
    
    write(sample_info, output_file)
    
  }
  
  read_sample_info <- function(ser_matrix_path, filter_features = NULL) {
    
    ## reads sample information from a given series matrix file
    
    get_sample_info(ser_matrix_path = ser_matrix_path, 
                    output_file = '_sample.tmp', 
                    filter_features = filter_features)
    
    sample_info <- read_tsv('_sample.tmp', col_names = F)
    
    col_names <- stri_replace(sample_info[[1]], fixed = '!', replacement = '') %>% unlist
    
    sample_info[[1]] <- NULL
    
    sample_info <- t(sample_info)
    
    colnames(sample_info) <- col_names
    
    sample_info <- sample_info %>% 
      as_tibble %>% 
      mutate(ID_REF = Sample_geo_accession)

    return(sample_info)
    
  }
  
# functions for reading expression data -----
  
  get_expression_data <- function(ser_matrix_path, output_file = 'exprs.txt') {
    
    ## extracts sample information form a given series_matrix file and saves it as
    ## a new Tab-delimited file
    
    mtx_file <- file(ser_matrix_path, 'r')
    
    mtx <- readLines(mtx_file)
    
    close(mtx_file)
    
    sample_info <- mtx[!startsWith(mtx, '!Sample') & 
                         !startsWith(mtx, '!Series') &
                         !startsWith(mtx, '!series')]
    
    write(sample_info, output_file)
    
  }
  
  read_exprs <- function(ser_matrix_path) {
    
    ## reads sample information from a given series matrix file
    
    get_expression_data(ser_matrix_path = ser_matrix_path, 
                        output_file = '_exprs.tmp')
    
    exprs <- read_tsv('_exprs.tmp', col_names = F)
    
    col_names <- stri_replace(exprs[[1]], fixed = '!', replacement = '') %>% unlist
    
    exprs[[1]] <- NULL
    
    exprs <- t(exprs)
    
    colnames(exprs) <- col_names
    
    exprs <- exprs %>% 
      as_tibble
    
    exprs_numbers <- exprs[names(exprs) != 'ID_REF'] %>% 
      map_dfc(as.numeric)
    
    exprs <- cbind(ID_REF = exprs[['ID_REF']], exprs_numbers)
    
    return(exprs %>% as_tibble)
    
  }
  
# functions for reading annotation -----
  
  get_annotation <- function(ann_path, output_file = 'ann.txt') {
    
    ## extracts annotation form a given series_matrix file and saves it as
    ## a new Tab-delimited file.
    
    mtx_file <- file(ann_path, 'r')
    
    mtx <- readLines(mtx_file)
    
    close(mtx_file)

    ann_file <- mtx[!startsWith(mtx, '#')]

    write(ann_file, output_file)
    
  }
  
  read_annotation <- function(ann_path, filter_features = NULL, cols_to_clear = NULL) {
    
    ## reads probe annotation. Filter features eneble selecting columns of interest
    ## The 'ID column is always selected. Enables clearing columns containing the ' /// '
    ## seperator.
    
    get_annotation(ann_path = ann_path, 
                   output_file = '_ann.temp')
    
    ann <- read_tsv('_ann.temp')

    if(!is.null(filter_features)) {
      
      ann <- ann[, unique(c('ID', filter_features))]
      
    }
    
    if(!is.null(cols_to_clear)) {
      
      ann <- clear_ann(ann, cols_to_clear = cols_to_clear)
      
    }
    
    return(ann)
    
  }
  
  clear_ann <- function(ann_table, cols_to_clear, separator = ' /// |,| // ') {
    
    ## clears an annotation table. The specified columns containing
    ## the ' /// ' separator and more than 1 ID will be splitted
    ## into lists
    
    cleared_table <- ann_table
    
    for (col in cols_to_clear) {
      
      cleared_table[[col]] <- cleared_table[[col]] %>% 
        stri_split_regex(pattern = separator)
      
    }
    
    return(cleared_table)
    
  }
  
# functions for reading an experiment. The experiment list consists of the following elements:
  # 1) sample information as tibble
  # 2) tibble with expression data
  # 3) tibble with annotation information
  
  read_exprmt <- function(ser_matrix_path = NULL, ann_path = NULL, exprmt_folder = NULL, 
                          sample_filter = NULL, ann_filter = NULL, ann_clear = NULL) {
    
    ## reads an experiment. If experiment folder is provided, tries to find 
    ## the series matrix and annotation on its own
    
    start_time <- Sys.time()
    
    curr_wd = getwd()
    
    if(!is.null(exprmt_folder)) {

      ann_path <- list.files(path = exprmt_folder, pattern = 'GPL.*')
      ser_matrix_path <- list.files(path = exprmt_folder, pattern = 'GSE.*')
      
      setwd(exprmt_folder)
      
    }
    
    message(rep('=', 60))
    message(paste('Parsing series matrix file:', ser_matrix_path))
    
    sample_info <- read_sample_info(ser_matrix_path = ser_matrix_path, 
                                    filter_features = sample_filter)
    
    exprs <- read_exprs(ser_matrix_path = ser_matrix_path)
    
    message(paste('Parsing annotation file:', ann_path))
    
    ann <- read_annotation(ann_path = ann_path, 
                           filter_features = ann_filter, 
                           cols_to_clear = ann_clear)
    
    setwd(curr_wd)
    
    stop_time <- Sys.time()
    
    message(paste('Time elapsed:', stop_time - start_time))
    
    return(list(sample = sample_info, 
                exprs = exprs, 
                annotation = ann))
    
  }
  
# other useful functions ----
  
  merge_sample_info <- function(geo_exprmt_list, id_var = 'ID_REF') {
    
    ## merges sample info with expression data within the experiment
    
    start_time <- Sys.time()
    
    message(rep('=', 60))
    message('Merging sample info with expression data')
    
    merged_experiment <- geo_exprmt_list
    
    merged_experiment$exprs <- left_join(merged_experiment$sample, merged_experiment$exprs, by = id_var)
    
    stop_time <- Sys.time()
    
    message(paste('Time elapsed:', stop_time - start_time))
    
    return(merged_experiment)
    
  }
  
  annot_reverter <- function(ann_tbl) {
    
    ## reverses the annotation table by listing all gene IDs with associated probe IDs
    
    rev_tbl <- list()
    
    for(i in 1:nrow(ann_tbl)) {
      
      for(j in unlist(ann_tbl[i, 'gene_ID'])) {
        
        rev_tbl <- c(rev_tbl, list(data.frame(gene_ID = j, 
                                              ID = ann_tbl[i, 'ID'])))
        
      }
      
    }
    
    return(rev_tbl %>% 
             do.call('rbind', .))
    
  }
  
# probe and expression extracting functions -----
  
  find_probe_id <- function(geo_exprmt_list, gene_id, id_var) {
    
    ## finds probe_ID's corresponding to a given gene_ID
    ## due to list format id_variable, a loop construct is needed to extract IDs
    
    probe_number <- nrow(geo_exprmt_list$annotation)
    probe_table <- tibble()
    
    for(i in 1:probe_number){
      
      if(all(is.na(geo_exprmt_list$annotation[[id_var]][[i]]))){
        
        next()
        
      }
      
      if(any(geo_exprmt_list$annotation[[id_var]][[i]] == gene_id)){
        
        probe_table <- rbind(probe_table, cbind(geo_exprmt_list$annotation[i, ], querry_id = gene_id))

      }
      
    }
    
    return(probe_table %>% as_tibble)
    
  }
  
  extract_probe_id <- function(geo_exprmt_list, gene_ids, id_var = 'Entrez_Gene_ID', .parallel = F) {
    
    ## extracts probe IDs or fragments of the annotation table corresponding to the given gene ID vector
   
    start_time <- Sys.time()
    
    gene_list <- gene_ids[!is.na(gene_ids)]
    probe_number <- nrow(geo_exprmt_list$annotation)
    
    message(rep('=', 60))
    message('Extracting probe IDs for n = ', length(gene_list), ' genes of interest')
    message('Total number of probes to scan: ', probe_number)
    
    if(.parallel) {
      
      plan('multiprocess')
      
      probe_table <- gene_list %>% 
        future_map(function(x) find_probe_id(geo_exprmt_list = geo_exprmt_list, 
                                             gene_id = x, 
                                             id_var = id_var))
      
      plan('sequential')
      
    } else {
      
      probe_table <- gene_list %>% 
        map(function(x) find_probe_id(geo_exprmt_list = geo_exprmt_list, 
                                      gene_id = x, 
                                      id_var = id_var))
      
    }
    
    
    
    message('Time elapsed: ', Sys.time() - start_time)
    
    return(probe_table %>% 
             do.call('rbind', .))
    
  }
  
  extract_exprs <- function(geo_exprmt_list, gene_ids = NULL, 
                            probe_vec = NULL, id_var = 'Entrez_Gene_ID', keep_sample_info = T) {
    
    ## extracts expression data for given gene IDs
    
    start_time = Sys.time()
    
    if(is.null(probe_vec)) {
      
      probe_vec <- extract_probe_id(geo_exprmt_list = geo_exprmt_list, 
                                    gene_ids = gene_ids, 
                                    id_var = id_var)$ID
      
    }
    
    message(rep())
    message('Extracting expression information for n = ', length(probe_vec), ' genes of interest')
    
    if(!keep_sample_info) {
      
      message('Time elapsed: ', Sys.time() - start_time)
      message(rep('=', 60))
      
      return(geo_exprmt_list$exprs[, names(geo_exprmt_list$exprs) %in% probe_vec])
      
    } else {
      
      sample_cols <- names(geo_exprmt_list$sample)
      
      message('Time elapsed: ', Sys.time() - start_time)
      message(rep('=', 60))
      
      return(geo_exprmt_list$exprs[, names(geo_exprmt_list$exprs) %in% c(sample_cols, probe_vec)])
      
    }
    
  }
  
# Jetset (http://www.cbs.dtu.dk/biotools/jetset/) probeset selection -----
  
  select_jetset_best <- function(geo_exprmt_list, jetset_score_tbl, rename = NULL) {
    
    ## selects the best probeset corresponding to each a gene using the Jetset scoring information
    ## optionally, the probesets can be renamed as either EntrezID or symbol. In this case
    ## the expression set is restricted to the probesets with an identifier assignment
    
    if(is.null(rename)) {
      
      jetset_id <- jetset_score_tbl %>% 
        filter(best) 
      
      naming_vec <- jetset_id$probeset %>% 
        set_names(jetset_id$probeset)
      
    } else if(rename == 'symbol') {
      
      jetset_id <- jetset_score_tbl %>% 
        filter(best, 
               !is.na(symbol), 
               symbol != '--')
      
      naming_vec <- jetset_id$symbol %>% 
        set_names(jetset_id$probeset)
      
      
    } else {
      
      jetset_id <- jetset_score_tbl %>% 
        filter(best, 
               !is.na(EntrezID), 
               EntrezID != '--')
      
      naming_vec <- jetset_id$EntrezID %>% 
        set_names(jetset_id$probeset)
      
    }
    
    exprs_subset <- geo_exprmt_list$exprs %>% 
      select(ID_REF, 
             all_of(jetset_id$probeset))
    
    exprs_probes <- names(exprs_subset)[names(exprs_subset) != 'ID_REF']
    
    exprs_subset <- exprs_subset %>% 
      set_names(c('ID_REF', 
                  naming_vec[exprs_probes]))
    
    return(list(sample = geo_exprmt_list$sample, 
                exprs = exprs_subset, 
                annotation = jetset_id))
    
    
  }
  
# END ----