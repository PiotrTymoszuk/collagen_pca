# TCGA import and cleaning toolbox

# libraries ----

  library(plyr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(readr)
  library(stringi)
  
# manifest retrieval ----

  read_manifest <- function(path) {
    
    ## reads the file manifest imported from the TCGA data portal
    
    return(read_tsv(path))
    
  }
  
  add_file_paths <- function(manifest_table, path_prefix) {
    
    ## adds a fixed prefix to the 'filename' column of the manifest
    
    out_manifest <- manifest_table %>% 
      mutate(filename = paste(path_prefix, filename, sep = '/'))
    
    return(out_manifest)
    
  }

# expression data import functions ----
  
  read_expression_file <- function(path, 
                                   file_id, 
                                   transpose = TRUE, 
                                   suppress_message = TRUE, ...) {
    
    ## reads expression data from a file defined by the given path
    ## If transpose, a tibble is returned with column names corresponding to transcript references 
    ## and an extra column id for file id is added
    ## If not transpose, a vector of file id and expression values named with transcript ids is returned
    ## the suppress_message_argument specifies, if the parsing info should be printed
    ## ... are additional arguments passed to the reader function
    
    if(suppress_message) {
      
      suppressMessages(expression_table <- read_tsv(path, 
                                                    col_names = FALSE, ...))
      
    } else {
      
      expression_table <- read_tsv(path, col_names = FAlse, ...)
      
    }
    
    names(expression_table) <- c('transcript_id', file_id) 
   
    if(transpose) {
      
      new_colnames <- expression_table$transcript_id
      
      expression_table <- expression_table[[file_id]] %>% 
        map_dfc(as.numeric) %>% 
        set_names(new_colnames)
      
      expression_table$file_id <- file_id
      
    }
    
    return(expression_table)
    
  }
  
  read_expression_data <- function(expression_data_folder = 'expression', 
                                   manifest_file = 'MANIFEST.txt', 
                                   tibble_output = FALSE, 
                                   parallel = FALSE, ...) {
    
    ## reads tcga rna seq expression data from a given folder containing the manifest file and subfolders with expression .tar files
    ## ... specifies additional arguments passed to the read_expression_file function
    
    print('Reading expression file manifest')
    
    expr_manifest <- read_manifest(print_path(paste(expression_data_folder, 
                                                    manifest_file, 
                                                    sep = '/'))) %>% 
      filter(id != '\\N') %>% ### filtering otu additional files not assigned to an assay
      add_file_paths(., expression_data_folder)
    
    print('Reading expression files')
    
    if(parallel) {
      
      start = Sys.time()
      
      plan('multiprocess')
      
      expression_data <- future_map2(expr_manifest$filename, 
                                     expr_manifest$id, 
                                     read_expression_file)
      
      print('Elapsed time:')
      print(Sys.time() - start)
      
    } else {
      
      start = Sys.time()
      
      expression_data <- map2(expr_manifest$filename, 
                              expr_manifest$id, 
                              read_expression_file)
      
      print('Elapsed time:')
      print(Sys.time() - start)
      
    }
    
    
    if(tibble_output) {
      
      print('Creating expression table')
      
      start = Sys.time()
      
      expression_data <- fast_rbind(expression_data) %>% 
        map_dfc(parse_guess) ### probably a faster alternative to the option retain_class in the fast_rbind function
      
      print('Elapsed time:')
      print(Sys.time() - start)  
      
    }
    
    return(expression_data)
    
  }
  
# assignment data import and append functions ----
  
  read_assignment <- function(assignment_data_folder = 'sample_assignment', 
                              expand_barcodes = TRUE, 
                              barcode_variable = 'entity_submitter_id') {
    
    ## reads assignment data (rna_seq file id with case_id and additional information)
    ## from a metadata json file, which can be downloaded from GDC parallel with expression data files
    ## it's assumed (and tested by hand) that the row order in the case information data set
    ## and the file data table match 100%
    ## optionally: expands the barcode in the output assignment file
    
    require(jsonlite)
    
    assignment_file <- list.files(assignment_data_folder, 
                                  pattern = '*.\\.json')
    
    assignment_file <- fromJSON(paste('sample_assignment', 
                                      assignment_file, 
                                      sep = '/'))
    
    case_information <- assignment_file$associated_entities %>% 
      reduce(rbind)
    
    file_data <- assignment_file[c('file_id', 'file_name', 'submitter_id')]
    
    assignment_file <- cbind(case_information, file_data) %>% 
      as_tibble
    
    if(expand_barcodes) {
      
      assignment_file <- expand_barcode(assignment_file, 
                                        barcode_variable = barcode_variable)
      
    }
    
    return(assignment_file)
    
  }
  
  add_assignment <- function(assignment_table, 
                             expression_table, 
                             index_column = 'file_id', 
                             cols_to_select = c('case_id', 
                                                'file_id', 
                                                'entity_submitter_id', 
                                                'sample_type'), 
                             join_mode = 'left') {
    
    ## The function merges selected columns (cols_to_select) from the given assignment table with the given expression table
    ## the merging key is defined by the index column. Joining modes available (LHS = assignment_table): left, right or inner
    ## for more info, see the dplyr package
    
    assign_new_tbl <- assignment_table %>%
      select(all_of(cols_to_select))
    
    join_fun <- list(left = left_join, 
                     right = right_join, 
                     inner = inner_join)
    
    return(join_fun[[join_mode]](assign_new_tbl, 
                                 expression_table, 
                                 by = index_column))
    
  }
  
# clinical data import functions ----
  
  read_clinical_file <- function(path, file_id, ...) {
    
    ## reads and clears a given clinical data file with a specified path and file_id
    ## ... are additional arguments passed to the xml parser
    
    require(XML)
    
    prim_clinical_file <- xmlParse(path)
    
    out_clinical_file <- xmlRoot(prim_clinical_file)['patient'] %>% 
      xmlToDataFrame %>% 
      as_tibble %>% 
      mutate(file_id = file_id) %>% 
      mutate(case_id = tolower(bcr_patient_uuid)) ### adding case_id in the same format as in the expression
    ### assignment file
    
    ## obtaining gleason scoring and staging data
    
    stage_info <- xmlToList(prim_clinical_file)
    
    stage_info <- stage_info$patient$stage_event
    
    psa <- tryCatch(stage_info$psa$psa_value$text, 
                    error = function(e) return(NA))
    
    gleason_score <- tryCatch(stage_info$gleason_grading$gleason_score$text, 
                              error = function(e) return(NA))
    
    primary_gleason <- tryCatch(stage_info$gleason_grading$primary_pattern$text, 
                                error = function(e) return(NA))
    
    secondary_gleason <- tryCatch(stage_info$gleason_grading$secondary_pattern$text, 
                                  error = function(e) return(NA))
    
    staging <- data.frame(psa = psa, 
                          gleason_score = gleason_score, 
                          primary_gleason = primary_gleason, 
                          secondary_gleason = secondary_gleason)
    
    return(cbind(out_clinical_file, staging) %>% as_tibble)
    
  }
  
  read_clinical_data <- function(clinical_data_folder = 'clinical', 
                                 manifest_file = 'MANIFEST.txt', 
                                 tibble_output = FALSE, 
                                 parallel = FALSE, 
                                 remove_omf = TRUE, 
                                 guess_format = TRUE, ...) {
    
    ## reads tcga clinical data from a given folder containing the manifest file and subfolders with clinical xml files
    ## ... specifies additional arguments passed to the read_clinical_file function
    
    print('Reading clinical file manifest')
    
    clinical_manifest <- 
      read_manifest(print_path(paste(clinical_data_folder, 
                                     manifest_file, sep = '/'))) %>% 
      filter(id != '\\N') %>% ### filtering additional files not assigned to the clinical information
      add_file_paths(., clinical_data_folder)
    
    ## clearing the manifest from txt files and xml files containing the omf (other malignancy form) data
    
    clinical_manifest <- clinical_manifest %>% 
      filter(!stri_detect(filename, fixed = '.txt'))
    
    if(remove_omf) {
      
      clinical_manifest <- clinical_manifest %>% 
        filter(!stri_detect(filename, fixed = 'omf'))
      
    }
    
    print('Reading clinical files')
    
    if(parallel) {
      
      start = Sys.time()
      
      plan('multiprocess')
      
      clinical_data <- future_map2(clinical_manifest$filename, 
                                   clinical_manifest$id, 
                                   read_clinical_file, ...)
      
      print('Elapsed time:')
      print(Sys.time() - start)
      
    } else {
      
      start = Sys.time()
      
      clinical_data <- map2(clinical_manifest$filename, 
                            clinical_manifest$id, 
                            read_clinical_file, ...)
      
      print('Elapsed time:')
      print(Sys.time() - start)
      
    }
    
    if(tibble_output) {
      
      print('Creating clinical data table')
      
      start = Sys.time()
      
      clinical_data <- do.call('rbind', clinical_data)
      
      print('Elapsed time:')
      print(Sys.time() - start)  
      
    } 
    
    if(guess_format) {
      
      clinical_data <- clinical_data %>% 
        map_dfc(parse_guess) ### trying to find the right format for each column
      
    }
    
    return(clinical_data)
    
  }

# experiment import and working tibble generation functions ----
  
  read_experiment <- function(assignment_data_folder = 'sample_assignment', 
                              expression_data_folder = 'expression', 
                              clinical_data_folder = 'clinical', 
                              manifest_file = 'MANIFEST.txt', 
                              parallel = TRUE) {
    
    ## A combi function reading the complete set of a TCGA expression experiment data.
    ## For info for the particular assignment, expression and clinical reading functions, please refere to
    ## the headers of the read_assignment, read_expression_data and read_clinical_data
    ## The function returns a list with three components: $assignment containing the data to match expression
    ## with clinical information, $expression containing a tibble with colnames corresponding to transcript IDs, file_id, 
    ## case_id, entity_submitter_id (TCGA barcode) and sample type (control, tumor, non-malignant) columns 
    ## and $clinical containing all available clinical information. 
    ## The indexing variable for relations between the clinical and assignment table is case_id corresponding to case_UUID.
    ## the indexing variable for relations between the expression and assignment table if file_id corresponding to file_UUID
    
    ### reading assignment data
    
    assignment_data <- 
      read_assignment(assignment_data_folder = assignment_data_folder, 
                      expand_barcodes = TRUE, 
                      barcode_variable = 'entity_submitter_id')
    
    ### reading expression data and appending it with 
    
    expression_data <- 
      read_expression_data(expression_data_folder = expression_data_folder, 
                           parallel = parallel, 
                           tibble_output = TRUE, 
                           manifest_file = manifest_file)
    
    print('Appending expression table with assignment data')
    start = Sys.time()
    
    expression_data <- 
      add_assignment(assignment_data, expression_data, 
                     index_column = 'file_id', 
                     cols_to_select = c('case_id', 
                                        'file_id', 
                                        'entity_submitter_id', 
                                        'sample_type'), 
                     join_mode = 'left')
    
    print('Elapsed time:')
    print(Sys.time() - start)
    
    #### reading clinical data
    
    clinical_data <- read_clinical_data(clinical_data_folder = clinical_data_folder, 
                                        manifest_file = manifest_file, 
                                        tibble_output = TRUE, 
                                        parallel = TRUE, 
                                        remove_omf = TRUE)
    
    
    return(list(assignment = assignment_data, 
                expression = expression_data, 
                clinical = clinical_data))
    
  }
  
  create_working_table <- function(experiment_object, 
                                   gene_expression_cols, 
                                   new_gene_expression_names = NULL) {
    
    ## a function creating an easy-to-handle tibble containing the entire clinical information, sample type data
    ## and expression for genes of interest (gene_expression_cols takes a vector with valid gene_id.version)
    # expression columns can be renamed for convenience with new_gene_expression_names vector
    
    new_assignment <- experiment_object$assignment %>% 
      select(case_id, entity_submitter_id, sample_type, file_id)
    
    new_clinics <- experiment_object$clinical %>% 
      select(- file_id) %>% 
      left_join(new_assignment, ., by = 'case_id')
    
    new_expression <- 
      experiment_object$expression[c('case_id', 'file_id', gene_expression_cols)]
    
    if(!is.null(new_gene_expression_names)) {
      
      new_expression <- 
        set_names(new_expression, 
                  c('case_id', 'file_id', new_gene_expression_names))
      
    }
    
    working_table <- left_join(new_clinics, 
                               new_expression, 
                               by = 'file_id') %>% 
      mutate(case_id = case_id.x) %>% 
      select(- file_id, - case_id.y, - case_id.x)
      
    
    return(working_table)
    
  }
  
  get_tumor_expression <- function(inp_table, 
                                   sample_type_var = 'sample_type', 
                                   tumor_value = 'primary_tumor', 
                                   duplicate_handling_fun = select_recent_records, ...) {
    
    ## a function extracting tumor-only expression from a given expression table
    ## applies a duplicate handling function to the data (... for its additional arguments), if provided
    
    tumor_only <- inp_table %>% 
      filter(.[[sample_type_var]] == tumor_value)
    
    if(!is.null(duplicate_handling_fun)) {
      
      tumor_only <- duplicate_handling_fun(tumor_only, ...)
      
    }
    
    return(tumor_only)
    
  }
  
# expression data duplicate handling ---- 
  
  select_recent_records <- function(inp_table, index_variable = 'case_id', barcode_variable = 'entity_submitter_id') {
    
    ## The function finds duplicates in the input expression table based on the index_variable,
    ## expands the barcodes and selects those samples which display a higher barcode-encoded sample_vial number
    ## as recommended by the issuing repository (Broad Institute)
    
    select_actual_record <- function(index_var_value) {
      
      return(duplicated_records %>% 
               filter(.[[index_variable]] == index_var_value) %>% filter(sample_vial == max(sample_vial)))
      
    }
    
    duplicated_records <- find_duplicates(inp_table, 
                                          index_variable = index_variable) %>% 
      select(- sample_type) %>% 
      expand_barcode(barcode_variable = barcode_variable)
    
    unique_records <- inp_table %>% 
      filter(!.[[index_variable]] %in% duplicated_records[[index_variable]])
    
    recent_records <- unique(duplicated_records[[index_variable]]) %>% 
      map_dfr(select_actual_record) %>% 
      select( - project, - tss, - participant, 
              - sample_vial, - portion_analyte, - plate, - center)
    
    return(rbind(unique_records, recent_records))
    
  }
  
# patient stratification ----
  
  stratify_samples <- function(inp_table, 
                               stratification_assignment, 
                               index_variable = 'case_id') {
    
    ## appends a given table with stratification data by merging by a given index variable
    
    return(left_join(inp_table, stratification_assignment, by = index_variable))
    
  }
  
  stratify_samples_expression <- function(inp_table, 
                                          expression_variable, 
                                          cutoff, 
                                          new_var_name = 'strata', 
                                          new_levels = c('low', 'high')) {
    
    ## stratifies samples in the given table in hi and low expressors by expression variable and cutoff
    ## The new stratification variable can be named with new_var_name
    
    strat_table <- inp_table
    
    strat_table[[new_var_name]] <- 
      ifelse(strat_table[[expression_variable]] < cutoff, 
             new_levels[1], 
             new_levels[2])
    
    return(strat_table)
    
  }
  
# varia ----

  print_path <- function(sub_path, parent_dir = getwd()) {
    
    ## returns a read to use path to a given file or folder
    
    return(paste(parent_dir, sub_path, sep = '/'))
    
  }
  
  convert_class <- function(vector, output_class) {
    
    function_list <- list(as.numeric, as.character) %>% 
      set_names('numeric', 'character')
    
    return(function_list[[output_class]](vector))
    
  }
  
  fast_rbind <- function(inp_list, 
                         retain_class = FALSE, 
                         parallel = FALSE) {
    
    new_colnames <- names(inp_list[[1]])
    
    output_tibble <- matrix(unlist(inp_list), 
                            nrow = length(inp_list), 
                            byrow = TRUE) %>% 
      as_tibble %>% 
      set_names(new_colnames)
    
    if(retain_class){
      
      class_vector <- map_dfc(inp_list[[1]], class)
      
      if(parallel) {
        
        plan('multiprocess')
        
        output_tibble <- future_map2_dfc(output_tibble, 
                                         class_vector, 
                                         convert_class)
        
      } else {
        
        output_tibble <- map2_dfc(output_tibble, 
                                  class_vector, 
                                  convert_class)
        
      }
      
    }
    
    return(output_tibble)
    
  }
  
  find_duplicates <- function(inp_table, index_variable) {
    
    ## a function seeking for records containing records with the duplicated index variable
    
    duplicate_table <- inp_table %>% 
      filter(duplicated(.[[index_variable]]))
    
    duplicate_table <- inp_table %>% 
      filter(.[[index_variable]] %in% duplicate_table[[index_variable]])
    
    return(duplicate_table %>% arrange(.[[index_variable]]))
    
  }
  
  expand_barcode <- function(inp_table, 
                             barcode_variable = 'entity_submitter_id') {
    
    ## expands the TCGA barcodes stored in the barcode_variable column of the given input table
    ## returns a copy of the inp_table appended with the crunched barcode information
    ## for more information on barcode reading, visit: 
    ## https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
    
    crunched_barcodes <- inp_table[[barcode_variable]] %>% 
      stri_split_fixed('-', simplify = TRUE) %>% 
      as_tibble %>% 
      set_names(c('project', 'tss', 'participant', 'sample_vial', 'portion_analyte', 'plate', 'center'))
    
    crunched_barcodes <- crunched_barcodes %>% 
      mutate(sample_type = as.numeric(stri_extract(sample_vial, regex = '^\\d\\d'))) %>% 
      mutate(sample_type = ifelse(sample_type %in% 0:9, 'primary_tumor', ### information on sample type added
                                  ifelse(sample_type %in% 10:19, 'normal_tissue', 'control_sample')))
    
    return(cbind(inp_table, crunched_barcodes) %>% as_tibble)
    
  }
  
  find_gene_version <- function(experiment_object, gene_id, max_version = 100) {
    
    ## a function checking which vesion of the given gene_id is present 
    ## in the expression tibble of the given experiment file
    ## returns a gene_id.version
    
    version_vector <- paste(gene_id, 1:max_version, sep = '.')
    
    for(i in version_vector) {
      
      if(i %in% names(experiment_object$expression)) {
        
        return(i)
        
      }
      
    }
    
  }
  
  log2_1 <- function(value) {
    
    ## calculates log2(value + 1)
    
    return(log2(value + 1))
    
  }
  
# benchmarking ----
  
  benchmark_exprs <- function(expression) {
    
    ## benchmarks a given expression (unquoted)
    
    start_time = Sys.time()
    
    eval(expression)
    
    end_time = Sys.time()
    
    print(end_time - start_time)
    
    return(end_time - start_time)
    
  }