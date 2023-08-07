# This is a collection of functions used for data import

# libraries ----

  require(plyr)
  require(tidyverse)
  require(tidyr)
  require(readxl)
  library(furrr)

# accessory functions ----

  geom_mean <- function(x) {
    
    ## calculates a geometric mean for a vector or numeric table
    
    if(any(class(x) == 'data.frame')) {
      
      if(ncol(x) > 1) {
        
        geom_vec <- 1:nrow(x) %>% 
          map(function(a) unlist(x[a, ])) %>% 
          map_dbl(geom_mean)
        
        return(tibble(gmean = geom_vec))
        
      } else {
        
        return(geom_mean(x[[1]]))
        
      }
      
    }
    
    if(!is.numeric(x)) {
      
      warning('Vector has to be numeric')
      
      return(NULL)
      
    }
    
    if(is.null(x)) {
      
      warning('The input has to be positive')
      
      return(NULL)
      
    }
    
    if(any(x < 0)) {
      
      warning('The input has to be positive')
      
      return(NULL)
      
    }
    
    geom <- x %>% 
      reduce(function(a, b) a * b)
    
    return(geom^(1/length(x)))
    
  }
  
  select_gene_exprs <- function(expression_tbl, annotation_tbl, gene_id, 
                                identifyier_type = 'gene_ID') {
    
    ## selects all probes corresponding to the given gene
    
    probes_interest <- annotation_tbl %>% 
      filter(.data[[identifyier_type]] == gene_id) %>% 
      .$probe_ID
    
    if(length(probes_interest) == 0) {
      
      warning('No probes corresponding to the gene number', call. = FALSE)
      
      return(NULL)
      
    }
    
    probes_present <- probes_interest[probes_interest %in% names(expression_tbl)]
    
    if(length(probes_present) < length(probes_interest)) {
      
      warning('Some probes absent from the expression table', call. = FALSE)
      
    }
    
    return(expression_tbl[, probes_present])
    
  }

# import functions ----

  transpose_tibble <- function(input_tibble) {
  
    ## a function wihich transposes a tibble and preserving the class
    ## colnames of the input tibble will be stored in a new column sample_ID
  
    output <- data.frame(input_tibble)
    
    rownames(output) <- output[, 1]
  
    output[, 1] <- NULL
  
    output <- data.frame(t(output))
  
    output$sample_ID <- rownames(output)
  
    output <- as_tibble(output)
  
    return(output)
  
  }

  read_expression <- function(folder_name, transpose = TRUE) {
    
    ## a function which reads the expression.xlsx file in the given folder in the working directory
    ## and, on request, transposes it and supplies with a sample_ID information
    
    require(stringi)
    
    path <- paste(getwd(), '/', folder_name, '/expression.xlsx', sep = '')
    
    output <- read_excel(path)
    
    for (i in colnames(output)) {
      
      ### replacing null with numeric 0
      
      output[[i]] <- stri_replace(output[[i]], replacement = 0, fixed = 'null')
      
      output[[i]] <- parse_guess(output[[i]])
      
      }
    
    
    if(transpose) {
      
      output <- transpose_tibble(output)
      
    }
    
    return(output)
    
  }
  
  read_design <- function(folder_name, clearing_phrases = NA) {
    
    ## a function which reads the design.xlsx file from a given folder in the working directory
    ## clearing_phrases is a vector of strings to be removed from the output table
    
    require(stringi)
    
    path <- paste(getwd(), '/', folder_name, '/design.xlsx', sep = '')
    
    output <- read_excel(path)
    
    if(any(!is.na(clearing_phrases))) {
      
      for(i in colnames(output)) {
        
        for(j in clearing_phrases) {
          
          output[[i]] <- stri_replace(output[[i]], replacement = '', fixed = j)
          
          output[[i]] <- parse_guess(output[[i]])
          
        }
        
      }
      
    }
    
    return(output)
    
  }
  
  read_annotation <- function(folder_name) {
    
    ## a function which reads the annotation.xlsx file from a given folder in the working directory
    
    path <- paste(getwd(), '/', folder_name, '/annotation.xlsx', sep = '')
    
    output <- read_excel(path)
    
    if(class(output$probe_ID) == 'numeric') {
      
      output$probe_ID <- paste('X', output$probe_ID, sep = '')
      
    }
    
    return(output)
    
  }
  
  untangle_gleason_sum <- function(design_dataset) {
    
    require(stringi)
    
    output <- design_dataset
    
    if(any('gleason_sum' %in% colnames(output))) {
      
      gleason_split <- stri_split_regex(output$gleason_sum, 
                                        pattern = "\\:|\\+|\\=", 
                                        simplify = TRUE)
      
      output$gleason <- stri_replace(gleason_split[, ncol(gleason_split) - 2], 
                                     replacement = 'NA', 
                                     fixed = 'N/A') 
      output$gleason <- as.numeric(output$gleason)
      
      output$major_gleason <- gleason_split[, ncol(gleason_split) - 1]
      output$major_gleason <- as.numeric(output$major_gleason )
      
      output$minor_gleason <- gleason_split[, ncol(gleason_split)]
      output$minor_gleason <- as.numeric(output$minor_gleason)
      
    }
    
    return(output)
    
  }
  
  read_study <- function(folder_name, merge_clinical = T, ...) {
    
    ## a function which reads the files annotation.xlsx, expression.xlsx and design.xlsx from a give folder
    ## of the working direcrtory. Optionally, clearing_phrases can be passed to the design reading function.
    ## The design data set is merged with the transosed expression data set by sample_ID (by left joining, if merge_clinical is T).
    ## a pair of data sets expression and annotation (after merging with the design) or a triplet of 
    ## design, expression and annotation data sets is returned as a named list (expression and annotation)
    
    
    annotation <- read_annotation(folder_name)
    
    expression <- read_expression(folder_name, TRUE)
    
    design <- read_design(folder_name, ...)
    
    design <- untangle_gleason_sum(design)
      
    if(merge_clinical) {
      
      expression <- left_join(design, expression, by = 'sample_ID')
      
      output <- list(expression, annotation)
      
      output <- set_names(output, c('expression', 'annotation'))
      
    } else {
      
      output <- list(expression, design, annotation)
      
      output <- set_names(output, c('expression', 'design', 'annotation'))
      
    }
    
   
    return(output)
    
  }
  
# expression integrating functions -----
  
  integrate_expression <- function(expression_tbl, 
                                   annotation_tbl, 
                                   gene_identifier = 'gene_ID', 
                                   int_fun = geom_mean, 
                                   trans_fun = function(x) x, 
                                   .parallel = FALSE) {
    
    ## integrates the expression from multiple probes corresponding to the same gene
    
    ## splitting the tables into fixed design/clinical columns
    ## and the columns with expression values
    
    start_time <- Sys.time()
    
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))
    
    fix_columns <- names(expression_tbl)[!names(expression_tbl) %in% annotation_tbl$probe_ID]
    
    fix_part <- expression_tbl[, fix_columns]
    
    ## unique genes with EntrezID or other identifier present in the expression set
    ## the user-provided integration function has to return a data frame
    
    unique_genes <- annotation_tbl[[gene_identifier]] %>% 
      unique

    message(paste('Integrating expression for', 
                  length(unique_genes), 
                  'genes'))
    
    integrator <- function(expr_tbl, int_fun) {
      
      if(is.null(expr_tbl)) {
        
        return(NULL)
        
      }
      
      if(ncol(expr_tbl) == 0) {
        
        return(NULL)
        
      } else if(ncol(expr_tbl) == 1) {
        
        return(expr_tbl)
        
      } else {
        
        output <- tryCatch(int_fun(expr_tbl), 
                           error = function(e) NULL)
        
        return(output)
        
      }
      
    }
    
    if(.parallel) {
      
      plan('multisession')
      
      int_expr <- unique_genes %>% 
        future_map(select_gene_exprs, 
                   annotation_tbl = annotation_tbl, 
                   expression_tbl = expression_tbl, 
                   identifyier_type = gene_identifier) %>% 
        set_names(unique_genes)
      
      int_expr <- int_expr %>% 
        future_map(integrator, 
                   int_fun = int_fun) %>% 
        future_map(function(x) tryCatch(trans_fun(x), 
                                        error = function(e) NULL))
      
      plan('sequential')
      
    } else {
      
      int_expr <- unique_genes %>% 
        map(select_gene_exprs, 
            annotation_tbl = annotation_tbl, 
            expression_tbl = expression_tbl, 
            identifyier_type = gene_identifier) %>% 
        set_names(unique_genes)
      
      int_expr <- int_expr %>% 
        map(integrator, 
            int_fun = int_fun) %>% 
        future_map(function(x) tryCatch(trans_fun(x), 
                                        error = function(e) NULL))
      
      
    }
    
    int_expr <- int_expr %>% 
      compact

    int_expr <- int_expr %>% 
      do.call('cbind', .) %>% 
      set_names(names(int_expr))

    mod_tbl <- cbind(fix_part,
                     int_expr) %>% 
      as_tibble
    
    correct_genes <- unique_genes[unique_genes %in% names(mod_tbl)]
    
    correct_genes <- annotation_tbl %>% 
      filter(.data[[gene_identifier]] %in% correct_genes)
    
    message(paste('Elapsed:', Sys.time() - start_time))
    
    return(list(expression = mod_tbl, 
                annotation = correct_genes))
    
  }
  
# duplicate handling ------
  
  rm_duplicates <- function(expression_tbl, 
                            patient_identifier = 'patient_id', 
                            reposit_identifier = 'sample_ID') {
    
    ## screens the expression table for duplicates and chooses the sample
    ## with the latest repository id
    
    start_time <- Sys.time()
    message(paste('Screening for duplicates within', 
                  nrow(expression_tbl), 
                  'observations'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))
    
    
    clean_tbl <- expression_tbl %>% 
      dlply(patient_identifier, as_tibble) %>% 
      map(arrange, 
          desc(.data[[reposit_identifier]])) %>% 
      map_dfr(function(x) x[1, ])

    return(clean_tbl %>% 
             as_tibble)
    
    
  }
  
# END ----