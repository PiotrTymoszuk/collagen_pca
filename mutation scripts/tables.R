# Computes and caches mutation indexes, counts and frequencies

  insert_head()
  
# container ------
  
  mut_tables <- list()
  
# mutation indexes ------
  
  insert_msg('Mutation indexes')
  
  ## raw matrix
  
  mut_tables$mtx <- 
    tab_mutations(xena_db = mut_globals$mutect_tbl, 
                  as_matrix = TRUE, 
                  sparse = FALSE, 
                  .parallel = TRUE)
  
  ## a data frame
  
  mut_tables$mut_tbl <- mut_tables$mtx %>% 
    as.matrix %>% 
    as.data.frame %>% 
    rownames_to_column('sample_id') 
  
# mutation counts -------
  
  insert_msg('Mutation counts per participants')
  
  mut_tables$count_tbl <- mut_tables$mtx %>% 
    rowSums %>% 
    compress(names_to = 'sample_id', 
             values_to = 'mut_count')
  
# CNV indexes --------
  
  insert_msg('CNV indexes')
  
  mut_tables$cnv_tbl <- mut_globals$gistic_tbl %>% 
    select(-ensembl_id, - `Gene Symbol`) %>% 
    column_to_rownames('gene_symbol') %>%
    t %>% 
    as.data.frame %>% 
    rownames_to_column('sample_id') %>% 
    as_tibble
  
  ## tables for amplifications and deletions
  
  mut_tables[c('ampl_tbl', 'del_tbl')] <- mut_tables[c('cnv_tbl', 'cnv_tbl')]
  
  mut_tables$ampl_tbl[mut_globals$gistic_tbl$gene_symbol] <- 
    mut_tables$ampl_tbl[mut_globals$gistic_tbl$gene_symbol] %>% 
    map_dfc(~ifelse(.x > 0, 'amplified', 'non-amplified')) %>% 
    map_dfc(factor, c('non-amplified', 'amplified'))
  
  mut_tables$del_tbl[mut_globals$gistic_tbl$gene_symbol] <- 
    mut_tables$del_tbl[mut_globals$gistic_tbl$gene_symbol] %>% 
    map_dfc(~ifelse(.x < 0, 'deleted', 'non-deleted')) %>% 
    map_dfc(factor, c('non-deleted', 'deleted'))
  
# CNV counts -------
  
  insert_msg('CNV counts')
  
  mut_tables$cnv_mtx <- mut_tables$cnv_tbl %>% 
    column_to_rownames('sample_id') %>% 
    as.matrix
   
  mut_tables$cnv_counts <- 
    list(cnv_count = mut_tables$cnv_mtx != 0, 
         ampl_count = mut_tables$cnv_mtx > 0, 
         del_count = mut_tables$cnv_mtx < 0) %>% 
    map(rowSums) %>% 
    map2(., names(.), 
         ~compress(.x, 
                   names_to = 'sample_id', 
                   values_to = .y)) %>% 
    reduce(left_join, by = 'sample_id')
  
  ## appending the general count table
  
  mut_tables$count_tbl <- 
    full_join(mut_tables$count_tbl, 
              mut_tables$cnv_counts, 
              by = 'sample_id')
  
# appending with the cluster assignment -----
  
  insert_msg('Appending with the cluster assignment')
  
  mut_tables[c("mut_tbl", "count_tbl", "cnv_tbl", "ampl_tbl", "del_tbl")] <- 
    mut_tables[c("mut_tbl", "count_tbl", "cnv_tbl", "ampl_tbl", "del_tbl")] %>% 
    map(mutate, 
        patient_id = stri_extract(sample_id, 
                                  regex = '\\w{4}-\\w{2}-\\w{4}'), 
        patient_id = stri_extract(patient_id, regex = '\\w{4}$')) %>% 
    map(inner_join, 
        mut_globals$assignment, 
        by = 'patient_id') %>% 
    map(as_tibble) %>% 
    map(relocate, clust_id) %>% 
    map(relocate, patient_id)
  
  ## and restricting to the primary tumor samples
  ## sample_id ending with the '-01A' string

  mut_tables[c("mut_tbl", "count_tbl", "cnv_tbl", "ampl_tbl", "del_tbl")] <- 
    mut_tables[c("mut_tbl", "count_tbl", "cnv_tbl", "ampl_tbl", "del_tbl")] %>% 
    map(filter, stri_detect(sample_id, regex = '-01A$'))
    
# mutation and copy number frequencies, per gene ------
  
  insert_msg('Mutation an CNV frequencies per gene')
  
  ## restricting to the samples with the known cluster assignment
  
  mut_tables$frequency_tbl[c('mut', 'cnv', 'ampl', 'del')] <- 
    mut_tables[c("mtx", "cnv_mtx", "cnv_mtx", "cnv_mtx")] %>% 
    map(as.matrix) %>% 
    map2(., 
         map(mut_tables[c("mut_tbl", "cnv_tbl", "ampl_tbl", "del_tbl")], 
             ~.x$sample_id), 
         ~.x[.y, ])
    
  mut_tables$frequency_tbl$cnv <- mut_tables$frequency_tbl$cnv != 0
  mut_tables$frequency_tbl$ampl <- mut_tables$frequency_tbl$ampl > 0
  mut_tables$frequency_tbl$del <- mut_tables$frequency_tbl$del < 0
  
  mut_tables$frequency_tbl <- mut_tables$frequency_tbl %>% 
    map(colSums) %>% 
    map2(., paste0('n_', names(.)), 
        ~compress(.x, 
                  names_to = 'gene_symbol', 
                  values_to = .y)) %>% 
    reduce(full_join, by = 'gene_symbol')
  
  mut_tables$frequency_tbl[c('n_mut', 'n_cnv', 'n_ampl', 'n_del')] <- 
    mut_tables$frequency_tbl[c('n_mut', 'n_cnv', 'n_ampl', 'n_del')] %>% 
    map_dfc(~ifelse(is.na(.x), 0, .x))
  
  ## computing the percentages
  
  mut_tables$frequency_tbl <- 
    mut_tables$frequency_tbl %>% 
    mutate(n_mut_total = nrow(mut_tables$mut_tbl), 
           n_cnv_total = nrow(mut_tables$cnv_tbl), 
           n_ampl_total = nrow(mut_tables$ampl_tbl), 
           n_del_total = nrow(mut_tables$del_tbl)) %>% 
    mutate(perc_mut = n_mut/n_mut_total * 100, 
           perc_cnv = n_cnv/n_cnv_total * 100, 
           perc_ampl = n_ampl/n_ampl_total * 100, 
           perc_del = n_del/n_del_total * 100)
  
# candidate outliers, 6 sigma threshold ------
  
  insert_msg('Candidate outliers, 6-sigma')
  
  ## 6 sigma of the counts and CNV
  
  mut_tables$outliers <- mut_tables$count_tbl
  
  mut_tables$outliers <- 
    mut_tables$outliers[c('mut_count', 'cnv_count', 'ampl_count', 'del_count')] %>% 
    map(~scale(.x)[, 1]) %>% 
    map(set_names, mut_tables$outliers$patient_id)
  
  mut_tables$outliers <- mut_tables$outliers %>% 
    map(~.x[!is.na(.x)]) %>% 
    map(~.x[abs(.x) > 6]) %>% 
    map(names) %>% 
    reduce(union)
  
# top alterations -----
  
  insert_msg('Top most frequent genetic alterations')
  
  ## mutations: present in at least 2% of participants
  ## CNV: at least 10% of participants
  ## amplifications: at least 4% of participants
  ## deletions: at least 10% of participants
  
  mut_tables$top[c('mutations', 
                   'cnv', 
                   'amplifications', 
                   'deletions')] <- 
    c('perc_mut', 'perc_cnv', 'perc_ampl', 'perc_del') %>% 
    map2(., c(2, 10, 4, 10), 
         ~filter(mut_tables$frequency_tbl, .data[[.x]] >= .y)) %>% 
    map(~.x$gene_symbol)

# removing redundant stuff ------
  
  insert_msg('Removal of the redundant stuff')
  
  mut_tables$mtx <- NULL
  mut_tables$cnv_mtx <- NULL
  mut_tables$cnv_counts <- NULL
  
  mut_tables <- compact(mut_tables)
  
# caching the results -----
  
  insert_msg('Caching the results')
  
  save(mut_tables, file = './cache/mut_tables.RData')
  
# END -----
  
  insert_tail()