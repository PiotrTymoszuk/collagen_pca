# This script clears the TCGA data -----

  insert_head()

# container ------

  tcga_cleared <- list()

# merging the clinical and expression data ------

  insert_msg('Merging the clinical and expression data')

  tcga_data$expression <- tcga_data$expression %>% 
    left_join(tcga_data$assignment[, c('entity_submitter_id', 'participant')] %>% 
                set_names(c('entity_submitter_id', 'patient_id')), 
              by = 'entity_submitter_id') %>% 
    left_join(tcga_data$clinical, 
              by = 'patient_id')
  
# clearing the data: uniform format for fup times, tissue etc -----
  
  insert_msg('Clearing the data')
  
  tcga_data$expression <- tcga_data$expression %>% 
    mutate(death = ifelse(vital_status == 'Alive', 0, 1), 
           vitality_fup = ifelse(death == 1, 
                                 days_to_death/30.436875, 
                                 days_to_last_followup/30.436875), 
           relapse = ifelse(biochemical_recurrence == 'NO', 0, 1), 
           relapse_fup = ifelse(relapse == 1, 
                                days_to_first_biochemical_recurrence/30.436875, 
                                days_to_last_followup/30.436875), 
           tissue = ifelse(sample_type == 'normal_tissue', 
                           'benign', 
                           'tumor') %>% 
             factor, 
           age = age_at_initial_pathologic_diagnosis, 
           gleason = gleason_score, 
           psa_at_diagnosis = psa, 
           pathology_stage_node = stri_extract(stage_event, 
                                               regex = 'N\\w{1}'), 
           pathology_stage_meta = stri_extract(stage_event, 
                                               regex = 'M\\w{1}'), 
           pathology_stage_tumor = stri_replace(stage_event, 
                                                regex = '^T', 
                                                replacement = ''), 
           pathology_stage_tumor = stri_extract(pathology_stage_tumor, 
                                                regex = 'T\\d{1}(a|b|c|d)'))
  
# renaming the expression dataset: the part after. is removed -----
  
  insert_msg('Renaming the expression dataset')
  
  names(tcga_data$expression) <- names(tcga_data$expression) %>% 
    stri_replace(regex = '\\..*', replacement = '')
  
  tcga_data$expression <- tcga_data$expression %>% 
    select( - file_id, 
            - case_id)
  
# integrating the expression data -----
  
  insert_msg('Integrating the TCGA data')
  
  tcga_cleared <- integrate_expression(expression_tbl = tcga_data$expression, 
                                          annotation_tbl = tcga_data$annotation %>% 
                                            set_names(c('probe_ID', 'symbol', 'protein')), 
                                          gene_identifier = 'symbol', 
                                          trans_fun = function(x) log2(x + 1), 
                                          .parallel = TRUE)
  
# removal of the duplicates: the most recent kept -----
  
  insert_msg('Duplicate removal')
  
  tcga_cleared$expression <- tcga_cleared$expression %>% 
    dlply(.(tissue), 
          rm_duplicates, 
          patient_identifier = 'patient_id', 
          reposit_identifier = 'entity_submitter_id') %>% 
    reduce(rbind)
  
# END ----
  
  insert_msg()