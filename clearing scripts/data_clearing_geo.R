# This script clears the GEO data -----

  insert_head()

# GSE40272: merging the sub-experiments ----

  insert_msg('GSE40272: merging the sub-experiments')

  geo_data$GSE40272$expression <- geo_data[c('GSE40272-GPL15971', 
                                             'GSE40272-GPL15972', 
                                             'GSE40272-GPL15973', 
                                             'GSE40272-GPL9497')] %>% 
    map(function(x) x$expression) %>% 
    map_dfr(mutate, 
            dfs_months = as.numeric(dfs_months))
  
  geo_data$GSE40272$annotation <- geo_data$`GSE40272-GPL15971`$annotation

# integrating the expression as geometric mean of multiple probes -----

  insert_msg('Integrating the gene expression')
  
  study_data <- geo_data[c('GSE16560', 
                           'GSE40272', 
                           'GSE70768', 
                           'GSE70769')] %>% 
    map(function(x) integrate_expression(expression_tbl = x$expression, 
                                         annotation_tbl = x$annotation, 
                                         gene_identifier = 'symbol', 
                                         .parallel = T))
  
# clearing the GSE16560 study metadata ----
  
  insert_msg('Clearing the GSE16560 metadata')
  
  study_data$GSE16560$expression <- study_data$GSE16560$expression %>% 
    mutate(death = ifelse(vitality == 'Alive', 0, 1), 
           vitality_fup = followup_months, 
           patient_id = provided_ID)
  
# clearing the GSE40272 study metadata ----
  
  insert_msg('Clearing the GSE40272 metadata')
  
  study_data$GSE40272$expression <- study_data$GSE40272$expression %>% 
    mutate(relapse = 1, 
           relapse_fup = dfs_months, 
           tissue = ifelse(stri_detect(individual, regex = 'T$'), 
                           'tumor', 
                           'benign') %>% 
             factor, 
           patient_id = stri_extract(individual, regex = '\\d+'))
  
  ## appending the expression data set with clinical data found in the paper supplement
  
  study_data$GSE40272$supplement <- read_excel('./input data/GSE40272-GPL9497/NIHMS349620-supplement-4.xls') %>% 
    mutate(individual = Sample, 
           age = Age, 
           ethnicity = Ethnicity, 
           psa_at_diagnosis = `Pre-PSA`, 
           pathology_stage_tumor = `T`,
           pathology_stage_node = N, 
           pathology_stage_meta = M0, 
           gleason_major = stri_split_fixed(`Path Gr`, pattern = '+', simplify = T)[, 1] %>% 
             as.numeric, 
           gleason_minor = stri_split_fixed(`Path Gr`, pattern = '+', simplify = T)[, 2] %>% 
             as.numeric, 
           gleason = gleason_major + gleason_minor, 
           positive_surgical_margins = ifelse(Margins == 'Positive', 
                                              'yes', 'no'), 
           relapse = ifelse(Recurr == 'Biochemical', 
                            1, 0), 
           relapse_fup = `DFS Months`) %>% 
    select(individual, 
           age, 
           ethnicity, 
           psa_at_diagnosis, 
           pathology_stage_tumor, 
           pathology_stage_meta, 
           pathology_stage_node, 
           gleason, 
           positive_surgical_margins, 
           relapse, 
           relapse_fup)
  
  ## merging with the expression data set
  
  study_data$GSE40272$expression <- left_join(study_data$GSE40272$expression %>% 
                                                select(- age, 
                                                       - ethnicity, 
                                                       - relapse, 
                                                       - relapse_fup), 
                                              study_data$GSE40272$supplement, 
                                              by = 'individual')
  ## removal of the duplicates
  
  study_data$GSE40272$expression  <- study_data$GSE40272$expression %>% 
    dlply(.(tissue), rm_duplicates) %>% 
    reduce(rbind)
  
# clearing the GSE70768 study metadata ----
  
  insert_msg('Clearing the GSE70768 metadata')
  
  study_data$GSE70768$expression <- study_data$GSE70768$expression %>% 
    mutate(age = age_at_diagnosis, 
           tissue = ifelse(sample_type == 'Tumour', 
                           'tumor', 
                           'benign') %>% 
             factor, 
           extra_capsular_extension = car::recode(extra_capsular_extension, 
                                                  "'unknown' = NA;
                                                  'N' = 'no'; 
                                                  'Y' = 'yes'") %>% 
             factor(c('no', 'yes')), 
           positive_surgical_margins = car::recode(positive_surgical_margins, 
                                                   "'unknown' = NA;
                                                  'N' = 'no'; 
                                                  'Y' = 'yes'") %>% 
             factor(c('no', 'yes')), 
           relapse = car::recode(biochemical_relapse, 
                                 "'N/A' = NA; 
                                 'N' = 0; 
                                 'Y' = 1") %>% 
             as.numeric, 
           relapse_fup = ifelse(relapse == 1, 
                                as.numeric(time_to_bcr_months), 
                                as.numeric(total_follow_up_months)), 
           clinical_stage = stri_split_fixed(clinical_stage, pattern = ' ', simplify = T)[, 1], 
           pathology_stage_tumor = ifelse(pathology_stage != 'UNKNOWN', 
                                          stri_split_fixed(pathology_stage, pattern = ' ', simplify = T)[, 1],
                                          NA), 
           pathology_stage_node = ifelse(pathology_stage != 'UNKNOWN', 
                                         stri_extract(pathology_stage, regex = 'N\\w{1}'), 
                                         NA), 
           pathology_stage_meta = ifelse(pathology_stage != 'UNKNOWN', 
                                         stri_extract(pathology_stage, regex = 'M\\w{1}'), 
                                         NA), 
           patient_id = stri_extract(provided_ID , regex = 'TB\\w+\\.\\w+'))
  
  ## removal of the duplicates
  
  study_data$GSE70768$expression  <- study_data$GSE70768$expression %>% 
    dlply(.(tissue), rm_duplicates) %>% 
    reduce(rbind)
  
# clearing the GSE70769 study metadata ----
  
  insert_msg('Clearing the GSE70769 metadata')
  
  study_data$GSE70769$expression <- study_data$GSE70769$expression %>% 
    mutate(extra_capsular_extension = car::recode(extra_capsular_extension, 
                                                  "'unknown' = NA;
                                                  'N' = 'no'; 
                                                  'Y' = 'yes'") %>% 
             factor(c('no', 'yes')), 
           positive_surgical_margins = car::recode(positive_surgical_margins, 
                                                   "'unknown' = NA;
                                                  'N' = 'no'; 
                                                  'Y' = 'yes'") %>% 
             factor(c('no', 'yes')), 
           relapse = car::recode(biochemical_relapse, 
                                 "'N/A' = NA; 
                                 'N' = 0; 
                                 'Y' = 1") %>% 
             as.numeric, 
           relapse_fup = ifelse(relapse == 1, 
                                as.numeric(time_to_bcr_months), 
                                as.numeric(total_follow_up)), 
           clinical_stage = ifelse(clinical_stage %in% c('UNKNOWN', ':'), 
                                   NA, 
                                   clinical_stage), 
           pathology_stage_tumor = stri_extract(pathology_stage, 
                                                regex = 'T\\d{1}(a|b|c|d)'), 
           pathology_stage_tumor = paste0('p', pathology_stage_tumor), 
           pathology_stage_node = stri_extract(pathology_stage, 
                                                regex = 'N\\w{1}'), 
           pathology_stage_meta = stri_extract(pathology_stage, 
                                               regex = 'M\\w{1}'), 
           patient_id = stri_extract(provided_ID, 
                                     regex = 'STKHLM.*$'))
  
  study_data$GSE70769$expression <- study_data$GSE70769$expression %>% 
    rm_duplicates

# END ----
  
  insert_msg()