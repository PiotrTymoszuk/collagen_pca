# This script clears the GEO data -----

  insert_head()

# container ------

  geo_cleared <- list()

# integrating the expression as geometric mean of multiple probes -----

  insert_msg('Integrating the gene expression')
  
  geo_cleared <- geo_data %>% 
    map(function(x) integrate_expression(expression_tbl = x$expression, 
                                         annotation_tbl = x$annotation, 
                                         gene_identifier = 'symbol', 
                                         .parallel = TRUE))
  
# clearing the GSE16560 study metadata ----
  
  insert_msg('Clearing the GSE16560 metadata')
  
  geo_cleared$GSE16560$clinic <- geo_data$GSE16560$clinic %>% 
    transmute(sample_ID = geo_accession, 
              patient_id = title, 
              gleason = extract_pfield(characteristics_ch1, 'numeric'), 
              major_gleason = extract_pfield(characteristics_ch1.1, 'numeric'), 
              minor_gleason = extract_pfield(characteristics_ch1.2, 'numeric'), 
              cancer_percent = extract_pfield(characteristics_ch1.5, 'numeric'), 
              age = extract_pfield(characteristics_ch1.6, 'numeric'), 
              vitality = extract_pfield(characteristics_ch1.8, 'factor', 
                                        c('Alive', 'Dead')), 
              death = ifelse(vitality == 'Alive', 0, 1), 
              followup_months = extract_pfield(characteristics_ch1.9, 'numeric'), 
              vitality_fup = followup_months, 
              tissue = factor('tumor', c('benign', 'tumor')))
  
  geo_cleared$GSE16560$expression <- 
    left_join(geo_cleared$GSE16560$clinic, 
              geo_cleared$GSE16560$expression, 
              by = 'sample_ID')
  
# clearing the GSE70768 study metadata ----
  
  insert_msg('Clearing the GSE70768 metadata')
  
  geo_cleared$GSE70768$clinic <- geo_data$GSE70768$clinic %>% 
    transmute(sample_ID = geo_accession, 
              patient_id = title, 
              patient_id = stri_split_fixed(patient_id, 
                                            pattern = '_', 
                                            simplify = TRUE)[, 3], 
              tissue = extract_pfield(characteristics_ch1, 'character'), 
              tissue = car::recode(tolower(tissue), "'tumour' = 'tumor'"), 
              tissue = factor(tolower(tissue), c('benign', 'tumor')), 
              gleason = extract_pfield(characteristics_ch1.1, 'character'), 
              gleason = ifelse(tissue == 'benign', NA, gleason), 
              gleason_minor = stri_extract(gleason, regex = '\\d{1}$'), 
              gleason_minor = as.numeric(gleason_minor), 
              gleason = stri_extract(gleason, regex = '^\\d{1,2}'), 
              gleason = as.numeric(gleason), 
              gleason_major = gleason - gleason_minor, 
              cancer_percent = extract_pfield(characteristics_ch1.2, 'character'), 
              cancer_percent = stri_extract(cancer_percent, regex = '\\d{1,3}'), 
              cancer_percent = as.numeric(cancer_percent), 
              extra_capsular_extension = extract_pfield(characteristics_ch1.4, 
                                                        'character'), 
              extra_capsular_extension = ifelse(tissue == 'benign', 
                                                NA, 
                                                extra_capsular_extension), 
              extra_capsular_extension = car::recode(extra_capsular_extension, 
                                                     "'N' = 'no'; 'Y' = 'yes'"), 
              extra_capsular_extension = factor(extra_capsular_extension, 
                                                c('no', 'yes')), 
              positive_surgical_margins = extract_pfield(characteristics_ch1.5,
                                                         'character'), 
              positive_surgical_margins = car::recode(positive_surgical_margins, 
                                                      "'N' = 'no'; 'Y' = 'yes'"), 
              positive_surgical_margins = factor(positive_surgical_margins, 
                                                 c('no', 'yes')), 
              biochemical_relapse = extract_pfield(characteristics_ch1.6, 
                                                   'character'), 
              biochemical_relapse = car::recode(biochemical_relapse, 
                                                "'N'= 'no'; 'Y' = 'yes'"), 
              biochemical_relapse = factor(biochemical_relapse, 
                                           c('no', 'yes')), 
              relapse = ifelse(biochemical_relapse == 'yes', 1, 0), 
              time_to_bcr_months = extract_pfield(characteristics_ch1.7, 'numeric'), 
              total_follow_up_months = extract_pfield(characteristics_ch1.13, 'numeric'), 
              relapse_fup = ifelse(relapse == 1, 
                                   as.numeric(time_to_bcr_months), 
                                   as.numeric(total_follow_up_months)), 
              age = extract_pfield(characteristics_ch1.9, 'numeric'), 
              psa_at_diagnosis = extract_pfield(characteristics_ch1.10, 'numeric'), 
              clinical_stage = extract_pfield(characteristics_ch1.11, 'character'), 
              clinical_stage = stri_extract(clinical_stage, regex = '^T\\d{1}'), 
              clinical_stage = factor(clinical_stage, c('T1', 'T2', 'T3', 'T4')), 
              pathology_stage_tumor = extract_pfield(characteristics_ch1.12, 'character'), 
              pathology_stage_node = stri_extract(pathology_stage_tumor, regex = 'N\\w{1}'), 
              pathology_stage_meta = stri_extract(pathology_stage_tumor, regex = 'M\\w{1}'), 
              pathology_stage_tumor = stri_extract(pathology_stage_tumor, regex = 'T\\d{1}'), 
              pathology_stage_tumor = factor(pathology_stage_tumor, c('T1', 'T2', 'T3', 'T4')), 
              pathology_stage_node = factor(pathology_stage_node, c('N0', 'N1')), 
              pathology_stage_meta = factor(pathology_stage_meta, c('M0', 'M1')))

  
  ## removal of the duplicates
  
  geo_cleared$GSE70768$expression  <- 
    left_join(geo_cleared$GSE70768$clinic, 
              geo_cleared$GSE70768$expression,
              by = 'sample_ID') %>% 
    blast(tissue) %>% 
    map_dfr(rm_duplicates)
  
# clearing the GSE70769 study metadata ----
  
  insert_msg('Clearing the GSE70769 metadata')
  
  geo_cleared$GSE70769$clinic <- geo_data$GSE70769$clinic %>% 
    transmute(sample_ID = geo_accession,
              patient_id = stri_split_fixed(title, pattern = '_', simplify = TRUE)[, 3], 
              tissue = factor('tumor', c('benign', 'tumor')), 
              gleason = extract_pfield(characteristics_ch1, 'character'), 
              gleason_minor = stri_extract(gleason, regex = '\\d{1}$'), 
              gleason_minor = as.numeric(gleason_minor), 
              gleason = stri_extract(gleason, regex = '^\\d{1,2}'), 
              gleason = as.numeric(gleason), 
              gleason_major = gleason - gleason_minor, 
              cancer_percent = extract_pfield(characteristics_ch1.1, 'character'), 
              cancer_percent = stri_extract(cancer_percent, regex = '\\d{1,3}'), 
              cancer_percent = as.numeric(cancer_percent), 
              extra_capsular_extension = extract_pfield(characteristics_ch1.2, 
                                                        'character'), 
              extra_capsular_extension = car::recode(extra_capsular_extension, 
                                                     "'N' = 'no'; 'Y' = 'yes'"), 
              extra_capsular_extension = factor(extra_capsular_extension, 
                                                c('no', 'yes')), 
              positive_surgical_margins = extract_pfield(characteristics_ch1.3,
                                                         'character'), 
              positive_surgical_margins = car::recode(positive_surgical_margins, 
                                                      "'N' = 'no'; 'Y' = 'yes'"), 
              positive_surgical_margins = factor(positive_surgical_margins, 
                                                 c('no', 'yes')), 
              biochemical_relapse = extract_pfield(characteristics_ch1.4, 
                                                   'character'), 
              biochemical_relapse = car::recode(biochemical_relapse, 
                                                "'N'= 'no'; 'Y' = 'yes'"), 
              biochemical_relapse = factor(biochemical_relapse, 
                                           c('no', 'yes')), 
              relapse = ifelse(biochemical_relapse == 'yes', 1, 0), 
              time_to_bcr_months = extract_pfield(characteristics_ch1.5, 'numeric'), 
              total_follow_up_months = extract_pfield(characteristics_ch1.10, 'numeric'), 
              relapse_fup = ifelse(relapse == 1, 
                                   as.numeric(time_to_bcr_months), 
                                   as.numeric(total_follow_up_months)), 
              psa_at_diagnosis = extract_pfield(characteristics_ch1.7, 'numeric'), 
              clinical_stage = extract_pfield(characteristics_ch1.8, 'character'), 
              clinical_stage = stri_extract(clinical_stage, regex = '^T\\d{1}'), 
              clinical_stage = factor(clinical_stage, c('T1', 'T2', 'T3', 'T4')), 
              pathology_stage_tumor = extract_pfield(characteristics_ch1.9, 'character'), 
              pathology_stage_node = stri_extract(pathology_stage_tumor, regex = 'N\\w{1}'), 
              pathology_stage_meta = stri_extract(pathology_stage_tumor, regex = 'M\\w{1}'), 
              pathology_stage_tumor = stri_extract(pathology_stage_tumor, regex = 'T\\d{1}'), 
              pathology_stage_tumor = factor(pathology_stage_tumor, c('T1', 'T2', 'T3', 'T4')), 
              pathology_stage_node = factor(pathology_stage_node, c('N0', 'N1')), 
              pathology_stage_meta = factor(pathology_stage_meta, c('M0', 'M1')))

  
  geo_cleared$GSE70769$expression <- 
    left_join(geo_cleared$GSE70769$clinic, 
              geo_cleared$GSE70769$expression, 
              by = 'sample_ID') %>% 
    rm_duplicates
  
# clearing the GSE116918 study metadata ----
  
  insert_msg('Clearing the GSE116918 metadata')
  
  geo_cleared$GSE116918$clinic <- geo_data$GSE116918$clinic %>% 
    transmute(sample_ID = geo_accession, 
              patient_id = title, 
              tissue = factor('tumor', c('benign', 'tumor')), 
              age = extract_pfield(characteristics_ch1.1, 'numeric'), 
              psa_at_diagnosis = extract_pfield(characteristics_ch1.2, 
                                                'numeric'), 
              pathology_stage_tumor = extract_pfield(characteristics_ch1.3, 
                                                     'character'), 
              pathology_stage_tumor = stri_extract(pathology_stage_tumor, 
                                                   regex = 'T\\d{1}'), 
              pathology_stage_tumor = factor(pathology_stage_tumor, 
                                             c('T1', 'T2', 'T3', 'T4')), 
              gleason = extract_pfield(characteristics_ch1.4, 'numeric'), 
              relapse = extract_pfield(characteristics_ch1.5, 'numeric'), 
              biochemical_relapse = ifelse(relapse == 1, 'yes', 'no'), 
              biochemical_relapse = factor(biochemical_relapse, c('no', 'yes')), 
              relapse_fup = extract_pfield(characteristics_ch1.6, 'numeric'))
  
  geo_cleared$GSE116918$expression <- 
    left_join(geo_cleared$GSE116918$clinic,
              geo_cleared$GSE116918$expression, 
              by = 'sample_ID') %>% 
    rm_duplicates

# END ----
  
  insert_msg()