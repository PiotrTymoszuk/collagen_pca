# Project globals 

  insert_head()

# container --------

  globals <- list()
  
# genes of interest --------

  ## classification of the genes by their function
  
  globals$genes_interest <-
    c('ADAMTS2' = 'ECM processing', 
      'ALDH18A1' = 'Pro', 
      'BMP1' = 'ECM processing', 
      'CD151' = 'adhesion', 
      'COL11A1' = 'ECM component', 
      'COL11A2' = 'ECM component', 
      'COL14A1' = 'ECM component', 
      'COL15A1' = 'ECM component', 
      'COL16A1' = 'ECM component', 
      'COL17A1' = 'ECM component', 
      'COL18A1' = 'ECM component', 
      'COL19A1' = 'ECM component', 
      'COL1A1' = 'ECM component', 
      'COL1A2' = 'ECM component', 
      'COL2A1' = 'ECM component', 
      'COL3A1' = 'ECM component', 
      'COL4A1' = 'ECM component', 
      'COL4A2' = 'ECM component', 
      'COL4A3' = 'ECM component', 
      'COL4A5' = 'ECM component', 
      'COL4A6' = 'ECM component', 
      'COL5A1' = 'ECM component', 
      'COL5A2' = 'ECM component', 
      'COL6A1' = 'ECM component', 
      'COL6A2' = 'ECM component', 
      'COL6A3' = 'ECM component', 
      'COL7A1' = 'ECM component', 
      'COL9A1' = 'ECM component', 
      'COL9A2' = 'ECM component', 
      'COL9A3' = 'ECM component', 
      'CTSS' = 'ECM processing', 
      'DST' = 'adhesion', 
      'ITGA6' = 'adhesion', 
      'ITGB4' = 'adhesion', 
      'LAMA3' = 'ECM component', 
      'LAMB3' = 'ECM component', 
      'LAMC2' = 'ECM component', 
      'LOX' = 'collagen modification', 
      'LOXL1' = 'collagen modification',
      'LOXL2' = 'collagen modification', 
      'MMP13' = 'ECM processing', 
      'MMP7' = 'ECM processing', 
      'MMP9' = 'ECM processing', 
      'P4HA1' = 'collagen modification', 
      'P4HA2' = 'collagen modification', 
      'P4HB' = 'collagen modification', 
      'PCOLCE' = 'ECM processing', 
      'PCOLCE2' = 'ECM processing', 
      'PEPD' = 'Pro', 
      'PLOD1' = 'collagen modification', 
      'PLOD2' = 'collagen modification', 
      'PLOD3' = 'collagen modification', 
      'PPIB' = 'collagen modification', 
      'PYCR1' = 'Pro', 
      'SERPINH1' = 'ECM processing') %>% 
    compress(names_to = 'gene_symbol', 
             values_to = 'gene_group')
  
# clinical variables -------
  
  globals$clinical_lexicon <- 
    c('age' = 'Age at diagnosis, years', 
      'psa_diagnosis' = 'PSA at diagnosis', 
      'ct_stage' = 'Clinical tumor stage', 
      'pt_stage' = 'Pathological tumor stage', 
      'pn_stage' = 'Pathological node stage', 
      'pm_stage' = 'Pathological metastasis stage', 
      'gleason_sum' = 'Gleason score', 
      'gleason_simple' = 'Gleason score', 
      'surgical_margins' = 'Surgical margins', 
      'ece' = 'Extracapsular extension', 
      'death' = 'Death', 
      'os_months' = 'Overall survival, months', 
      'relapse' = 'Biochemical relapse', 
      'rfs_months' = 'Biochemical relapse-free survival, months') %>% 
    compress(names_to = 'variable', 
             values_to = 'label') %>% 
    mutate(format = ifelse(variable %in% c('age', 
                                           'psa_diagnosis', 
                                           'gleason', 
                                           'vitality_fup', 
                                           'relapse_fup'), 
                           'numeric', 'factor'))
  
# graphics --------
  
  ## ggplot theme
  
  globals$common_text <- element_text(size = 8, 
                                      face = 'plain', 
                                      color = 'black')
  
  globals$common_margin <- ggplot2::margin(t = 3, 
                                           l = 3, 
                                           r = 3, 
                                           unit = 'mm')
  
  globals$common_theme <- theme_classic() + 
    theme(axis.text = globals$common_text, 
          axis.title = globals$common_text, 
          plot.title = element_text(size = 8, 
                                    face = 'bold'), 
          plot.subtitle = globals$common_text, 
          plot.tag = element_text(size = 8, 
                                  face = 'plain', 
                                  color = 'black', 
                                  hjust = 0, 
                                  vjust = 1), 
          plot.tag.position = 'bottom', 
          legend.text = globals$common_text, 
          legend.title = globals$common_text, 
          strip.text = globals$common_text,
          strip.background = element_rect(fill = 'white'), 
          plot.margin = globals$common_margin, 
          panel.grid.major = element_line(color = 'gray90'))
  
# Study labels and colors --------
  
  ## GSE116918 is excluded fromm the analysis
  ## due to missing genes (LOXL1)
  
  globals$study_labels <- c('gse16560' = 'GSE16560', 
                            'gse54460' = 'GSE54460', 
                            'gse70768' = 'GSE70768', 
                            'gse70769' = 'GSE70769', 
                            #'gse116918' = 'GSE116918',
                             'gse220095' = 'GSE220095', 
                            'tcga' = 'TCGA', 
                            'dkfz' = 'DKFZ')
  
  globals$study_colors <- c('gse16560' = 'cornflowerblue', 
                            'gse54460' = 'firebrick', 
                            'gse70768' = 'darkseagreen4', 
                            'gse70769' = 'darkorange2', 
                            #'gse116918' = 'darkgreen', 
                            'gse220095' = 'cornsilk4', 
                            'tcga' = 'coral3', 
                            'dkfz' = 'plum4')
  
  globals$study_arrays <- c('gse16560' = TRUE, 
                            'gse54460' = FALSE, 
                            'gse70768' = TRUE, 
                            'gse70769' = TRUE, 
                            #'gse116918' = TRUE,
                            'gse220095' = FALSE, 
                            'tcga' = FALSE, 
                            'dkfz' = FALSE)
  
  ## a ready to use expression for calling a cohort list
  
  globals$study_exprs <- names(globals$study_labels) %>% 
    map_chr(~paste(.x, .x, sep = ' = ')) %>% 
    paste(collapse = ', ') %>% 
    paste0('list(', ., ')') %>% 
    parse_expr

# Cluster colors -------
  
  globals$cluster_colors <-
    c('Collagen hi' = 'coral3', 
      'Collagen int' = 'gray60', 
      'Collagen low' = 'steelblue3')
  
# END -----
  
  insert_tail()