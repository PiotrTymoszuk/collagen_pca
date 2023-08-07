# Project globals 

  insert_head()

# tools -------

  library(plyr)
  library(tidyverse)

# container --------

  globals <- list()
  
# genes of interest --------
  
  ## genes of interest: Kocher et al and Reactome R-HSA-1474290
  
  globals$genes_interest$kocher <- './data/genes_interest.xlsx' %>% 
    read_excel
  
  globals$genes_interest$reactome <- 
    './data/Participating Molecules [R-HSA-1474290].tsv' %>% 
    read_tsv %>% 
    mutate(gene_symbol = stri_split_fixed(MoleculeName, 
                                          pattern = ' ', 
                                          simplify = T)[, 2], 
           entrez_id = NA) %>% ## will be filled in if needed
    select(gene_symbol,
           entrez_id)
  
  globals$genes_interest <- globals$genes_interest %>% 
    reduce(rbind) %>% 
    filter(!duplicated(gene_symbol))
  
  ## filtering the gene list for the features present in all data sets
  
  globals$present_genes <- study_data %>% 
    map(~.x$expression) %>% 
    map(names) %>% 
    map(intersect, 
        globals$genes_interest$gene_symbol) %>% 
    reduce(intersect)
  
  globals$genes_interest <- globals$genes_interest %>% 
    filter(gene_symbol %in% globals$present_genes)
  
  ## classification of the genes by their function
  
  globals$gene_classes <-
    c('ALDH18A1' = 'Pro', 
      'PYCR1' = 'Pro', 
      'PLOD2' = 'collagen modification', 
      'P4HA1' = 'collagen modification', 
      'P4HA2' = 'collagen modification', 
      'PEPD' = 'Pro', 
      'COL3A1' = 'ECM component', 
      'COL2A1' = 'ECM component', 
      'COL11A1' = 'ECM component', 
      'COL11A2' = 'ECM component', 
      'COL5A2' = 'ECM component', 
      'COL5A1' = 'ECM component', 
      'COL1A2' = 'ECM component', 
      'COL7A1' = 'ECM component', 
      'BMP1' = 'ECM processing', 
      'LAMA3' = 'ECM component', 
      'LAMC2' = 'ECM component', 
      'LAMB3' = 'ECM component', 
      'COL4A2' = 'ECM component', 
      'COL4A1' = 'ECM component', 
      'COL4A5' = 'ECM component', 
      'COL4A3' = 'ECM component', 
      'COL4A6' = 'ECM component', 
      'MMP13' = 'ECM processing', 
      'MMP9' = 'ECM processing', 
      'CTSS' = 'ECM processing', 
      'MMP7' = 'ECM processing', 
      'COL18A1' = 'ECM component', 
      'COL15A1' = 'ECM component', 
      'COL14A1' = 'ECM component', 
      'COL6A3' = 'ECM component', 
      'COL6A1' = 'ECM component', 
      'COL6A2' = 'ECM component', 
      'DST' = 'adhesion', 
      'COL17A1' = 'ECM component', 
      'ITGA6' = 'adhesion', 
      'ITGB4' = 'adhesion', 
      'CD151' = 'adhesion', 
      'COL9A2' = 'ECM component', 
      'COL9A1' = 'ECM component', 
      'COL9A3' = 'ECM component', 
      'PCOLCE' = 'ECM processing', 
      'LOXL2' = 'collagen modification', 
      'LOXL1' = 'collagen modification', 
      'LOX' = 'collagen modification', 
      'COL16A1' = 'ECM component', 
      'COL19A1' = 'ECM component', 
      'PLOD3' = 'collagen modification', 
      'PLOD1' = 'collagen modification', 
      'P4HB' = 'collagen modification', 
      'ADAMTS2' = 'ECM processing', 
      'PPIB' = 'collagen modification', 
      'SERPINH1' = 'ECM processing', 
      'PCOLCE2' = 'ECM processing', 
      'COL1A1' = 'ECM component') %>% 
    compress(names_to = 'gene_symbol', 
             values_to = 'gene_group')
  
  globals$genes_interest <- globals$genes_interest %>% 
    left_join(globals$gene_classes, by = 'gene_symbol')
  
# clinical variables -------
  
  globals$clinical_lexicon <- 
    c('age' = 'Age at diagnosis, years', 
      'psa_at_diagnosis' = 'PSA at diagnosis', 
      'clinical_stage' = 'Clinical stage', 
      'pathology_stage_tumor' = 'Pathological tumor stage', 
      'pathology_stage_node' = 'Pathological node stage', 
      'pathology_stage_meta' = 'Pathological metastasis stage', 
      'gleason' = 'Gleason score', 
      'gleason_factor' = 'Gleason score', 
      'positive_surgical_margins' = 'Positive surgical margins', 
      'extra_capsular_extension' = 'Extracapsular extension', 
      'death' = 'Death', 
      'vitality_fup' = 'Overall survival, months', 
      'relapse' = 'Biochemical relapse', 
      'relapse_fup' = 'Biochemical relapse-free survival, months') %>% 
    compress(names_to = 'variable', 
             values_to = 'label') %>% 
    mutate(format = ifelse(variable %in% c('age', 
                                           'psa_at_diagnosis', 
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
  
  globals$study_labels <- c('GSE16560' = 'GSE16560', 
                            'GSE70768' = 'GSE70768', 
                            'GSE70769' = 'GSE70769', 
                            'GSE116918' = 'GSE116918', 
                            'tcga' = 'TCGA')
  
  globals$study_colors <- c('GSE16560' = 'cornflowerblue', 
                            'GSE70768' = 'darkseagreen4', 
                            'GSE70769' = 'darkorange2', 
                            'GSE116918' = 'cornsilk4', 
                            'tcga' = 'coral3')
  
# Cluster colors -------
  
  globals$cluster_colors <-
    c('Collagen hi' = 'coral3', 
      'Collagen int' = 'gray60', 
      'Collagen low' = 'steelblue3')
  
# END -----
  
  insert_tail()