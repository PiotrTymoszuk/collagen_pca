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
  
  globals$genes_interest <- globals$genes_interest %>% 
    mutate(gene_group = ifelse(gene_symbol %in% c('ALDH18A1', 
                                                  'PYCR1', 
                                                  'PYCR2', 
                                                  'PEPD'), 
                               'proline pathway', 
                               ifelse(stri_detect(gene_symbol, regex = 'COL\\d+') | 
                                        gene_symbol %in% c('LAMB3', 'ITGA6', 'CD151'), 
                                      'ECM component', 
                                      'ECM processing')))
  
# clinical variables -------
  
  globals$clinical_lexicon <- 
    c('age' = 'Age at diagnosis, years', 
      'psa_at_diagnosis' = 'PSA at diagnosis', 
      'clinical_stage' = 'Clinical stage', 
      'pathology_stage_tumor' = 'Pathological tumor stage', 
      'pathology_stage_node' = 'Pathological node stage', 
      'pathology_stage_meta' = 'Pathological metastasis stage', 
      'gleason' = 'Gleason sum score', 
      'gleason_factor' = 'Gleason sum score', 
      'positive_surgical_margins' = 'Positive surgical margins', 
      'extra_capsular_extension' = 'Extracapsular extension', 
      'death' = 'Death', 
      'vitality_fup' = 'Overall survival, months', 
      'relapse' = 'Relapse', 
      'relapse_fup' = 'Relapse-free survival, months') %>% 
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
                            'GSE40272' = 'GSE40272', 
                            'GSE70768' = 'GSE70768', 
                            'GSE70769' = 'GSE70769', 
                            'tcga' = 'TCGA')
  
  globals$study_colors <- c('GSE16560' = 'cornflowerblue', 
                            'GSE40272' = 'cornsilk4', 
                            'GSE70768' = 'darkseagreen4', 
                            'GSE70769' = 'darkorange2', 
                            'tcga' = 'coral3')
  
# Cluster colors -------
  
  globals$cluster_colors <-
    c('Collagen hi' = 'coral3', 
      'Collagen int' = 'gray60', 
      'Collagen low' = 'steelblue3')
  
# END -----
  
  insert_tail()