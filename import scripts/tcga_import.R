# This script imports RNA and clinical data for the TCGA PRAD (prostate adenocarcinoma) project:
# (1) An experiment list is generated containing the clinical and expression data
# (2) a table with gene and transcript identifiers for the genes of interest: CAND1, SKP1 and SKP2
# (3) a working table containing the entire clinical information, sample type data (see: $assignment)
# and expression data for the transcripts of interest

  insert_head()

# reading the experiment -----

  enter_directory('./input data/tcga')

  tcga_data <- read_experiment(assignment_data_folder = 'sample_assignment', 
                               expression_data_folder = 'expression', 
                               clinical_data_folder = 'clinical')
  
# cleaning clinical data with parse guess -----
  
  tcga_data$clinical <- tcga_data$clinical %>% 
    map_dfc(as.character) %>% 
    map_dfc(parse_guess)
  
# reading gene symbol annotation ----
  
  tcga_data$annotation <- read_tsv('annotation.txt') %>% 
    set_names(c('gene_id', 
                'gene_symbol', 
                'protein_id')) %>% 
    filter(!duplicated(gene_id))
  
# END ----
  
  go_proj_directory()
  
  insert_tail()