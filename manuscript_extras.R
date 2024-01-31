# Extra calculations and plots for the manuscript

# tools --------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)

  library(caret)
  library(caretExtra)
  library(bootStat)

  library(furrr)
  library(soucer)

  library(writexl)
  
  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  set_rownames <- trafo::set_rownames
  extract <- clustTools::extract
  
  insert_head()
  
  c('./tools/globals.R', 
    './tools/functions.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# Calculation of kappas for the proteomic models, heat maps ------
  
  insert_msg('Kappas for the proteome models, heat maps of confusion matrices')
  
  c('./extra scripts/kappas.R', 
    './extra scripts/pirads_hm.R') %>% 
    source_all(message = TRUE, crash = TRUE)

# END ------
  
  insert_tail()