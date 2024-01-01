# Analysis report and paper scripts 

# tools ----

  library(tidyverse)
  library(exda)
  library(figur)
  library(rmarkdown)
  library(knitr)
  library(bookdown)
  library(flextable)
  library(writexl)
  library(soucer)
  library(stringi)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(pathview)
  library(coxExtensions)
  library(rlang)
  library(survival)

  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  set_rownames <- trafo::set_rownames
  extract <- clustTools::extract
  
  insert_head()
  
  c('./tools/globals.R', 
    './tools/functions.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# sourcing the scripts -----
  
  insert_msg('Paper scripts')

  c('./report scripts/links.R', 
    './report scripts/paper_tables.R', 
    './report scripts/paper_figures.R', 
    './report scripts/paper_supplement.R') %>% 
    source_all(message = TRUE, crash = TRUE) %>% 
    print
  
  ## rendering parts of the manuscript
  
  c('./report scripts/render.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END -----
  
  insert_tail()