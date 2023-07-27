# Analysis report and paper scripts 

# tools ----

  library(plyr)
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

  select <- dplyr::select
  reduce <- purrr::reduce
  width <- flextable::width

  insert_head()
  
  source_all('./tools/project_tools.R', 
             message = TRUE, crash = TRUE)
  
# sourcing the scripts -----
  
  insert_msg('Paper scripts')
  
  ## files for the analysis report
  
  c('./report scripts/links.R', 
    './report scripts/report_tables.R', 
    './report scripts/report_figures.R', 
    './report scripts/report_supplement.R') %>% 
    source_all(message = TRUE, crash = TRUE) %>% 
    print
  
  ## files for the manuscript
  
  c('./report scripts/links.R', 
    './report scripts/paper_tables.R', 
    './report scripts/paper_figures.R', 
    './report scripts/paper_supplement.R') %>% 
    source_all(message = TRUE, crash = TRUE) %>% 
    print
  
  ## rendering the analysis report and parts of the manuscript
  
  c('./report scripts/render.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END -----
  
  insert_tail()