# Renders the manuscript and the supplementary material

  insert_head()
  
# reading the bibliography -------
  
  insert_msg('Raeding the bibliography')
  
  coll_bib <- read_bib('./report/markdown/coll_biblio.bib')
  
# supplementary material, analysis report ------
  
  insert_msg('Rendering the supplements for the analysis report')
  
  render('./report/markdown/report_supplement.Rmd', 
         output_format = word_document2(number_sections = FALSE, 
                                        reference_docx = 'ms_template.docx'), 
         output_dir = './report')
  
# analysis report -----
  
  insert_msg('Rendering the analysis report')
  
  my_word <- function(...) {
    
    form <- word_document2(number_sections = FALSE, 
                           reference_docx = 'ms_template.docx')
    
    form$pandoc$lua_filters <- c(form$pandoc$lua_filters, 
                                 'scholarly-metadata.lua', 
                                 'author-info-blocks.lua')
    
    form
    
  }
  
  render('./report/markdown/report.Rmd', 
         output_format = my_word(), 
         output_dir = './report')
  
# END -----
  
  insert_tail()