# Renders the manuscript and the supplementary material

  insert_head()
  
# reading the bibliography -------
  
  insert_msg('Raeding the bibliography')
  
  coll_bib <- read_bib('./report/markdown/coll_biblio.bib')
  
# supplementary material, analysis report ------
  
  insert_msg('Rendering the analysis report')
  
  ## pending, I'll include the new figures and tables
  
 # render('./report/markdown/report_supplement.Rmd', 
  #       output_format = word_document2(number_sections = FALSE, 
   #                                     reference_docx = 'ms_template.docx'), 
    #     output_dir = './report')
  
#  render('./report/markdown/report.Rmd', 
 #        output_format = my_word(), 
  #       output_dir = './report')

# analysis report -----
  
  insert_msg('Rendering manuscript parts')
  
  render('./report/markdown/manuscript_figures_tables.Rmd', 
         output_format = word_document2(number_sections = FALSE, 
                                        reference_docx = 'ms_template.docx'), 
         output_dir = './report')
  

  
# END -----
  
  insert_tail()