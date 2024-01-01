# Source the complete analysis pipeline

  library(soucer)
  
  print(source_all(c('import.R', 
                     'exploration.R', 
                     'gleason.R', 
                     'survival.R', 
                     'clustering.R', 
                     'characteristic.R', 
                     'report.R'), 
                   message = TRUE, 
                   crash = TRUE))
  
  #save.image()