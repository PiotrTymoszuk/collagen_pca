# Report tables 

  insert_head()
  
# containers ------
  
  report_tbl <- list()
  suppl_report_tbl <- list()
  
# Table 1: cohort characteristic --------
  
  insert_msg('Table 1: cohort characteristic')
  
  report_tbl$cohorts <- cohorts$stats %>% 
    set_names(c('Variable', globals$study_labels)) %>% 
    mdtable(label = 'table_1_study_cohorts', 
            ref_name = 'cohorts', 
            caption = paste('Characteristic of the analyzed cohorts.', 
                            'Numeric variables are presented as medians with', 
                            'interquartile rnges (IQR) and ranges.', 
                            'Qualitative variables are presented as', 
                            'percentages of categories within the', 
                            'complete observation set.'))
  
# Table 2: collagen genes of interest -----
  
  insert_msg('Table 2: collagen genes of interest')
  
  report_tbl$genes <- globals$genes_interest %>% 
    mutate(entrez_id = mapIds(org.Hs.eg.db, 
                              keys = gene_symbol, 
                              column = 'ENTREZID', 
                              keytype = 'SYMBOL')) %>% 
    set_names(c('Gene symbol', 'Entrez ID', 'Gene group')) %>% 
    mdtable(label = 'table_2_genes', 
            ref_name = 'genes', 
            caption = paste('Collagen genes of interest', 
                            'and their classification.'))
  
# Table S1 and S2: differentially regulated genes ------
  
  insert_msg('Table S1 and S2: DGE, collagen clusters')
  
  suppl_report_tbl[c("dge_int", "dge_hi")] <- 
    dge[c("dge_collagen_int", "dge_collagen_hi")] %>% 
    map(compress, names_to = 'cohort') %>% 
    map(mutate, cohort = globals$study_labels[cohort]) %>% 
    map(select, 
        cohort, 
        gene_symbol, 
        entrez_id, 
        regulation, 
        estimate, 
        lower_ci, 
        upper_ci, 
        p_adjusted) %>% 
    map(set_names, 
        c('Cohort', 
          'Gene symbol', 
          'Entrez ID', 
          'Regulation', 
          'log2 fold-regulation', 
          'lower 95% CI', 
          'upper 95% CI', 
          'p value')) %>% 
    map(map_dfc, function(x) if(is.numeric(x)) signif(x, 3) else x)
  
  suppl_report_tbl[c("dge_int", "dge_hi")] <- 
    suppl_report_tbl[c("dge_int", "dge_hi")] %>% 
    list(x = ., 
         label = c('table_s1_dge_collagen_int', 
                   'table_s2_dge_collagen_hi'), 
         ref_name = names(.), 
         caption = c(paste('Genes differentially regulated between', 
                           'Collagen intermediate and', 
                           'Collagen low cluster tumors.', 
                           'The table is available as a supplementary', 
                           'Excel sheet.'), 
                     paste('Genes differentially regulated between', 
                           'Collagen high and', 
                           'Collagen low cluster tumors.', 
                           'The table is available as a supplementary', 
                           'Excel sheet.'))) %>% 
    pmap(mdtable)

# Table S3 and S4: regulated biochemical reactions -------
  
  insert_msg('Table S3 and S4: reactions, collagen clusters')
  
  suppl_report_tbl[c("react_int", "react_hi")] <- 
    meta$models[c("int", "hi")] %>% 
    map(map, components, 'regulation') %>% 
    map(map, filter, p_adjusted < 0.05) %>% 
    map(map, function(x) if(nrow(x) == 0) NULL else x) %>% 
    map(compact) %>%  
    map(map, 
        mutate, 
        react_name = annotate_bigg(react_id, 
                                   value = 'name', 
                                   annotation_db = Recon2D), 
        react_name = unlist(react_name)) %>% 
    map(compress, names_to = 'cohort') %>% 
    map(map_dfc, function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    map(mutate, cohort = globals$study_labels[cohort]) %>% 
    map(select, 
        cohort, 
        subsystem, 
        react_id, 
        react_name, 
        fold_reg, 
        lower_ci, 
        upper_ci, 
        p_value) %>% 
    map(set_names, 
        c('Cohort', 
          'Reaction group', 
          'BiGG ID', 
          'Name', 
          'Fold regulation', 
          'lower CI', 
          'upper CI', 
          'p value'))

  suppl_report_tbl[c("react_int", "react_hi")] <-
    suppl_report_tbl[c("react_int", "react_hi")] %>% 
    list(x = ., 
         label = c('table_s3_differentially_regulated reactions', 
                   'table_s4_differentially_regulated_reactions'), 
         ref_name = names(.), 
         caption = c(paste('Biochemical reactions whose activity was predicted', 
                           'to be significantly regulated in', 
                           'Collagen intermediate versus Collagen low', 
                           'cluster cancers.', 
                           'The table is available as a supplementary', 
                           'Excel sheet.'), 
                     paste('Biochemical reactions whose activity was predicted', 
                           'to be significantly regulated in', 
                           'Collagen high versus Collagen low', 
                           'cluster cancers.', 
                           'The table is available as a supplementary', 
                           'Excel sheet.'))) %>% 
    pmap(mdtable)
  
# Saving the supplementary tables ------
  
  insert_msg('Saving the supplementary tables')
  
  suppl_report_tbl$cover <- 
    tibble(Table = paste0('Supplementary Table S', 1:length(suppl_report_tbl)), 
           Caption = map_chr(suppl_report_tbl, attr, 'caption'))
  
  suppl_report_tbl <- 
    c(suppl_report_tbl["cover"], 
      suppl_report_tbl[names(suppl_report_tbl) != "cover"])
  
  suppl_report_tbl %>% 
    set_names(c('Cover', paste0('Table S', 1:(length(suppl_report_tbl) - 1)))) %>% 
    write_xlsx(path = './report/report_supplementary_tables.xlsx')
  
# END ------
  
  insert_tail()