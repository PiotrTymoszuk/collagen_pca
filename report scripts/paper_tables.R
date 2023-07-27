# Main and supplementary tables for the manuscript 

  insert_head()
  
# containers ------
  
  paper_tbl <- list()
  suppl_paper_tbl <- list()
  
# Table 1: cohort characteristic --------
  
  insert_msg('Table 1: cohort characteristic')
  
  paper_tbl$cohorts <- cohorts$stats %>% 
    set_names(c('Variable', globals$study_labels)) %>% 
    mdtable(label = 'table_1_study_cohorts', 
            ref_name = 'cohorts', 
            caption = paste('Characteristic of the analyzed cohorts.', 
                            'Numeric variables are presented as medians with', 
                            'interquartile ranges (IQR) and ranges.', 
                            'Qualitative variables are presented as', 
                            'percentages of categories within the', 
                            'complete observation set.'))
  
# Table 2: collagen genes of interest -----
  
  insert_msg('Table 2: collagen genes of interest')
  
  paper_tbl$genes <- globals$genes_interest %>% 
    mutate(entrez_id = mapIds(org.Hs.eg.db, 
                              keys = gene_symbol, 
                              column = 'ENTREZID', 
                              keytype = 'SYMBOL')) %>% 
    set_names(c('Gene symbol', 'Entrez ID', 'Gene group')) %>% 
    mdtable(label = 'table_2_genes', 
            ref_name = 'genes', 
            caption = paste('Collagen genes of interest', 
                            'and their classification.'))
  
# Table S1: normal-tumor ------
  
  insert_msg('Table S1: normal - tumor')
  
  suppl_paper_tbl$norm_tumor <- 
    map2(norm_tumor$stats, 
         map(norm_tumor$test_results, 
             ~.x[c('variable', 'significance', 'est_lab')]), 
         left_join, by = 'variable') %>% 
    map(format_summ_tbl) %>% 
    map2(., norm_tumor$n_numbers, 
         ~full_rbind(tibble(variable = 'Samples, n', 
                            tumor = .y, 
                            benign = .y), 
                     .x)) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    select(cohort, variable, tumor, benign, significance, est_lab) %>% 
    set_names(c('Cohort', 'Variable', 'Tumor', 'Benign', 
                'Significance', 'Effect size')) %>% 
    mdtable(label = 'table_s1_normal_tumor', 
            ref_name = 'norm_tumor', 
            caption = paste("Expression of the collagen pathway genes", 
                            "in the malignant and benign tissue compared", 
                            "by paired T test with Cohen's d", 
                            "effect size statistic.", 
                            "P values were corrected for multiple testing", 
                            "with the false discovery rate method.", 
                            "log2-transformed expression values are presented", 
                            "as medians with interquartile ranges (IQR)", 
                            "and ranges.", 
                            "The table is available as a supplementary", 
                            "Excel file."))
  
# Table S2: expression of the collagen pathway genes in the clusters -------
  
  insert_msg('Table S2: collagen gene expression in the clusters')
  
  suppl_paper_tbl$cluster <- coll_clust$result_tbl %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    select(cohort, variable, 
           `Collagen low`, `Collagen int`, `Collagen hi`, 
           significance, eff_size) %>% 
    set_names(c('Cohort', 'Variable', 
                'Collagen low', 'collagen intermediate', 'Collagen high', 
                'Significance', 'Effect size')) %>% 
    mdtable(label = 'table_s2_cluster_collagen_genes', 
            ref_name = 'cluster', 
            caption = paste("Expression of the cluster-defining", 
                            "collagen pathway genes in the collagen clusters", 
                            "of prostate cancer.", 
                            "Statistical significance was assessed by one-way", 
                            "ANOVA with eta-squared effect size statistic.", 
                            "P values were corrected for multiple testing", 
                            "with the false discovery rate method.", 
                            "log2-transformed expression values are presented", 
                            "as medians with interquartile ranges (IQR)", 
                            "and ranges.", 
                            "The table is available as a supplementary", 
                            "Excel file."))
  
# Table S3: clinical characteristic of the collagen clusters -----
  
  insert_msg('Table S3: Clinical characteristic of the clusters')
  
  suppl_paper_tbl$clinic <- cs_cluster$result_tbl %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    filter(variable != 'Collagen Score') %>% 
    select(cohort, variable, low, int, hi, significance, eff_size) %>% 
    set_names(c('Cohort', 'Variable', 
                'Collagen low', 'Collagen intermediate', 'Collagen high', 
                'Significance', 'Effect size')) %>% 
    mdtable(label = 'table_s3_cluster_clinical_characteristic', 
            ref_name = 'clinic',
            caption = paste('Clinical characteristic of the collagen clusters.', 
                            'Numeric variables are presented as medians', 
                            'with interquartile ranges (IQR) and ranges.', 
                            'Nominal variables are presented as percentages', 
                            'and counts of categories within the cluster.'))
  
# Table S4 - S5: infiltration -------
  
  insert_msg('Table S4 - S5: MCP and xcell counter infiltration')
  
  for(i in names(clust_infil$stats)) {
    
    suppl_paper_tbl[[i]] <- 
      map2(clust_infil$stats[[i]], 
           map(clust_infil$test[[i]], 
               ~.x[c('variable', 'significance', 'eff_size')]), 
           left_join, by = 'variable') %>% 
      map2(., clust_infil$n_numbers[[i]], 
           ~full_rbind(tibble(variable = 'Samples, n', 
                              low = .y[['n']][1], 
                              int = .y[['n']][2], 
                              hi = .y[['n']][3]), 
                       .x)) %>% 
      compress(names_to = 'cohort') %>% 
      mutate(cohort = globals$study_labels[cohort]) %>% 
      format_summ_tbl %>% 
      select(cohort, variable, 
             low, int, hi, significance, eff_size) %>% 
      set_names(c('Cohort', 'Variable', 
                  'Collagen low', 'Collagen int', 'Collagen high', 
                  'Significance', 'Effect size'))
    
  }

  suppl_paper_tbl[c("mcp", "xcell")] <- 
    suppl_paper_tbl[c("mcp", "xcell")] %>% 
    list(x = ., 
         label = c('table_s4_mcp_counter_infiltration', 
                   'table_s5_xcell_infiltration'), 
         ref_name = names(.),
         caption = paste(c('Non-malignant cell numbers predicted for the collagen', 
                           'clusters by the MCP counter algorithm.'), 
                         c('Non-maignant cell fractions predicted for the collagen', 
                           'clusters by the xCell algorithm.'), 
                         'Statistical significance was assessed by Kruskal-Wallis',
                         'test with eta-squared effect size statistic.', 
                         'P values were corrected for multiple testing with the', 
                         'false discovery method.', 
                         'The table is available as a supplementary Excel file.')) %>% 
    pmap(mdtable)
  
# Table S6: GSVA ------
  
  insert_msg('Table S6: GSVA')
  
  suppl_paper_tbl$gsva <- 
    dge_gsva$lm_signif %>% 
    map_dfr(map_dfr, 
            compress, 
            names_to = 'cohort') %>% 
    filter(response %in% unique(unlist(dge_gsva$cmm_signatures))) %>% 
    left_join(compress(dge_gsva$eff_size, 
                       names_to = 'cohort'), 
              by = c('response', 'cohort')) %>% 
    left_join(set_names(dge_gsva$hm_plot$sign_classes, 
                        c('response', 'class')), 
              by = 'response') %>% 
    mutate(cohort = globals$study_labels[cohort], 
           response = exchange(response, 
                               reactome$lexicon), 
           level = as.character(level),
           level = paste('Collagen', level, 'vs low'),
           class = stri_replace_all(class, 
                                    fixed = '\n', 
                                    replacement = ', ')) %>% 
    select(cohort, class, response, level, regulation, 
           estimate, lower_ci, upper_ci, 
           p_adjusted, eta_sq) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    arrange(cohort, class, response) %>% 
    set_names(c('Cohort', 'Signature class', 
                'Signature', 'Comparison', 
                "Regulation", 
                'Fold-regulation', 
                'lower 95% CI', 
                'upper 95% CI', 
                'pFDR for fold-regulation', 
                'Effect size')) %>% 
    mdtable(label = 'table_s6_cluster_gsva_reactome', 
            ref_name = 'gsva', 
            caption = paste("Gene set variation", 
                            "analysis with the Reactome pathway gene signatures.", 
                            "Differences between the collagen intermediate", 
                            "or high clusters versus collagen low cancers", 
                            "were investigated by one-way ANOVA with", 
                            "eta-squared effect size statistic and linear", 
                            "modeling.", 
                            "Results for signatures significantly regulated",
                            "with moderate-to-large effect size", 
                            "(eta-squared at leat 0.06) in at least four", 
                            "cohorts are presented.", 
                            "P values were corrected for multiple testing with",
                            "the false discovery rate method (FDR).", 
                            "The table is available as a supplementary", 
                            "Excel file."))
  
# Table S7: differential gene expression -------
  
  insert_msg('Table S7: differential gene expression')
  
  suppl_paper_tbl$dge <- dge[c("dge_collagen_int", "dge_collagen_hi")] %>% 
    map_dfr(compress, names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort], 
           level = paste(level, 'vs low')) %>% 
    select(cohort, 
           gene_symbol, 
           entrez_id, 
           regulation, 
           estimate, 
           lower_ci, 
           upper_ci, 
           p_adjusted) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    set_names(c('Cohort', 
                'Gene symbol', 
                'Entrez ID', 
                'Regulation', 
                'log2 fold-regulation', 
                'lower 95% CI', 
                'upper 95% CI', 
                'pFDR')) %>% 
    mdtable(label = 'table_s4_clusters_differential_gene_expression', 
            ref_name = 'dge', 
            caption = paste('Genes differentially expressed in the collagen', 
                            'intermediate or high cluster as compared', 
                            'with collagen low cancers were identified', 
                            'by one-way ANOVA and linear modeling', 
                            'with the 1.25-fold regulation cutoff', 
                            'P values were corrected for multiple testing with', 
                            'the false discovery rate method (FDR).', 
                            'The table is available as a supplementary', 
                            'Excel file.'))
  
# Table S8: signaling -------
  
  insert_msg('Table S8: signaling')
  
  suppl_paper_tbl$signaling <- dge_spia$test %>% 
    map(compress, names_to = 'cohort') %>% 
    map(filter, Name %in% unique(unlist(dge_spia$common))) %>% 
    compress(names_to = 'level') %>% 
    mutate(cohort = globals$study_labels[cohort], 
           level = paste('Collagen', level, 'vs low'), 
           pGFDr = ifelse(pGFdr < 0.05, 
                          paste('p =', signif(pGFdr, 2)), 
                          paste0('ns (p = ', signif(pGFdr, 2), ')'))) %>% 
    select(cohort, Name, ID, level, Status, tA, pGFdr) %>% 
    set_names(c('Cohort', 'Pathway name', 'KEGG ID', 
                'Comparison', 'Activation status', 
                'Regulation, tA', 'Significance, pGFDR')) %>% 
    mdtable(label = 'table_s8_clusters_signaling_pathways', 
            ref_name = 'signaling',
            caption = paste('Signaling pathway activity in the collagen', 
                            'clusters investigated by the SPIA algorithm.', 
                            'Resulat for signaling pathways significantly', 
                            'activated or inhibited in at least four cohorts', 
                            'are shown.', 
                            'The table is available as a supplementary', 
                            'Excel file.'))
  
# Table S9: biochemical reactions -------
  
  insert_msg('Table S9: Biochemical reactions')
  
  suppl_paper_tbl$reactions <- 
    meta$models %>% 
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
    compress(names_to = 'level') %>% 
    filter(p_adjusted < 0.05, 
           fold_reg != 1) %>% 
    mutate(cohort = globals$study_labels[cohort], 
           status = ifelse(fold_reg > 1, 'activated', 'inhibited'), 
           level = paste('Collagen', level, 'vs low')) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    select(cohort, 
           subsystem,
           react_id, 
           react_name,
           level, 
           status, 
           fold_reg, 
           lower_ci, 
           upper_ci, 
           p_adjusted) %>% 
    set_names(c('Cohort', 
                'Recon subsystem', 
                'BIGG ID', 
                'Reaction name', 
                'Comparison', 
                'Activation status', 
                'Fold-regulation', 
                'lower 95% CI', 
                'upper 95% CI', 
                'pFDR')) %>% 
    mdtable(label = 'table_s9_clusters_ciochamical_reactions', 
            ref_name = 'reactions',
            caption = paste('Biochemical reactions predicted to be', 
                            'significantly activated in the collagen high', 
                            'or collagen low cluster as compared', 
                            'with collagen low cancers.', 
                            'Statistical significance was determined by', 
                            'Monte Carlo simulation.', 
                            'P values were corrected for multiple testing with', 
                            'the false discovery rate method.', 
                            'The table is available as a supplementary', 
                            'Excel file.'))
  
# Table S10: biochemical subsystems --------
  
  insert_msg('Table S10: biochemical subsystems')
  
  suppl_paper_tbl$subsystems <- meta_sub$test %>% 
    map(map, filter, status != 'regulated') %>% 
    map(map, select, -eff_size) %>% 
    map(compress, names_to = 'cohort') %>% 
    compress(names_to = 'level') %>% 
    mutate(cohort = globals$study_labels[cohort], 
           OR = signif(OR, 2), 
           p_adjusted = signif(p_adjusted, 2), 
           level = paste('Collagen', level, 'vs low')) %>% 
    select(cohort, 
           level, 
           subsystem, 
           status, 
           n, 
           n_total, 
           n_all, 
           n_all_total, 
           OR, 
           p_adjusted) %>% 
    set_names(c('Cohort', 
                'Comparison', 
                'Recon subsystem', 
                'Activation status', 
                'Regulated reaction in the subsystem', 
                'Total reactions in the subsystem', 
                'Total regulated reactions', 
                'Total reactions', 
                'Odds ratio', 
                'pFDR')) %>% 
    mdtable(label = 'table_s10_metabolic_subsystem_enrichment', 
            ref_name = 'subsystems', 
            caption = paste('Results of enrichment analysis', 
                            'for significantly activated and inhibited', 
                            'biochemical reaction within the Recon', 
                            'metabolism subsystem.', 
                            'Statistical significance was determined by', 
                            "Fisher's exact test corrected for multiple", 
                            'testing with the false discovery', 
                            'rate method (FDR).', 
                            'The table is available as a supplementary', 
                            'Excel file.'))
  
# Saving the supplementary tables ------
  
  insert_msg('Saving the supplementary tables')
  
  suppl_paper_tbl$cover <- 
    tibble(Table = paste0('Supplementary Table S', 1:length(suppl_paper_tbl)), 
           Caption = map_chr(suppl_paper_tbl, attr, 'caption'))
  
  suppl_paper_tbl <- 
    c(suppl_paper_tbl["cover"], 
      suppl_paper_tbl[names(suppl_paper_tbl) != "cover"])
  
  suppl_paper_tbl %>% 
    set_names(c('Cover', paste0('Table S', 1:(length(suppl_paper_tbl) - 1)))) %>% 
    write_xlsx(path = './report/manuscript_supplementary_tables.xlsx')
  
# END ------
  
  insert_tail()