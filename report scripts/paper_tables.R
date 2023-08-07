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
  
# Table S1: collagen genes of interest -----
  
  insert_msg('Table S1: collagen genes of interest')
  
  suppl_paper_tbl$genes <- globals$genes_interest %>% 
    mutate(entrez_id = mapIds(org.Hs.eg.db, 
                              keys = gene_symbol, 
                              column = 'ENTREZID', 
                              keytype = 'SYMBOL'),
           gene_group = car::recode(gene_group, "'Pro' = 'proline turnover'"), 
           gene_group = factor(gene_group, 
                               c('proline turnover', 
                                 'collagen modification', 
                                 'ECM component', 
                                 'ECM processing', 
                                 'adhesion'))) %>% 
    arrange(gene_group) %>% 
    select(gene_group, gene_symbol, entrez_id) %>% 
    set_names(c('Gene group', 'Gene symbol', 'Entrez ID')) %>% 
    mdtable(label = 'table_s2_genes', 
            ref_name = 'genes', 
            caption = paste('Investigated collagen pathway genes and', 
                            'their classification.'))
  
# Table S2: normal-tumor ------
  
  insert_msg('Table S2: normal - tumor')
  
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
    mdtable(label = 'table_s2_normal_tumor', 
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
  
# Table S3: expression of the collagen pathway genes in the clusters -------
  
  insert_msg('Table S3: collagen gene expression in the clusters')
  
  suppl_paper_tbl$cluster <- coll_clust$result_tbl %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    select(cohort, variable, 
           `Collagen low`, `Collagen hi`, 
           significance, eff_size) %>% 
    set_names(c('Cohort', 'Variable', 
                'Collagen low', 'Collagen high', 
                'Significance', 'Effect size')) %>% 
    mdtable(label = 'table_s3_cluster_collagen_genes', 
            ref_name = 'cluster', 
            caption = paste("Expression of the cluster-defining", 
                            "collagen pathway genes in the collagen clusters", 
                            "of prostate cancer.", 
                            "Statistical significance was assessed by two-tailed", 
                            "T test with Cohen's d effect size statistic.", 
                            "P values were corrected for multiple testing", 
                            "with the false discovery rate method.", 
                            "log2-transformed expression values are presented", 
                            "as medians with interquartile ranges (IQR)", 
                            "and ranges.", 
                            "The table is available as a supplementary", 
                            "Excel file."))
  
# Table S4: clinical characteristic of the collagen clusters -----
  
  insert_msg('Table S4: Clinical characteristic of the clusters')
  
  suppl_paper_tbl$clinic <- cs_cluster$result_tbl %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    filter(variable != 'Collagen Score') %>% 
    select(cohort, variable, low, hi, significance, eff_size) %>% 
    set_names(c('Cohort', 'Variable', 
                'Collagen low', 'Collagen high', 
                'Significance', 'Effect size')) %>% 
    mdtable(label = 'table_s4cluster_clinical_characteristic', 
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
             low, hi, significance, eff_size) %>% 
      set_names(c('Cohort', 'Variable', 
                  'Collagen low', 'Collagen high', 
                  'Significance', 'Effect size'))
    
  }

  suppl_paper_tbl[c("mcp", "xcell")] <- 
    suppl_paper_tbl[c("mcp", "xcell")] %>% 
    list(x = ., 
         label = c('table_s4_mcp_counter_infiltration', 
                   'table_s5_xcell_infiltration'), 
         ref_name = names(.),
         caption = map(list(paste('Non-malignant cell numbers predicted for the collagen', 
                                  'clusters by the MCP counter algorithm.'), 
                            paste('Non-maignant cell fractions predicted for the collagen', 
                                  'clusters by the xCell algorithm.')), 
                       ~paste(.x, 
                              'Statistical significance was assessed by Mann-Whitney',
                              'test with r effect size statistic.', 
                              'P values were corrected for multiple testing with the', 
                              'false discovery method.', 
                              'The table is available as a supplementary Excel file.'))) %>% 
    pmap(mdtable)
  
# Table S6: GSVA ------
  
  insert_msg('Table S6: GSVA')
  
  suppl_paper_tbl$gsva <- dge_gsva$test %>% 
    map(filter, response %in% unique(unlist(dge_gsva$cmm_signatures))) %>% 
    map2(., names(.), ~mutate(.x, cohort = globals$study_labels[.y])) %>% 
    do.call('rbind', .) %>% 
    left_join(set_names(dge_gsva$hm_plot$sign_classes, 
                        c('response', 'class')), 
              by = 'response') %>% 
    mutate(cohort = globals$study_labels[cohort], 
           response = exchange(response, 
                               reactome$lexicon),
           class = stri_replace_all(class, 
                                    fixed = '\n', 
                                    replacement = ', ')) %>% 
    select(cohort, class, response, regulation, 
           estimate, lower_ci, upper_ci, 
           p_adjusted, eff_size) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    arrange(cohort, class, response) %>% 
    set_names(c('Cohort', 
                'Signature class', 
                'Signature', 
                "Regulation, collagen high vs low", 
                'Fold-regulation', 
                'lower 95% CI', 
                'upper 95% CI', 
                'pFDR for fold-regulation', 
                'Effect size')) %>% 
    mdtable(label = 'table_s6_cluster_gsva_reactome', 
            ref_name = 'gsva', 
            caption = paste("Gene set variation", 
                            "analysis with the Reactome pathway gene signatures.", 
                            "Differences in ssGSEA scores between collagen high", 
                            "and collagen low cancers were investigated by", 
                            "two-tailed T test with Cohen's d effect", 
                            "size statistic.", 
                            "Results for signatures significantly regulated",
                            "with moderate-to-large effect size", 
                            "(d at least 0.5) in at least four", 
                            "cohorts are presented.", 
                            "P values were corrected for multiple testing with",
                            "the false discovery rate method (FDR).", 
                            "The table is available as a supplementary", 
                            "Excel file."))
  
# Table S7: differential gene expression -------
  
  insert_msg('Table S7: differential gene expression')
  
  suppl_paper_tbl$dge <- dge$signif_results %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    select(cohort, 
           gene_symbol, 
           entrez_id, 
           regulation, 
           estimate, 
           lower_ci, 
           upper_ci, 
           p_adjusted, 
           eff_size) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    set_names(c('Cohort', 
                'Gene symbol', 
                'Entrez ID', 
                'Regulation', 
                'log2 fold-regulation, collagen high vs low', 
                'lower 95% CI', 
                'upper 95% CI', 
                'pFDR', 
                'Effect size')) %>% 
    mdtable(label = 'table_s4_clusters_differential_gene_expression', 
            ref_name = 'dge', 
            caption = paste('Genes differentially expressed in the collagen', 
                            'high cluster as compared', 
                            'with collagen low cancers were identified', 
                            'by two-tailed T test', 
                            "with the 1.25-fold regulation cutoff and", 
                            "the Cohen's d effect size statistic of 0.2", 
                            'P values were corrected for multiple testing with', 
                            'the false discovery rate method (FDR).', 
                            'The table is available as a supplementary', 
                            'Excel file.'))
  
# Table S8: signaling -------
  
  insert_msg('Table S8: signaling')
  
  suppl_paper_tbl$signaling <- dge_spia$test %>% 
    map(filter, Name %in% unique(unlist(dge_spia$common))) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort], 
           pGFDr = ifelse(pGFdr < 0.05, 
                          paste('p =', signif(pGFdr, 2)), 
                          paste0('ns (p = ', signif(pGFdr, 2), ')'))) %>% 
    select(cohort, Name, ID, Status, tA, pGFdr) %>% 
    set_names(c('Cohort', 
                'Pathway name', 
                'KEGG ID', 
                'Activation status', 
                'Regulation, collagen high vs low, tA', 
                'Significance, pGFDR')) %>% 
    mdtable(label = 'table_s8_clusters_signaling_pathways', 
            ref_name = 'signaling',
            caption = paste('Signaling pathway activity in the collagen', 
                            'clusters investigated by the SPIA algorithm.', 
                            'Results for signaling pathways significantly', 
                            'activated or inhibited in at least four cohorts', 
                            'are shown.', 
                            'The table is available as a supplementary', 
                            'Excel file.'))
  
# Table S9: biochemical reactions -------
  
  insert_msg('Table S9: Biochemical reactions')
  
  suppl_paper_tbl$reactions <- meta$models %>% 
    map(components, 'regulation') %>% 
    map(filter, 
        p_adjusted < 0.05, 
        fold_reg != 1) %>% 
    map(mutate, 
        react_name = annotate_bigg(react_id, 
                                   value = 'name', 
                                   annotation_db = Recon2D), 
        react_name = unlist(react_name)) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort], 
           status = ifelse(fold_reg > 1, 'activated', 'inhibited')) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    select(cohort, 
           subsystem,
           react_id, 
           react_name,
           status, 
           fold_reg, 
           lower_ci, 
           upper_ci, 
           p_adjusted) %>% 
    set_names(c('Cohort', 
                'Recon subsystem', 
                'BIGG ID', 
                'Reaction name', 
                'Activation status, collagen high vs low', 
                'Fold-regulation, collagen high vs low', 
                'lower 95% CI', 
                'upper 95% CI', 
                'pFDR')) %>% 
    mdtable(label = 'table_s9_clusters_ciochamical_reactions', 
            ref_name = 'reactions',
            caption = paste('Biochemical reactions predicted to be', 
                            'significantly activated in the collagen high', 
                            'as compared with collagen low cancers.', 
                            'Statistical significance was determined by', 
                            'a Monte Carlo simulation.', 
                            'P values were corrected for multiple testing with', 
                            'the false discovery rate method.', 
                            'The table is available as a supplementary', 
                            'Excel file.'))
  
# Table S10: biochemical subsystems --------
  
  insert_msg('Table S10: biochemical subsystems')
  
  suppl_paper_tbl$subsystems <- meta_sub$test %>% 
    map(filter, status != 'regulated') %>% 
    map(select, -eff_size) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = globals$study_labels[cohort], 
           OR = signif(OR, 2), 
           p_adjusted = signif(p_adjusted, 2)) %>% 
    select(cohort, 
           subsystem, 
           status, 
           n, 
           n_total, 
           n_all, 
           n_all_total, 
           OR, 
           p_adjusted) %>% 
    set_names(c('Cohort', 
                'Recon subsystem', 
                'Activation status, collagen high vs low', 
                'Regulated reaction in the subsystem', 
                'Total reactions in the subsystem', 
                'Total regulated reactions, collagen high vs low', 
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