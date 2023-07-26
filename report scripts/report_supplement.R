# supplementary figures

  insert_head()
  
# container ------
  
  suppl_report_fig <- list()
  
# Figure S1: genes upregulated in prostate cancer ------
  
  insert_msg('Figure S1: genes upregulated in PC')
  
  suppl_report_fig$norm_tumor_up <- norm_tumor$plots[c('tcga', 
                                                       'GSE40272', 
                                                       'GSE70768')] %>% 
    map(~.x[norm_tumor$common$upregulated]) %>% 
    transpose %>% 
    unlist(recursive = FALSE) %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              labels = c('A', '', '', 
                         'B', '', '', 
                         'C', '', '', 
                         'D', '', ''), 
              label_size = 10) %>% 
    as_figure(label = 'figure_s1_normal_tumor_up', 
              ref_name = 'norm_tumor_up', 
              caption = paste('Expression of ALDH18A1, PYCR1, P4HB and PPIB', 
                              'in the normal prostate and prostate', 
                              'cancer tissue.'), 
              w = 180, 
              h = 230)
  
# Figure S2: genes downregulated in prostate cancer -----
  
  insert_msg('Figure S2: genes downregulated in prostate cancer')
  
  suppl_report_fig$norm_tumor_down <- norm_tumor$plots[c('tcga', 
                                                         'GSE40272', 
                                                         'GSE70768')] %>% 
    map(~.x[norm_tumor$common$downregulated]) %>% 
    transpose %>% 
    unlist(recursive = FALSE) %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              labels = c('A', '', '', 
                         'B', '', ''), 
              label_size = 10) %>% 
    as_figure(label = 'figure_s2_normal_tumor_down', 
              ref_name = 'norm_tumor_down', 
              caption = paste('Expression of LAMB3 and COL4A2 in',
                              'the normal prostate and prostate', 
                              'cancer tissue.'), 
              w = 180, 
              h = 120)
  
# Figure S3 - S5: correlation of collagen gene expression ------- 
  
  insert_msg('Figures S3 - S5: correlation of the collagen pathway genes')
  
  suppl_report_fig[c('corr1', 
                     'corr2', 
                     'corr3')] <- list(corr$bubble_plots[c('tcga', 'GSE16560')], 
                                       corr$bubble_plots[c('GSE40272', 'GSE70768')], 
                                       corr$bubble_plots['GSE70769']) %>% 
    map(~map(.x, ~.x + theme(legend.position = 'none'))) %>% 
    map2(., c(2, 2, 1), 
         ~plot_grid(plotlist = .x, 
                    nrow = .y, 
                    align = 'hv', 
                    labels = if(.y > 1) LETTERS else NULL, 
                    label_size = 10) %>% 
           plot_grid(get_legend(corr$bubble_plots[[1]]), 
                     ncol = 2, 
                     rel_widths = c(0.9, 0.1)))
  
  suppl_report_fig[c('corr1', 
                     'corr2', 
                     'corr3')] <- suppl_report_fig[c('corr1', 
                                                     'corr2', 
                                                     'corr3')] %>% 
    list(x = ., 
         label = c('figure_s3_correlation', 
                   'figure_s4_correlation', 
                   'figure_s5_correlation'), 
         ref_name = c('corr1', 
                      'corr2', 
                      'corr3'), 
         caption = paste('Pairwise correlation of collagen pathway gene expression', 
                         c('TCGA and GSE16560 cohorts.', 
                           'GSE40272 and GSE70768 cohorts.', 
                           'GSE70769 cohort.'), 
                         sep = ', '), 
         h = c(230, 230, 135)) %>% 
    pmap(as_figure, 
         w = 180)
  
# Figure S6: development of the collagen clusters -------
  
  insert_msg('Figure S6: development of the collagen clusters')
  
  suppl_report_fig$clust_dev$upper_panel <- 
    plot_grid(clust_dev$plot + 
                theme(legend.position = 'none'), 
              plot_grid(get_legend(clust_dev$plot), 
                        ggdraw() + 
                          draw_text(coll_clust$diagnostic_plots$wss$labels$tag, 
                                    size = 8, 
                                    x = 0.1, 
                                    hjust = 0), 
                        nrow = 2), 
              ncol = 2, 
              rel_widths = c(0.7, 0.3), 
              labels = c('A', ''), 
              label_size = 10)
  
  suppl_report_fig$clust_dev$bottom_panel <- 
    c(coll_clust$diagnostic_plots[c('wss', 'silhouette')], 
      coll_clust[c("n_plot", "variance_plot")]) %>% 
    list(x = ., 
         y = c('none', 'none', 'bottom', 'bottom'), 
         z = c(8, 8, 7, 8)) %>% 
    pmap(function(x, y, z) x + 
           theme(legend.position = y, 
                 plot.subtitle = element_blank(), 
                 plot.tag = element_blank(), 
                 axis.text.x = element_text(size = z))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'v', 
              axis = 'tblr', 
              rel_heights = c(1, 1.15), 
              labels = c('B', '', 'C'), 
              label_size = 10)
  
  suppl_report_fig$clust_dev <- 
    plot_grid(suppl_report_fig$clust_dev$upper_panel, 
              suppl_report_fig$clust_dev$bottom_panel, 
              nrow = 2,
              rel_heights = c(1.25, 2)) %>% 
    as_figure(label = 'figure_s6_cluster_development', 
              ref_name = 'clust_dev', 
              caption = paste('Semi-supervised clustering of prostate', 
                              'cancer samples in respect to expression of', 
                              'the collagen pathway genes.'), 
              w = 180, 
              h = 230)
  
# Figure S7: heat maps of the collagen clusters -------
  
  insert_msg('Figure S7: collagen clusters heat map')
  
  suppl_report_fig$clust_hm <- 
    coll_clust$hm_plots[c('tcga', 
                          'GSE16560', 
                          'GSE40272', 
                          'GSE70768', 
                          'GSE70769')] %>% 
    map(~.x + 
          labs(subtitle = .x$labels$tag) + 
          theme(legend.position = 'none', 
                plot.tag = element_blank(), 
                axis.text.y = element_markdown(size = 6))) %>% 
    c(list(legend = get_legend(coll_clust$hm_plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_s7_heat_maps_clusters', 
              ref_name = 'clust_hm', 
              caption = paste('Expression levels of collagen pathway genes', 
                              'in the collagen clusters.'), 
              w = 180, 
              h = 230)
  
# Figure S8: clusters and Gleason score ------
  
  insert_msg('Figure S8: Gleason score and clusters')
  
  suppl_report_fig$clust_gleason <- cs_cluster$plots[c('tcga', 
                                                       'GSE16560', 
                                                       'GSE40272', 
                                                       'GSE70768', 
                                                       'GSE70769')] %>% 
    map(~.x$gleason + 
          theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    as_figure(label = 'figure_s8_cluster_gleason', 
              ref_name = 'clust_gleason', 
              caption = paste('Distribution of Gleason scores in', 
                              'the collagen clusters.'), 
              w = 180, 
              h = 210)
  
# Figure S9: clusters and the MCP counter infiltration ------
  
  insert_msg('Figure S9: clusters and the MCP counter infiltration estimates')
  
  suppl_report_fig$mcp_infil <- 
    clust_infil$ribbon_panels$mcp[c("tcga", 
                                    "GSE16560", 
                                    "GSE40272", 
                                    "GSE70768", 
                                    "GSE70769")] %>% 
    map2(., 
         globals$study_labels[c("tcga", 
                                "GSE16560", 
                                "GSE40272", 
                                "GSE70768", 
                                "GSE70769")], 
         ~.x + 
           labs(title = .y) + 
           theme(legend.position = 'none', 
                 plot.title.position = 'panel', 
                 axis.text = element_markdown(size = 6))) %>% 
    c(list(legend = get_legend(clust_infil$ribbon_panels$mcp[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_s9_clusters_mcp', 
              ref_name = 'mcp_infil', 
              caption = paste('MCP Counter estimates of non-malignant', 
                              'cell content in the collagen clusters.'), 
              w = 180, 
              h = 230)
  
# Figure S10: collagen clusters, infiltration of CAFS and endothelial cells, xCell -----
  
  insert_msg('Figure S10: collagen clusters and infiltration, xCell')
  
  ## upper panel: fibroblasts
  
  suppl_report_fig$xcell_infil$upper_panel <- 
    clust_infil$plots$xcell[c("tcga", 
                              "GSE16560", 
                              "GSE40272", 
                              "GSE70768", 
                              "GSE70769")] %>% 
    map2(., 
         paste0('CAF, ', 
                globals$study_labels[c('tcga', 
                                       'GSE16560', 
                                       'GSE40272', 
                                       'GSE70768', 
                                       'GSE70769')]), 
         ~.x[['Cancer associated fibroblast']] + 
           labs(title = .y) + 
           theme(legend.position = 'none'))
  
  ## bottom panel: endothelial cells
  
  suppl_report_fig$xcell_infil$bttom_panel <- 
    clust_infil$plots$xcell[c("tcga", 
                              "GSE16560", 
                              "GSE40272", 
                              "GSE70768", 
                              "GSE70769")] %>% 
    map2(., 
         paste0('Endothelial cells, ', 
                globals$study_labels[c('tcga', 
                                       'GSE16560', 
                                       'GSE40272', 
                                       'GSE70768', 
                                       'GSE70769')]), 
         ~.x[['Endothelial cell']] + 
           labs(title = .y) + 
           theme(legend.position = 'none'))
  
  ## the entire figure
  
  suppl_report_fig$xcell_infil <- 
    suppl_report_fig$xcell_infil %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 3, 
                   align = 'hv')) %>% 
    plot_grid(plotlist = ., 
              nrow = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s10_clusters_xcell', 
              ref_name = 'infil', 
              caption = paste('xCell estimates of cancer-associated', 
                              'fibroblast and enbothelial cell', 
                              'infiltration in the collagen clusters.'), 
              w = 180, 
              h = 230)

# Figure S11 - S13: differential gene expression -----
  
  insert_msg('Figures S11 - S13: Differential gene expression')
  
  suppl_report_fig[c('dge1', 'dge2', 'dge3')] <- 
    list(map(dge_plots$volcano_plots[c("dge_collagen_int", "dge_collagen_hi")], 
             ~.x[c('tcga', 'GSE16560')]), 
         map(dge_plots$volcano_plots[c("dge_collagen_int", "dge_collagen_hi")], 
             ~.x[c('GSE40272', 'GSE70768')]), 
         map(dge_plots$volcano_plots[c("dge_collagen_int", "dge_collagen_hi")], 
             ~.x[c('GSE70769')])) %>% 
    map(transpose) %>% 
    map(unlist, recursive = FALSE) %>% 
    map(map, ~.x + theme(legend.position = 'none')) %>% 
    map2(., c(2, 2, 1), 
         ~plot_grid(plotlist = .x, 
                    ncol = 2, 
                    nrow = .y, 
                    align = 'hv', 
                    axis = 'tblr', 
                    labels = if(.y > 1) c('A', '', 'B') else NULL, 
                    label_size = 10)) %>% 
    map(~plot_grid(.x, 
                   get_legend(dge_plots$volcano_plots$dge_collagen_hi$tcga + 
                                theme(legend.position = 'bottom')), 
                   nrow = 2, 
                   rel_heights = c(0.9, 0.1)))
  
  suppl_report_fig[c('dge1', 'dge2', 'dge3')] <-  
    suppl_report_fig[c('dge1', 'dge2', 'dge3')] %>% 
    list(x = ., 
         label = c('figure_s11_dge', 
                   'figure_s12_dge', 
                   'figure_s13_dge'), 
         ref_name = c('dge1', 
                      'dge2', 
                      'dge3'), 
         caption = paste('Genes differentially expressed in the collagen clusters', 
                         c('TCGA and GSE16560 cohorts.', 
                           'GSE40272 and GSE70768 cohorts.', 
                           'GSE70769 cohort.'), 
                         sep = ', '), 
         h = c(190, 190, 100)) %>% 
    pmap(as_figure, 
         w = 190)
  
# Figure S14 - S15: GO enrichment -------
  
  insert_msg('Figure S14 - s15: GO enrichment')
  
  suppl_report_fig[c('go1', 'go2')] <- 
    list(map(dge_go$top_go$plots[c("int", "high")], 
             ~.x[c('tcga', 'GSE16560', 'GSE40272')]), 
         map(dge_go$top_go$plots[c("int", "high")], 
             ~.x[c('GSE70768', 'GSE70769')])) %>% 
    map(transpose) %>% 
    map(unlist, recursive = FALSE) %>% 
    map(map, 
        ~.x + 
          theme(plot.title.position = 'plot', 
                plot.title = element_text(hjust = 0.5), 
                axis.text.y = element_text(size = 7))) %>% 
    map2(., c(3, 2), 
         ~plot_grid(plotlist = .x, 
                    nrow = .y, 
                    ncol = 2, 
                    align = 'hv', 
                    labels = if(.y == 3) c('A', '', 'B', '', 'C') else c('A', '', 'B'), 
                    label_size = 10))
  
  suppl_report_fig[c('go1', 'go2')] <- 
    suppl_report_fig[c('go1','go2')] %>% 
    list(x = ., 
         label = c('figure_s14_go', 
                   'figure_s15_go'), 
         ref_name = c('go1', 'go2'), 
         caption = paste('Top most enriched biological process GO terms in the collagen clusters', 
                         c('TCGA, GSE16560 and GSE40272 cohorts.', 
                           'GSE70768 and GSE70769 cohorts.'), 
                         sep = ', '), 
         h = c(220, 150)) %>% 
    pmap(as_figure, 
         w = 180)
  
# Figure S16: regulation of ECM and focal adhesion pathways, schemes -----
  
  insert_msg('Figure S16 Regulation of focal adhesion pathway components')
  
  suppl_report_fig$pathview1 <- 
    plot_grid(plot_grid(ggdraw()  +
                          draw_text('Focal adhesion pathway, Collagen int vs low tumors', 
                                    size = 8, 
                                    fontface = 'bold'), 
                        ggdraw()  +
                          draw_image('./report/kegg pathviews int/hsa04510.pathview.png'), 
                        nrow = 2, 
                        ncol = 1, 
                        rel_heights = c(0.1, 0.9)), 
              plot_grid(ggdraw()  +
                          draw_text('Focal adhesion pathway, Collagen high vs low tumors', 
                                    size = 8, 
                                    fontface = 'bold'), 
                        ggdraw()  +
                          draw_image('./report/kegg pathviews hi/hsa04510.pathview.png'), 
                        nrow = 2, 
                        ncol = 1, 
                        rel_heights = c(0.1, 0.9)), 
              nrow = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s16_focal_adhesion_view', 
              ref_name = 'pathview1', 
              caption = paste('Components of the focal adhesion pathway', 
                              'regulated in the collagen clusters.'), 
              w = 180, 
              h = 220)
  
# Figure S17: regulation of the ECM pathway ------ 
  
  insert_msg('Figure S17: regulation of the ECM pathway')
  
  suppl_report_fig$pathview2 <- 
    plot_grid(plot_grid(ggdraw()  +
                          draw_text('ECM-receptor interaction pathway, Collagen int vs low tumors', 
                                    size = 8, 
                                    fontface = 'bold'), 
                        ggdraw()  +
                          draw_image('./report/kegg pathviews int/hsa04512.pathview.png'), 
                        nrow = 2, 
                        ncol = 1, 
                        rel_heights = c(0.1, 0.9)), 
              plot_grid(ggdraw()  +
                          draw_text('ECM-receptor interaction pathway, Collagen high vs low tumors', 
                                    size = 8, 
                                    fontface = 'bold'), 
                        ggdraw()  +
                          draw_image('./report/kegg pathviews hi/hsa04512.pathview.png'), 
                        nrow = 2, 
                        ncol = 1, 
                        rel_heights = c(0.1, 0.9)), 
              nrow = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s17_ecm_view', 
              ref_name = 'pathview2', 
              caption = paste('Components of the ECM receptor interaction', 
                              'pathway regulated in the collagen clusters.'), 
              w = 180, 
              h = 220)
  
# Figure S18: regulation of the actin cytoskeleton pathway -----
  
  insert_msg('Figure S18: actin cytoskeleton pathway')
  
  suppl_report_fig$pathview3 <- 
    plot_grid(ggdraw()  +
                draw_text('Actin cytoskeleton regulation pathway, Collagen high vs low tumors', 
                          size = 8, 
                          fontface = 'bold'), 
              ggdraw()  +
                draw_image('./report/kegg pathviews hi/hsa04810.pathview.png'), 
              nrow = 2, 
              ncol = 1, 
              rel_heights = c(0.1, 0.9)) %>%  
    as_figure(label = 'figure_s18_actin_view', 
              ref_name = 'pathview3', 
              caption = paste('Components of the actin cytoskeleton', 
                              'pathway regulated in the collagen clusters.'), 
              w = 180, 
              h = 120)
  
# Figure S19: components of the small cell lung cancer pathway ----
  
  insert_msg('Figure S19: SCLC pathway')
  
  suppl_report_fig$pathview4 <- 
    plot_grid(ggdraw()  +
                draw_text('Small cell lung cancer pathway, Collagen high vs low tumors', 
                          size = 8, 
                          fontface = 'bold'), 
              ggdraw()  +
                draw_image('./report/kegg pathviews hi/hsa05222.pathview.png'), 
              nrow = 2, 
              ncol = 1, 
              rel_heights = c(0.1, 0.9)) %>% 
    as_figure(label = 'figure_s19_small_cell_lc_view', 
              ref_name = 'pathview4', 
              caption = paste('Components of the small cell lung cancer', 
                              'pathway pathway regulated in the', 
                              'collagen clusters.'), 
              w = 180, 
              h = 120)
  
# Figure S20: components of the pathways in cancer ------
  
  insert_msg('Figure 20: components of the pathways in cancer')
  
  suppl_report_fig$pathview5 <- 
    plot_grid(ggdraw()  +
                draw_text('Pathways in cancer, Collagen high vs low tumors', 
                          size = 8, 
                          fontface = 'bold'), 
              ggdraw()  +
                draw_image('./report/kegg pathviews hi/hsa05200.pathview.png'), 
              nrow = 2, 
              ncol = 1, 
              rel_heights = c(0.1, 0.9)) %>% 
    as_figure(label = 'figure_s20_pathways_cancer_view', 
              ref_name = 'pathview4', 
              caption = paste('Components of pathways in cancer', 
                              'regulated in the', 
                              'collagen clusters.'), 
              w = 180, 
              h = 120)
  
# Figure S21 - S23: signaling pathway regulation, volcano plots ------
  
  insert_msg('Figures S21 - S23: volcano plots of SPIA results')
  
  suppl_report_fig[c('spia1', 'spia2', 'spia3')] <- 
    list(map(dge_spia$volcano_plots[c("int", "hi")], 
             ~.x[c('tcga', 'GSE16560')]), 
         map(dge_spia$volcano_plots[c("int", "hi")], 
             ~.x[c('GSE40272', 'GSE70768')]), 
         map(dge_spia$volcano_plots[c("int", "hi")], 
             ~.x[c('GSE70769')])) %>% 
    map(transpose) %>% 
    map(unlist, recursive = FALSE) %>% 
    map(map, ~.x + theme(legend.position = 'none')) %>% 
    map2(., c(2, 2, 1), 
         ~plot_grid(plotlist = .x, 
                    nrow = .y, 
                    ncol = 2, 
                    align = 'hv', 
                    axis = 'tblr', 
                    labels = if(.y == 2) c('A', '', 'B') else NULL, 
                    label_size = 10)) %>% 
    map(~plot_grid(.x, 
                   get_legend(dge_spia$volcano_plots$hi$tcga + 
                                theme(legend.position = 'bottom')), 
                   nrow = 2, 
                   rel_heights = c(0.9, 0.1)))
  
  suppl_report_fig[c('spia1', 'spia2', 'spia3')] <- 
    suppl_report_fig[c('spia1', 'spia2', 'spia3')] %>% 
    list(x = ., 
         label = c('figure_s21_pathways', 
                   'figure_s22_pathways', 
                   'figure_s23_pathways'), 
         ref_name = c('spia1', 'spia2', 'spia3'), 
         caption = paste('Modulation of signaling pathways in the collagen clusters', 
                         c('TCGA and GSE16560 cohorts.', 
                           'GSE40272 and GSE70768 cohorts.', 
                           'GSE70769 cohort.'), 
                         sep = ', '), 
         h = c(200, 200, 110)) %>% 
    pmap(as_figure, 
         w = 180)
  
# Figure S24: biochemical reactions, subsystem enrichment ------
  
  insert_msg('Figure S24: metabolic subsystem enrichment')
  
  suppl_report_fig$react_enrich <- meta_sub$regulation_plots$plots %>% 
    map(~.x + 
          scale_radius(limits = c(0, 6.4), 
                       range = c(0.5, 5)))
  
  suppl_report_fig$react_enrich <- 
    suppl_report_fig$react_enrich %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              nrow = 2, 
              align = 'hv', 
              axis = 'tblr', 
              rel_heights = c(0.3, 1), 
              labels = LETTERS, 
              label_size = 10) %>% 
    plot_grid(., 
              get_legend(suppl_report_fig$react_enrich$hi), 
              ncol = 2,
              rel_widths = c(0.85, 0.15)) %>% 
    as_figure(label = 'figure_s24_metabolic_subsystems', 
              ref_name = 'react_enrich', 
              caption = paste('Metabolic subsystems predicted to be regulated in', 
                              'the collagen clusters.'), 
              h = 150, 
              w = 180)

# Figure S25: FAOx reaction regulation -------
  
  insert_msg('Figure S25: FAOx reaction regulation')
  
  suppl_report_fig$faox <- 
    meta_hg$regulation_plots[c("int", "hi")] %>% 
    map(~.x[['Fatty acid oxidation']]) %>% 
    map(map, ~.x + theme(legend.position = 'none')) %>% 
    map(~c(.x, 
           list(get_legend(meta_hg$regulation_plots$hi$`Fatty acid oxidation`$tcga)))) %>% 
    unlist(recursive = FALSE) %>% 
    map(~.x + theme(plot.title.position = 'plot')) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr', 
              labels = c('A', '', '', 
                         '', '', '', 
                         'B', '', '', 
                         '', '', ''), 
              label_size = 10) %>% 
    as_figure(label = 'figure_s25_faox_regulation', 
              ref_name = 'faox', 
              caption = paste('Predicted regulation of the fatty acid oxidation', 
                              'reactions in the Collagen clusters.'), 
              w = 180, 
              h = 220)
    
# Figure S26: TCA regulation --------  
  
  insert_msg('Figure S26: regulation of TCA reactions')
  
  suppl_report_fig$tca <- 
    meta_hg$regulation_plots[c("int", "hi")] %>% 
    map(~.x[['Citric acid cycle']]) %>% 
    map(map, ~.x + theme(legend.position = 'none')) %>% 
    map(~c(.x, 
           list(get_legend(meta_hg$regulation_plots$hi$`Citric acid cycle`$tcga)))) %>% 
    unlist(recursive = FALSE) %>% 
    map(~.x + theme(plot.title.position = 'plot')) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr', 
              labels = c('A', '', '', 
                         '', '', '', 
                         'B', '', '', 
                         '', '', ''), 
              label_size = 10) %>% 
    as_figure(label = 'figure_s26_tca_regulation', 
              ref_name = 'faox', 
              caption = paste('Predicted regulation of the citric acid cycle', 
                              'reactions in the Collagen clusters.'), 
              w = 180, 
              h = 220)

# Figure S27: OxPhos -------
  
  insert_msg('Figure S27: oxidative phosphorylation')
  
  suppl_report_fig$oxphos <- 
    meta_hg$regulation_plots[c("int", "hi")] %>% 
    map(~.x[['Oxidative phosphorylation']]) %>% 
    map(map, ~.x + theme(legend.position = 'none')) %>% 
    map(~c(.x, 
           list(get_legend(meta_hg$regulation_plots$hi$`Oxidative phosphorylation`$tcga)))) %>% 
    unlist(recursive = FALSE) %>% 
    map(~.x + theme(plot.title.position = 'plot')) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr', 
              labels = c('A', '', '', 
                         '', '', '', 
                         'B', '', '', 
                         '', '', ''), 
              label_size = 10) %>% 
    as_figure(label = 'figure_s27_oxphos_regulation', 
              ref_name = 'faox', 
              caption = paste('Predicted regulation of the oxidative', 
                              'phosphorylation reactions in the', 
                              'Collagen clusters.'), 
              w = 180, 
              h = 220)

# Figure S28: collagen score, Gleason score and the cluster assignment ------
  
  insert_msg('Figure S28: collagen score, Gleason score and the cluster assignment')
  
  ## upper panel: Gleason and Collagen score
  
  suppl_report_fig$cs_gleason$upper_panel <- 
    coll_clinic$gleason_plots[c('tcga', 
                                'GSE16560', 
                                'GSE40272', 
                                'GSE70768', 
                                'GSE70769')] %>% 
    map(~.x$gleason_collagen_score) %>% 
    map(~.x + 
          theme(plot.title.position = 'plot', 
                plot.subtitle = element_text(size = 7)) + 
          labs(subtitle = stri_replace(.x$labels$subtitle, 
                                       fixed = 'rho', 
                                       replacement = '\u03C1')))
  
  ## bottom panel: Gleason score and collagen clusters
  
  suppl_report_fig$cs_gleason$bottom_panel <- 
    cs_cluster$plots[c('tcga', 
                       'GSE16560', 
                       'GSE40272', 
                       'GSE70768', 
                       'GSE70769')] %>% 
    map(~.x$collagen_score) %>% 
    map(~.x + theme(legend.position = 'none'))
  
  ## the entire figure
  
  suppl_report_fig$cs_gleason <- 
    suppl_report_fig$cs_gleason %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 3, 
                   align = 'hv', 
                   axis = 'tblr')) %>% 
    plot_grid(plotlist = ., 
              nrow = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s28_cs_gleason_clusters', 
              ref_name = 'cs_gleason', 
              caption = paste('Association of the collagen score with', 
                              'the Gleason score and collagen cluster', 
                              'assignment.'), 
              w = 180, 
              h = 220)
  
# Figure S29: collagen score as an independent prognostic factor -----
  
  insert_msg('Figure S29: Collagen score as an independent prognostic factor')
  
  suppl_report_fig$cs_multi <- plot_grid(surv_multi$c_index_plot) %>% 
    as_figure(label = 'figure_s29_cs_independent', 
              ref_name = 'cs_multi', 
              caption = paste('Prediction of prostate cancer survival', 
                              'by the collagen score with Gleason score', 
                              'and age confounder variables.'), 
              w = 160, 
              h = 120)
  
# saving the figures ------
  
  insert_msg('Saving the supplementary figures')
  
  suppl_report_fig %>% 
    walk(pickle, 
         path = './report/report supplementary figures', 
         format = 'pdf', 
         device = cairo_pdf)
  
# END ------
  
  insert_tail()