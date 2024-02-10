# supplementary figures for the transcriptome part of the manuscript

  insert_head()
  
# container ------
  
  suppl_paper_fig <- list()
  
# globals -------
  
  insert_msg('Globals')
  
  suppl_paper_fig$cohorts <- c("gse54460", "gse70768", 
                               "gse70769", "gse220095", 
                               "dkfz", "tcga")
  
# Figure S1: normal - tumor comparison ------
  
  insert_msg('Figure S1: comparison of expression, tumor versus normal')
  
  ## upper panel, Venn plots, genes regulated in both two cohort
  ## are listed
  
  suppl_paper_fig$norm_tumor$upper <- 
    plot_grid(plotlist = norm_tumor$venn_plots, 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')
  
  ## bottom panel, ribbon plots
  
  suppl_paper_fig$norm_tumor$bottom <- 
    norm_tumor$ribbon_plots %>% 
    map(~.x + 
          theme(legend.position = 'none', 
                axis.text.y = element_markdown(size = 7))) %>% 
    c(list(legend = get_legend(norm_tumor$ribbon_plots$plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              rel_widths = c(1, 1, 0.25), 
              align = 'hv', 
              axis = 'tblr')
  
  ## the entire figures
  
  suppl_paper_fig$norm_tumor <- 
    plot_grid(suppl_paper_fig$norm_tumor$upper, 
              suppl_paper_fig$norm_tumor$bottom, 
              nrow = 2, 
              align = 'hv', 
              axis = 'tblr', 
              rel_heights = c(0.25, 1), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s1_normal_tumor', 
              ref_name = 'norm_tumor', 
              caption = paste('Expression of collagen pathway genes in', 
                              'the normal prostate and prostate cancer', 
                              'tissue.'), 
              w = 180, 
              h = 230)
  
# Figure S2: development of the collagen clusters -------
  
  insert_msg('Figure S2: Development of the collagen clusters')
  
  ## upper panel: comparison of the algorithms, only CV is presented
  ## since CV performance was decisive for choice of the clustering algorithm
  
  suppl_paper_fig$clust_dev$upper <- clust_dev$plot + 
    theme(legend.position = 'none', 
          plot.subtitle = element_blank(), 
          strip.text = element_text(size = 7), 
          axis.text.y = element_markdown(size = 7)) + 
    labs(title = 'Choice of clustering algorithm, TCGA, cross-validation')
  
  suppl_paper_fig$clust_dev$upper$data <- 
    suppl_paper_fig$clust_dev$upper$data %>% 
    filter(dataset == 'cv')
  
  suppl_paper_fig$clust_dev$upper <- 
    plot_grid(suppl_paper_fig$clust_dev$upper, 
              labels = 'A', 
              label_size = 10)
  
  ## middle panel: determination of the cluster number (silhouette)
  ## and size of the clusters
  
  suppl_paper_fig$clust_dev$middle <- 
    list(clust_semi$diagnostic_plots$silhouette, 
         clust_semi$n_plot + 
           scale_y_discrete(limits = suppl_paper_fig$cohorts, 
                            labels = globals$study_labels)) %>% 
    map2(., c('Choice of cluster number, TCGA', 'Collagen cluster size'), 
         ~.x + 
           theme(plot.subtitle = element_blank(), 
                 plot.tag = element_blank(), 
                 legend.position = 'right') + 
           labs(title = .y)) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              rel_widths = c(1, 1.3), 
              labels = c('B', 'C'), 
              label_size = 10)
  
  ## bottom panel: numeric stats in the training and test cohorts

  suppl_paper_fig$clust_dev$bottom <- 
    clust_semi$stat_plots[c("sil_width", "frac_misclassified", 
                            "frac_var", "frac_np")] %>% 
    map(~.x + 
          theme(legend.position = 'none') +
          scale_y_discrete(labels = globals$study_labels))
  
  for(i in names(suppl_paper_fig$clust_dev$bottom)) {
    
    suppl_paper_fig$clust_dev$bottom[[i]]$data <- 
      suppl_paper_fig$clust_dev$bottom[[i]]$data %>% 
      filter(cohort %in% suppl_paper_fig$cohorts)
    
  }
  
  suppl_paper_fig$clust_dev$bottom <- suppl_paper_fig$clust_dev$bottom %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr', 
              labels = 'D', 
              label_size = 10)
  
  ## the entire figure
  
  suppl_paper_fig$clust_dev <- 
    plot_grid(suppl_paper_fig$clust_dev$upper, 
              suppl_paper_fig$clust_dev$middle, 
              suppl_paper_fig$clust_dev$bottom, 
              nrow = 3,
              rel_heights = c(1.3, 1, 2)) %>% 
    as_figure(label = 'figure_s2_cluster_development', 
              ref_name = 'clust_dev', 
              caption = paste('Semi-supervised clustering of prostate', 
                              'cancer samples in respect to expression of', 
                              'collagen-related genes.'), 
              w = 180, 
              h = 230)
  
# Figure S3: cluster-defining factors ------
  
  insert_msg('Figure S3: cluster defining factors')
  
  ## left panel: UMAP plots
  
  suppl_paper_fig$clusters$left <- 
    clust_semi$data_umap_plots[rev(suppl_paper_fig$cohorts)] %>% 
    map(~.x + 
          theme(legend.position = 'bottom', 
                plot.subtitle = element_blank(), 
                legend.text = element_text(size = 7))) %>% 
    plot_grid(plotlist = .,
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')
  
  ## right panel: heat map of Z score means of the collagen-related transcripts
  
  suppl_paper_fig$clusters$right <- 
    plot_grid(clust_ft$mean_hm + 
                scale_x_discrete(limits = suppl_paper_fig$cohorts, 
                                 labels = globals$study_labels) + 
                theme(legend.position = 'bottom'))
  
  ## the entire figure 
  
  suppl_paper_fig$clusters <- 
    plot_grid(plot_grid(suppl_paper_fig$clusters$left, 
                        nrow = 2, 
                        rel_heights = c(0.8, 0.2)), 
              suppl_paper_fig$clusters$right, 
              ncol = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s3_clustering_genes', 
              ref_name = 'clusters', 
              caption = paste("Expression of the collagen-related genes", 
                              "in the collagen clusters of prostate cancers."), 
              w = 180, 
              h = 230)
  
# Figure S4: clusters, ISUP and tumor stage --------
  
  insert_msg('Figure S4: clusters and tumor pathology')
  
  suppl_paper_fig$clust_tumor <- 
    ana_clinic$plot_panels[c("gleason_simple", "pt_stage")]
  
  for(i in names(suppl_paper_fig$clust_tumor)) {
    
    suppl_paper_fig$clust_tumor[[i]]$data <- 
      suppl_paper_fig$clust_tumor[[i]]$data %>% 
      filter(cohort %in% suppl_paper_fig$cohorts) %>% 
      mutate(cohort = factor(cohort, suppl_paper_fig$cohorts))
    
  }
  
  suppl_paper_fig$clust_tumor <- suppl_paper_fig$clust_tumor %>% 
    map(~.x + 
          theme(strip.text.y = element_text(angle = 0, 
                                            hjust = 0, 
                                            size = 7), 
                legend.title = element_blank())) %>% 
    plot_grid(plotlist = ., 
              nrow = 2, 
              align = 'hv', 
              axis = 'tblr', 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s4_collagen_pathology', 
              ref_name = 'clust_tumor', 
              caption = paste('ISUP and pathological tumor stage', 
                              'in the collagen clusters.'), 
              w = 160, 
              h = 180)
  
# Figure S5: clusters and infiltration ------
  
  insert_msg('Figure S5: clusters and infiltration, MCP counter')
  
  suppl_paper_fig$infiltration <- 
    list(ana_mcp$box_plots[c("Cancer associated fibroblast", "Endothelial cell")], 
         ana_xcell$box_plots[c("Cancer associated fibroblast", "Endothelial cell")]) %>% 
    transpose %>% 
    unlist(recursive = FALSE)
  
  for(i in names(suppl_paper_fig$infiltration)) {
    
    suppl_paper_fig$infiltration[[i]]$data <- 
      suppl_paper_fig$infiltration[[i]]$data %>% 
      filter(cohort %in% suppl_paper_fig$cohorts) %>% 
      mutate(cohort = factor(cohort, suppl_paper_fig$cohorts))
    
  }
  
  suppl_paper_fig$infiltration <- 
    suppl_paper_fig$infiltration %>% 
    map(~.x + 
          theme(plot.subtitle = element_blank(), 
                legend.position = 'none') + 
          scale_x_continuous(trans = 'pseudo_log')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr', 
              labels = c('A', '', 'B', ''), 
              label_size = 10) %>% 
    plot_grid(get_legend(ana_mcp$box_plots[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.9, 0.1)) %>% 
    as_figure(label = 'figure_s5_infiltration', 
              ref_name = 'intiltration', 
              caption = paste('Cancer-associated fibroblasts and endothelial', 
                              'cell infiltration in the collagen clusters.'), 
              w = 180, 
              h = 180)
  
# Figure S6: biological processes ------
  
  insert_msg('Figure S6: Biological processes')
  
  ## left panel: Reactome signatures
  
  suppl_paper_fig$biology$left <- 
    plot_grid(ana_reactome$mean_heat_map + 
                labs(title = 'Reactome pathways') + 
                scale_x_discrete(limits = suppl_paper_fig$cohorts, 
                                 labels = globals$study_labels) + 
                theme(legend.position = 'bottom', 
                      plot.subtitle = element_blank(), 
                      strip.text = element_text(size = 7)))
  
  ## right panel: GO enrichment
  
  suppl_paper_fig$biology$right <- ana_go$plots %>% 
    map2(., c('GO enrichment, Collagen hi', 
              'GO enrichment, Collagen low'), 
         ~.x + 
           labs(title = .y) + 
           theme(plot.subtitle = element_blank(), 
                 plot.tag = element_blank(), 
                 legend.position = 'bottom',
                 legend.text = element_text(size = 7))) %>% 
    plot_grid(plotlist = ., 
              nrow = 2, 
              align = 'hv', 
              axis = 'tblr')
  
  ## the final figure
  
  suppl_paper_fig$biology <- 
    plot_grid(suppl_paper_fig$biology$left, 
              suppl_paper_fig$biology$right, 
              ncol = 2, 
              rel_widths = c(1, 1.1), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s6_collagen_cluster_biology', 
              ref_name = 'biology', 
              caption = paste('Gene set variation analysis of Reactome pathway', 
                              'gene signatures and gene ontology term', 
                              'enrichment analysis for the collagen clusters.'), 
              w = 180, 
              h = 220)
  
# Figure S7: signaling and regulons ---------
  
  insert_msg('Figure S7: signaling and regulons')
  
  suppl_paper_fig$signaling <- 
    list(regulons = ana_collectri$cmm_plot, 
         signaling = ana_progeny$cmm_plot)
  
  for(i in names(suppl_paper_fig$signaling)) {
    
    suppl_paper_fig$signaling[[i]]$data <- 
      suppl_paper_fig$signaling[[i]]$data %>% 
      filter(cohort %in% suppl_paper_fig$cohorts) %>% 
      mutate(cohort = factor(cohort, suppl_paper_fig$cohorts))
    
  }
  
  suppl_paper_fig$signaling <- suppl_paper_fig$signaling %>% 
    map(~.x + 
          theme(plot.subtitle = element_blank(), 
                axis.text = element_text(size = 7)) +
          guides(x = guide_axis(angle = 90)) + 
          scale_size_area(max_size = 4.5, 
                          limits = c(0, 15)))
  
  suppl_paper_fig$signaling <- 
    plot_grid(suppl_paper_fig$signaling[['regulons']] + 
                theme(legend.position = 'none'), 
              plot_grid(suppl_paper_fig$signaling[['signaling']] + 
                          theme(legend.position = 'none'), 
                        get_legend(suppl_paper_fig$signaling[[1]] + 
                                     theme(legend.position = 'bottom', 
                                           legend.box = 'vertical') + 
                                     labs(fill = 'Regulation', 
                                          size = 'abs(LM score)')), 
                        nrow = 2, 
                        rel_heights = c(0.6, 0.4)), 
              ncol = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s7_regulons_signaling', 
              ref_name = 'signaling', 
              caption = paste('Transcriptional regulons and signaling pathway', 
                              'activity in the collagen clusters.'), 
              w = 180, 
              h = 180)
  
# Figure S8: metabolism -------
  
  insert_msg('Figure S8: metabolism')
  
  suppl_paper_fig$metabolism <- ana_metaplots$cmm_plot
  
  suppl_paper_fig$metabolism$data <- suppl_paper_fig$metabolism$data %>% 
    filter(cohort != 'gse16560') %>% 
    mutate(cohort = factor(cohort, suppl_paper_fig$cohorts))

  suppl_paper_fig$metabolism <- 
    plot_grid(suppl_paper_fig$metabolism + 
                theme(plot.subtitle = element_blank()) + 
                guides(x = guide_axis(angle = 90))) %>% 
    as_figure(label = 'figure_s8_metabolism', 
              ref_name = 'metabolism', 
              caption = paste('Metabolism in the collagen clusters.'), 
              w = 160, 
              h = 120)
  
# Figure S9: survival, algorithm comparison ---------
  
  insert_msg('Figure S9: survival, algorithm comparison')
  
  suppl_paper_fig$surv_algos <- surv_plots$algorithm_stat_plots %>% 
    map(~.x + 
          scale_x_continuous(limits = c(0.5, 1)) + 
          scale_y_continuous(limits = c(0.62, 1)) + 
          theme(legend.position = 'none')) %>% 
    c(., list(get_legend(surv_plots$algorithm_stat_plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_s9_survival_algorithm_comparison', 
              ref_name = 'surv_algos', 
              caption = paste('Comparison of performance of various machine', 
                              'learning algorithms at prediction of biochemical', 
                              'relapse-free survival with expression of the', 
                              'collagen-related genes.'), 
              w = 180, 
              h = 180)
  
# Figure S10: univariable survival: Volcano plots --------
  
  insert_msg('Figure S10: univariable survival, volcano plots')
  
  suppl_paper_fig$surv_uni_volcano <- rfs_cut$volcano_plots %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    c(., list(get_legend(rfs_cut$volcano_plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_s10_univariable_survival_volcano', 
              ref_name = 'surv_uni_volcano', 
              caption = paste('Differences in biochemical relapse-free', 
                              'survival between high and low expressors of the', 
                              'collagen-related genes.'), 
              w = 180, 
              h = 180)
  
# Figure S11: Forest plot with HR for the common survival markers --------
  
  insert_msg('Figure 11: univariable survival, common genes, Forest plot')
  
  suppl_paper_fig$surv_uni_forest <- 
    plot_grid(rfs_cut$common_forest_plot + 
                theme(legend.position = 'none')) %>% 
    as_figure(label = 'figure_s11_univariable_survival_common_markers', 
              ref_name = 'surv_uni_forest', 
              caption = paste('Collagen-related transcriptional markers of', 
                              'biochemical relapse risk shared by all', 
                              'investigated cohorts.'), 
              w = 180, 
              h = 230)
  
# Figures S12 - S14: Kaplan-Meier plots for the common survival markers ---------
  
  insert_msg('Figures S12 - S14: Kaplan-Meier plots, common survival markers')
  
  suppl_paper_fig[c('surv_uni_detail1', 
                    'surv_uni_detail2', 
                    'surv_uni_detal3')] <- 
    list(rfs_cut$km_plots[c("COL1A1", "COL3A1", "COL5A2")], 
         rfs_cut$km_plots[c("COL11A1", "COL11A2", "LOXL1", "SERPINH1")], 
         rfs_cut$km_plots[c("COL4A6", "DST", "LAMB3", "PCOLCE2")]) %>% 
    map(unlist, recursive = FALSE) %>% 
    map(map, 
        ~.x + 
          theme(legend.position = 'none',
                axis.text.x = element_text(size = 7), 
                axis.text.y = element_text(size = 7), 
                axis.title = element_text(size = 7), 
                plot.subtitle = element_text(size = 7))) %>% 
    map2(., list(c('A', '', '', 'B', '', '', 'C', '', ''), 
                 c('A', '', '', 'B', '', '', 'C', '', '', 'D', '', ''), 
                 c('A', '', '', 'B', '', '', 'C', '', '', 'D', '', '')), 
         ~plot_grid(plotlist = .x, 
                   ncol = 3, 
                   align = 'hv', 
                   labels = .y, 
                   label_size = 10)) %>% 
    map(~plot_grid(.x, 
                   get_legend(rfs_cut$km_plots[[1]][[1]] + 
                                labs(color = 'Expressors') + 
                                theme(legend.position = 'bottom')), 
                   nrow = 2, 
                   rel_heights = c(0.9, 0.1)))
  
  suppl_paper_fig[c('surv_uni_detail1', 
                    'surv_uni_detail2', 
                    'surv_uni_detal3')] <- 
    suppl_paper_fig[c('surv_uni_detail1', 
                      'surv_uni_detail2', 
                      'surv_uni_detal3')] %>% 
    list(x = ., 
         label = c('figure_s12_univariable_survival_detail', 
                   'figure_s13_univariable_survival_detail', 
                   'figure_s14_univariable_survival_detail'), 
         ref_name = names(.), 
         caption = paste('Differences in biochemical relaspse-free survival', 
                         'between high and low gene expressors:', 
                         c('_COL1A1_, _COL3A1_, and _COL5A2_.', 
                           '_COL11A1_, _COL11A2_, _LOXL1_, and _SERPINH1_.', 
                           '_COL4A6_, _DST_, _LAMB3_, and _PCOLCE2_.')), 
         h = c(175, 230, 230)) %>% 
    pmap(as_figure, 
         w = 180)
  
# saving the figures ------
  
  insert_msg('Saving the figures')
  
  suppl_paper_fig$cohorts <- NULL
  
  suppl_paper_fig <- compact(suppl_paper_fig)
  
  suppl_paper_fig %>% 
    walk(pickle, 
         path = './report/manuscript supplementary figures', 
         format = 'pdf', 
         device = cairo_pdf)
  
# END ----
  
  rm(i)
  
  insert_tail()