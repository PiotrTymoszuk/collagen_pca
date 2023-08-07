# supplementary figures for the transcriptome part of the manuscript

  insert_head()
  
# container ------
  
  suppl_paper_fig <- list()
  
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
    norm_tumor$ribbon_plots$plots %>% 
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
  
  ## upper panel: comparison of the algorithms
  
  suppl_paper_fig$clust_dev$upper <- 
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
  
  ## bottom panels: determination of the cluster number
  ## and quality of semi-supervised clustering
  
  suppl_paper_fig$clust_dev$bottom <- 
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
  
  ## the entire figure
  
  suppl_paper_fig$clust_dev <- 
    plot_grid(suppl_paper_fig$clust_dev$upper, 
              suppl_paper_fig$clust_dev$bottom, 
              nrow = 2,
              rel_heights = c(1.25, 2)) %>% 
    as_figure(label = 'figure_s2_cluster_development', 
              ref_name = 'clust_dev', 
              caption = paste('Semi-supervised clustering of prostate', 
                              'cancer samples in respect to expression of', 
                              'the collagen pathway genes.'), 
              w = 180, 
              h = 230)
  
# Figure S3: cluster-defining factors ------
  
  insert_msg('Figure S3: cluster defining factors')
  
  suppl_paper_fig$clusters <- 
    coll_clust$ribbon_plots[c("tcga", 
                              "GSE16560", 
                              "GSE70768", 
                              "GSE70769", 
                              "GSE116918")] %>% 
    map2(., 
         c('TCGA, training', 
           'GSE16560, test', 
           'GSE70768, test', 
           'GSE70769, test', 
           'GSE116918, test'), 
         ~.x + 
           labs(title = .y) + 
           theme(legend.position = 'none', 
                 plot.title.position = 'plot', 
                 strip.text = element_text(size = 6), 
                 axis.text.y = element_markdown(size = 6),
                 plot.subtitle = element_text(size = 7), 
                 axis.text.x = element_text(size = 7), 
                 axis.title.x = element_text(size = 7), 
                 plot.margin = ggplot2::margin(t = 2, 
                                               r = 2, 
                                               b = 2, 
                                               l = 2, 
                                               unit = 'mm'))) %>% 
    c(list(legend = get_legend(coll_clust$ribbon_plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_s3_clusters', 
              ref_name = 'clusters', 
              caption = paste("Expression of the collagen pathway genes", 
                              "in the collagen clusters of prostate cancers."), 
              w = 180, 
              h = 230)

# Figure S4: clusters and infiltration, MCP counter ------
  
  insert_msg('Figure S4: clusters and infiltration, MCP counter')
  
  suppl_paper_fig$mcp_infil <- 
    list(x = clust_infil$ribbon_panels$mcp, 
         y = globals$study_labels, 
         z = clust_infil$test$mcp) %>% 
    pmap(function(x, y, z) x + 
           labs(title = y) + 
           scale_y_discrete(limits = rev(c('Cancer associated fibroblast', 
                                           'Endothelial cell', 
                                           'Macrophage/Monocyte', 
                                           'Monocyte', 
                                           'Neutrophil', 
                                           'Myeloid dendritic cell', 
                                           'B cell', 
                                           'NK cell', 
                                           'T cell', 
                                           'T cell CD8+')),
                            labels = set_names(z[['plot_lab']], 
                                               z[['variable']])) + 
           theme(legend.position = 'none', 
                 plot.title.position = 'panel', 
                 axis.text = element_markdown(size = 6), 
                 axis.title.x = element_text(size = 7), 
                 plot.subtitle = element_text(size = 7))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_s4_clusters_mcp', 
              ref_name = 'mcp_infil', 
              caption = paste('MCP Counter estimates of non-malignant', 
                              'cell content in the collagen clusters.'), 
              w = 180, 
              h = 230)
  
# Figure S5: CAF and EC, MCP counter, details ------
  
  insert_msg('Figure S5: MCP counter CAF and EC')
  
  ## upper panel: fibroblasts
  
  suppl_paper_fig$infil_mcp$upper <- 
    clust_infil$plots$mcp %>% 
    map2(., 
         paste0('CAF, ', globals$study_labels), 
         ~.x[['Cancer associated fibroblast']] + 
           scale_y_continuous(trans = 'sqrt') + 
           labs(title = .y, 
                y = 'cells, MCP counter') + 
           theme(legend.position = 'none'))
  
  ## bottom panel: EC
  
  suppl_paper_fig$infil_mcp$bottom <- 
    clust_infil$plots$mcp %>% 
    map2(., 
         paste0('Endothelial cells, ', globals$study_labels), 
         ~.x[['Endothelial cell']] + 
           scale_y_continuous(trans = 'sqrt') + 
           labs(title = .y, 
                y = 'cells, MCP counter') + 
           theme(legend.position = 'none'))
  
  ## the entire figure
  
  suppl_paper_fig$infil_mcp <- suppl_paper_fig$infil_mcp %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 3, 
                   align = 'hv', 
                   axis = 'tblr')) %>% 
    plot_grid(plotlist = ., 
              nrow = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s5_caf_ec_mcpcounter', 
              ref_name = 'infil_mcp', 
              caption = paste('Counts of cancer-associated fibroblast', 
                              'and endothelial cell in tumor samples in the ', 
                              'collagen clusters predicted by', 
                              'the MCP counter algorithm.'), 
              w = 180, 
              h = 230)
  
# Figure S6: CAF and EC, xCell counter, details ------
  
  insert_msg('Figure S6: xCell, CAF and EC')
  
  ## upper panel: fibroblasts
  
  suppl_paper_fig$infil_xcell$upper <- 
    clust_infil$plots$xcell %>% 
    map2(., 
         paste0('CAF, ', globals$study_labels), 
         ~.x[['Cancer associated fibroblast']] + 
           scale_y_continuous(trans = 'sqrt') + 
           labs(title = .y, 
                y = 'cells, MCP counter') + 
           theme(legend.position = 'none'))
  
  ## bottom panel: EC
  
  suppl_paper_fig$infil_xcell$bottom <- 
    clust_infil$plots$xcell %>% 
    map2(., 
         paste0('Endothelial cells, ', globals$study_labels), 
         ~.x[['Endothelial cell']] + 
           scale_y_continuous(trans = 'sqrt') + 
           labs(title = .y, 
                y = 'cells, MCP counter') + 
           theme(legend.position = 'none'))
  
  ## the entire figure
  
  suppl_paper_fig$infil_xcell <- suppl_paper_fig$infil_xcell %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 3, 
                   align = 'hv', 
                   axis = 'tblr')) %>% 
    plot_grid(plotlist = ., 
              nrow = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s6_caf_ec_xcell', 
              ref_name = 'infil_xcell', 
              caption = paste('Fractions of cancer-associated fibroblast', 
                              'and endothelial cell in tumor samples in the ', 
                              'collagen clusters predicted by', 
                              'the xCell algorithm.'), 
              w = 180, 
              h = 230)
  
  
  
# Figure S7: representative Reactome signatures ------
  
  insert_msg('Figure S7: representative signatures')
  
  ## top five signatures from the common ones for each signature
  ## category (e.g. ECM/collagen, GF signaling etc.)
  
  suppl_paper_fig$top_gsva <- dge_gsva$hm_plot$plot
  
  suppl_paper_fig$top_gsva$data <- suppl_paper_fig$top_gsva$data %>% 
    filter(variable %in% dge_gsva$top_signatures$top_variables$variable)
  
  suppl_paper_fig$top_gsva <- 
    plot_grid(suppl_paper_fig$top_gsva + 
                scale_y_discrete(labels = function(x) map_chr(x, gsva_labeller)) + 
                theme(legend.position = 'bottom', 
                      plot.subtitle = element_blank())) %>% 
    as_figure(label = 'figure_s7_gsva_reactome_representative', 
              ref_name = 'top_gsva', 
              caption = paste('Gene set variation analysis results', 
                              'for representative Reactome pathway gene', 
                              'signatures differentiating between the collagen', 
                              'clusters.'), 
              w = 180, 
              h = 210)
  
# Figure S8: differential gene expression, numbers -------
  
  insert_msg('Figure S8: Differential gene expression, gene numbers')

  suppl_paper_fig$dge_numbers <- dge_plots$dge_numbers$plot %>% 
    as_figure(label = 'figure_s8_diff_gene_expression_numbers', 
              ref_name = 'dge_numbers', 
              caption = paste('Percentages of the analyzed transcriptome', 
                              'significantly up- and downregulated in',
                              'collagen high vs collagen low cluster cancers.'), 
              w = 150, 
              h = 90)
  
# Figure S9: top differentially expressed genes -----
  
  insert_msg('Figure S9: top differentially expressed genes')
  
  suppl_paper_fig$top_dge <- dge_plots$top_plots %>% 
    map(~.x + 
          labs(subtitle = .x$labels$subtitle %>% 
                 stri_replace_all(regex = '(Collagen|clusters:)\\s{1}', 
                                  replacement = '')) + 
          theme(legend.position = 'none', 
                plot.title.position = 'plot', 
                axis.text.y = element_text(size = 7))) %>% 
    c(list(get_legend(dge_plots$top_plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv',
              axis = 'tblr') %>% 
    as_figure(label = 'figure_s9_top_diff_gene_expression', 
              ref_name = 'top_dge', 
              caption = paste('Top strongest differentially expressed genes in', 
                              'the collagen high cluster as compared with', 
                              'collagen low cancers.'), 
              w = 180, 
              h = 230)

# Figure S10: numbers of differentially regulated reactions ------
  
  insert_msg('Figure S10: numbers of differentially regulated reactions')
  
  suppl_paper_fig$reaction_counts <- meta_plots$counts$plot %>% 
    as_figure(label = 'figure_s10_clusters_reaction_numbers', 
              ref_name = 'reaction_counts', 
              caption = paste('Numbers of significantly activated and', 
                              'inhibited biochemical reactions predicted', 
                              'for the collagen high cluster', 
                              'as compared with collagen low cancers.'), 
              w = 150, 
              h = 90)

# Figure S11: regulation of fatty acid oxidation and TCA reactions ------
  
  insert_msg('Figure S11: Regulation of FAOX and TCA reactions')
  
  suppl_paper_fig$faox <- c('Fatty acid oxidation', 
                            'Citric acid cycle') %>% 
    map(~meta_hg$regulation_plots[[.x]]) %>% 
    list(react = ., 
         label = c('FAOX reactions', 'TCA reactions'), 
         span = list(c(-1.1, 1.5), 
                  c(-0.82, 0.15))) %>% 
    pmap(function(react, label, span) react %>% 
            map(~.x + 
                  labs(title = .x$labels$title %>% 
                         stri_replace(regex = '.*,\\s{1}', 
                                      replacement = ''), 
                       x = .x$labels$x %>% 
                         stri_replace(fixed = 'Reactions', 
                                      replacement = label)) + 
                  theme(legend.position = 'none') + 
                  scale_y_continuous(limits = span))) %>% 
    map(~c(.x, list(get_legend(meta_hg$regulation_plots[[1]][[1]])))) %>% 
    unlist(recursive = FALSE) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr', 
              labels = c('A', '', '', 
                         '', '', '', 
                         'B', '', '', 
                         '', '', ''), 
              label_size = 10) %>% 
    as_figure(label = 'figure_s11_clusters_faox_and_tca_activity', 
              ref_name = 'faox', 
              caption = paste('Activity of fatty acid oxidation and', 
                              'citric acid cycle reactions', 
                              'predicted for the collagen high cluster', 
                              'as compared with collagen low cancers.'), 
              w = 180,
              h = 230)
  
# Figure S12: regulation of extracellular transport reactions ------
  
  insert_msg('Figure S12: regulation of extracellular transport reactions')
  
  suppl_paper_fig$extrans <- c('Chondroitin synthesis', 
                               'Transport, extracellular') %>% 
    map(~meta_hg$regulation_plots[[.x]]) %>% 
    list(react = ., 
         label = c('CS reactions', 'ExTrans reactions'), 
         span = list(c(-0.4, 1.5), 
                     c(-1.75, 2.5))) %>% 
    pmap(function(react, label, span) react %>% 
           map(~.x + 
                 labs(title = .x$labels$title %>% 
                        stri_replace(regex = '.*,\\s{1}', 
                                     replacement = ''), 
                      x = .x$labels$x %>% 
                        stri_replace(fixed = 'Reactions', 
                                     replacement = label)) + 
                 theme(legend.position = 'none') + 
                 scale_y_continuous(limits = span))) %>% 
    map(~c(.x, list(get_legend(meta_hg$regulation_plots[[1]][[1]])))) %>% 
    unlist(recursive = FALSE) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr', 
              labels = c('A', '', '', 
                         '', '', '', 
                         'B', '', '', 
                         '', '', ''), 
              label_size = 10) %>% 
    as_figure(label = 'figure_s12_clusters_chondroitin_and_transport_activity', 
              ref_name = 'extrans', 
              caption = paste('Activity of chondroitin synthesis and', 
                              'extracellular transport reactions', 
                              'predicted for the collagen high cluster', 
                              'as compared with collagen low cancers.'), 
              w = 180,
              h = 230)
  
# saving the figures ------
  
  insert_msg('Saving the figures')
  
  suppl_paper_fig %>% 
    walk(pickle, 
         path = './report/manuscript supplementary figures', 
         format = 'pdf', 
         device = cairo_pdf)
  
# END ----
  
  insert_tail()