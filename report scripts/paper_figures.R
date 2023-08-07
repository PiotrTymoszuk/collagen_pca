# Figures for the manuscript

  insert_head()
  
# container ------
  
  paper_fig <- list()
  
# Figure 1: MRI and tissue comparison ------
  
  insert_msg('Figure 1: MRI and tissue comparison')
  
  ## upper panel: MRI results
  
  paper_fig$norm_tumor$upper <- 
    plot_grid(ggdraw(clip = 'on') + 
                draw_image('./aux files/prostate_mri.png') + 
                theme(plot.margin = globals$common_margin), 
              plot_grid(mri$plot + 
                          theme(legend.position = 'none'), 
                        nrow = 2, 
                        rel_heights = c(0.7, 0.3)), 
              ncol = 2, 
              rel_widths = c(1.7, 1),
              align = 'hv', 
              axis = 'tblr')
  
  ## bottom panel: heat map of mean normalized expression 
  ## for the common regulated collagen pathway genes
  
  paper_fig$norm_tumor$bottom <- norm_tumor$hm_plot$plot

  ## the entire figure
  
  paper_fig$norm_tumor <- 
    plot_grid(paper_fig$norm_tumor$upper, 
              paper_fig$norm_tumor$bottom, 
              nrow = 2, 
              rel_heights = c(1.4, 1),
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_1_normal_tumor', 
              ref_name = 'norm_tumor', 
              caption = paste('Differences in tissue diffusion capacity', 
                              'and expression of the collagen pathway genes', 
                              'between the prostate cancer and benign tissue.'), 
              w = 180, 
              h = 160)
  
# Figure 2: collagen clusters ------
  
  insert_msg('Figure 2: Collagen clusters')
  
  ## left panel: heat map with the mean levels of clustering factors
  ## and cluster distribution
  
  paper_fig$coll_clusters$left <- 
    list(coll_clust$umap_layouts$tcga + 
           labs(x = 'UMAP 1', 
                y = 'UMAP 2'), 
         coll_clust$mean_hm$plot) %>% 
    map2(., c('right', 'bottom'), 
         ~.x + 
          theme(legend.position = .y, 
                plot.subtitle = element_blank(),
                axis.text = element_text(size = 7), 
                strip.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))) %>%
    plot_grid(plotlist = ., 
              nrow = 2, 
              rel_heights = c(1, 2.8), 
              labels = LETTERS, 
              label_size = 10)
  
  ## clinical plots 
  
  paper_fig$coll_clusters$gleason <- cs_cluster$plots %>% 
    map(~.x$gleason_factor) %>% 
    map2(., names(.), 
         ~.x + 
           labs(title = globals$study_labels[.y], 
                y = 'GS, % of cluster') + 
           scale_fill_manual(values = c('5 - 6' = 'bisque', 
                                        '7' = 'coral1', 
                                        '8+' = 'coral4'), 
                             name = 'Gleason\nscore', 
                             drop = FALSE))
  
  paper_fig$coll_clusters$patho <- 
    cs_cluster$plots[c("GSE70768", "GSE70769", "GSE116918", "tcga")] %>% 
    map(~.x$pathology_stage_tumor) %>% 
    map2(., names(.), 
         ~.x + 
           labs(title = globals$study_labels[.y], 
                y = 'pT, % of cluster') + 
           scale_fill_manual(values = c(T1 = 'bisque', 
                                        T2 = 'coral1', 
                                        T3 = 'coral4', 
                                        T4 = 'gray60'), 
                             name = 'Pathologic\ntumor stage', 
                             drop = FALSE))
  
  ## legend panel
  
  paper_fig$coll_clusters$legends <- 
    list(paper_fig$coll_clusters$gleason$GSE116918, 
         paper_fig$coll_clusters$patho$GSE116918) %>% 
    map(~.x + 
          theme(legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))) %>% 
    map(get_legend) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv')

  ## right panel: Gleason score and stages
  
  paper_fig$coll_clusters$right <- 
    c(paper_fig$coll_clusters$gleason, 
      paper_fig$coll_clusters["legends"], 
      paper_fig$coll_clusters$patho) %>% 
    map(~.x + 
          labs(title = stri_replace(.x$labels$title, 
                                    regex = '.*,', 
                                    replacement = '')) + 
          theme(legend.position = 'none', 
                plot.subtitle = element_text(size = 7), 
                axis.text = element_text(size = 7),
                axis.title.x = element_blank(), 
                axis.title.y = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr', 
              labels = c('C', '', 
                         '', '', 
                         '', '', 
                         'D', ''), 
              label_size = 10)
  
  ## the entire figure
  
  paper_fig$coll_clusters <- 
    plot_grid(paper_fig$coll_clusters$left, 
              paper_fig$coll_clusters$right, 
              ncol = 2, 
              rel_widths = c(1, 1.1)) %>% 
    as_figure(label = 'figure_2_collagen_clusters_pcarcinoma', 
              ref_name = 'coll_clusters', 
              caption = paste('Collagen clusters of prostate carcinoma', 
                              'and their clinical characteristic.'), 
              w = 180,
              h = 230)
  
  
  
# Figure 3: collagen cluster biology ------
  
  insert_msg('Figure 3: collagen cluster biology')
  
  ## left panel: infiltration and GSVA
  
  paper_fig$cluster_biology$left_upper <- 
    clust_infil$plots$mcp$tcga[c("Cancer associated fibroblast", "Endothelial cell")] %>% 
    map2(., c('CAF, TCGA', 'EC, TCGA'), 
         ~.x + 
           labs(title = .y) + 
           scale_y_continuous(trans = 'sqrt', 
                              labels = function(x) x/10000) + 
           labs(y = '10<sup>4</sup> \u00D7cells, MCP counter') + 
           theme(legend.position = 'none', 
                 axis.text = element_text(size = 7), 
                 axis.title.y = element_markdown(size = 7), 
                 plot.subtitle = element_text(size = 7))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')
  
  paper_fig$cluster_biology$left_bottom <- 
    plot_grid(dge_gsva$hm_plot$plot + 
                labs(fill = 'mean\nssGSEA') + 
                theme(legend.position = 'bottom', 
                      axis.text.x = element_text(size = 7), 
                      plot.subtitle = element_blank(), 
                      axis.text.y = element_blank(),
                      axis.line.y = element_blank(), 
                      axis.ticks.y = element_blank(), 
                      strip.text.x = element_text(size = 7), 
                      strip.text.y = element_text(size = 7, 
                                                angle = 0, 
                                                hjust = 0)))
  
  paper_fig$cluster_biology$left <- 
    plot_grid(paper_fig$cluster_biology$left_upper, 
              paper_fig$cluster_biology$left_bottom, 
              nrow = 2, 
              rel_heights = c(0.27, 1), 
              labels = LETTERS, 
              label_size = 10)
  
  ## right panel: signaling and metabolism
  
  paper_fig$cluster_biology$right <- 
    list(x = list(dge_spia$cmm_plot, 
                  meta_sub$regulation_plots$plot), 
         y = c('Signaling, collagen hi vs low', 
               'Metabolism, collagen hi vs low'), 
         v = c('KEGG pathways, SPIA', 
               'Recon metabolis subsystem'), 
         z = list(element_text(size = 7, hjust = 1, angle = 45), 
                  element_text(size = 7, hjust = 1, angle = 45))) %>% 
    pmap(function(x, y, z, v) x + 
           labs(title = y, 
                subtitle = v) + 
           scale_y_discrete(labels = function(x) map_chr(x, biol_labeller)) + 
           theme(legend.text = element_text(size = 7), 
                 legend.title = element_text(size = 7),
                 axis.text.y = element_text(size = 7), 
                 axis.text.x = z, 
                 legend.position = 'none', 
                 strip.background = element_blank(), 
                 strip.text = element_blank())) %>% 
    plot_grid(plotlist = ., 
              nrow = 2, 
              align = 'v', 
              axis = 'tblr', 
              rel_heights = c(1, 1), 
              labels = c('C', 'D'), 
              label_size = 10)

  ## the entire figure
  
  paper_fig$cluster_biology <- 
    plot_grid(paper_fig$cluster_biology$left, 
              paper_fig$cluster_biology$right, 
              ncol = 2, 
              rel_widths = c(1.2, 1)) %>% 
    as_figure(label = 'figure_3_collagen_cluster_biology', 
              ref_name = 'cluster_biology', 
              caption = paste('Infiltration, biological processes, signaling', 
                              'and metabolic pathways differentially', 
                              'regulated in the collagen clusters.'), 
              w = 180, 
              h = 230)
    
# Saving the figures ------
  
  insert_msg('Saving the figures')
  
  paper_fig %>% 
    walk(pickle, 
         path = './report/manuscript figures', 
         format = 'pdf', 
         device = cairo_pdf)
  
# END ------
  
  insert_tail()
  