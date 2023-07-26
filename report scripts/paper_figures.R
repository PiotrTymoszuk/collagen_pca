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
                        get_legend(norm_tumor$ribbon_plots$plots[[1]] + 
                                     theme(legend.position = 'bottom')), 
                        nrow = 2, 
                        rel_heights = c(0.7, 0.3)), 
              ncol = 2, 
              rel_widths = c(1.7, 1),
              align = 'hv', 
              axis = 'tblr')
  
  ## bottom panel: ribbon plots with the expression comparison
  ## common regulated collagen pathway genes
  
  paper_fig$norm_tumor$bottom <- norm_tumor$ribbon_plots$plots %>% 
    map(~.x + 
          labs(title = paste(.x$labels$title, .x$labels$subtitle, sep = ', ')) + 
          theme(legend.position = 'none', 
                plot.subtitle = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv',
              axis = 'tblr')
  
  ## the entire figure
  
  paper_fig$norm_tumor <- 
    plot_grid(paper_fig$norm_tumor$upper, 
              paper_fig$norm_tumor$bottom, 
              nrow = 2, 
              rel_heights = c(1, 1.4),
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_1_normal_tumor', 
              ref_name = 'norm_tumor', 
              caption = paste('Differences in tissue diffusion capacity', 
                              'and expression of the collagen pathway genes', 
                              'between the prostate cancer and benign tissue.'), 
              w = 180, 
              h = 220)
  
# Figure 2: collagen clusters ------
  
  insert_msg('Figure 2: Collagen clusters')
  
  ## left panel: heat map with the mean levels of clustering factors
  ## and cluster distribution
  
  paper_fig$coll_clusters$left <- 
    list(coll_clust$mean_hm$plot, 
         coll_clust$n_plot) %>% 
    map(~.x + 
          theme(legend.position = 'bottom', 
                plot.subtitle = element_blank(),
                axis.text = element_text(size = 7), 
                strip.text = element_text(size = 7), 
                legend.text = element_text(size = 7), 
                legend.title = element_text(size = 7))) %>%
    plot_grid(plotlist = ., 
              nrow = 2, 
              rel_heights = c(1.9, 1), 
              labels = LETTERS, 
              label_size = 10)
  
  ## legend panel
  
  paper_fig$coll_clusters$legends <- 
    cs_cluster$plots$GSE40272[c("gleason_factor", "pathology_stage_tumor")] %>% 
    map2(., c('Gleason\nscore', 'Pathologic\ntumor stage'), 
         ~.x + 
           labs(fill = .y) + 
           theme(legend.text = element_text(size = 7), 
                 legend.title = element_text(size = 7))) %>% 
    map(get_legend) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv')

  ## right panel: Gleason score and stages
  
  paper_fig$coll_clusters$right <- 
    c(map(cs_cluster$plots, 
          ~.x[["gleason_factor"]] + 
            labs(y = 'GS, % of cluster')), 
      paper_fig$coll_clusters["legends"], 
      map(cs_cluster$plots[c("GSE40272", "GSE70768", "GSE70769", "tcga")], 
          ~.x[["pathology_stage_tumor"]] + 
            labs(y = 'pT, % of cluster'))) %>% 
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
                      strip.text.y= element_text(size = 7, 
                                                angle = 0, 
                                                hjust = 0)))
  
  paper_fig$cluster_biology$left <- 
    plot_grid(paper_fig$cluster_biology$left_upper, 
              paper_fig$cluster_biology$left_bottom, 
              nrow = 2, 
              rel_heights = c(0.25, 0.75), 
              labels = LETTERS, 
              label_size = 10)
  
  ## right panel: signaling and metabolism
  
  paper_fig$cluster_biology$right <-
    list(x = c(dge_spia$cmm_plots, 
               meta_sub$regulation_plots$plots), 
         y = c('Signaling, Collagen int vs low', 
               'Signaling, Collagen hi vs low', 
               'Metabolism, Collagen int vs low', 
               'Metabolism, Collagen hi vs low'), 
         z = list(element_text(size = 7, hjust = 1, angle = 45), 
                  element_text(size = 7, hjust = 1, angle = 45), 
                  element_text(size = 7, hjust = 1, angle = 45), 
                  element_text(size = 7, hjust = 1, angle = 45))) %>% 
    pmap(function(x, y, z) x + 
           labs(title = y) + 
           scale_y_discrete(labels = function(x) map_chr(x, biol_labeller)) + 
           theme(legend.text = element_text(size = 7), 
                 legend.title = element_text(size = 7),
                 axis.text.y = element_text(size = 7), 
                 axis.text.x = z, 
                 plot.subtitle = element_blank(), 
                 legend.position = 'none', 
                 strip.background = element_blank(), 
                 strip.text = element_blank())) %>% 
    plot_grid(plotlist = ., 
              nrow = 4, 
              align = 'v', 
              axis = 'tblr', 
              rel_heights = c(0.4, 1, 0.3, 0.85), 
              labels = c('C', '', 'D', ''), 
              label_size = 10)
  
  ## the entire figure
  
  paper_fig$cluster_biology <- 
    plot_grid(paper_fig$cluster_biology$left, 
              paper_fig$cluster_biology$right, 
              ncol = 2, 
              rel_widths = c(1.25, 1)) %>% 
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
  