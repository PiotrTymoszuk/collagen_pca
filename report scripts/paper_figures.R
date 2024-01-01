# Figures for the manuscript. The GSE16560 cohort is excluded from all figures
# since the cohort does not fulfill the selection criteria for the studies (BCR
# data availability)

  insert_head()
  
# container ------
  
  paper_fig <- list()
  
# Figure 1: MRI, tissue and Gleason score comparison ------
  
  insert_msg('Figure 1: MRI, tissue and Gleason')
  
  ## upper panel: MRI results
  
  paper_fig$norm_tumor$upper <- 
    plot_grid(ggdraw(clip = 'on') + 
                draw_image('./aux files/prostate_mri.png') + 
                theme(plot.margin = globals$common_margin), 
              plot_grid(mri$plot + 
                          theme(legend.position = 'none'), 
                        nrow = 2, 
                        rel_heights = c(0.9, 0.1)), 
              ncol = 2, 
              rel_widths = c(1.7, 1),
              align = 'hv', 
              axis = 'tblr')
  
  ## bottom panel: tissues and Gleason
  
  paper_fig$norm_tumor$bottom_plots <- 
    list(tissue = norm_tumor$hm_plot, 
         gleason = gs_uni$mean_heat_map) %>% 
    map2(., c('Tumor vs benign', 'Gleason score'), 
         ~.x + 
           labs(title = .y) + 
           scale_fill_gradient2(low = 'steelblue', 
                                mid = 'black', 
                                high = 'firebrick', 
                                midpoint = 0, 
                                limits = c(-1, 1), 
                                oob = scales::squish, 
                                name = 'Mean Z-score') + 
           theme(axis.line = element_blank()))
  
  paper_fig$norm_tumor$bottom_plots$gleason <- 
    paper_fig$norm_tumor$bottom_plots$gleason + 
    scale_x_discrete(limits = c("gse54460", "gse70768", 
                                "gse70769", "gse220095", 
                                "dkfz", "tcga"), 
                     labels = globals$study_labels)
  
  paper_fig$norm_tumor$bottom <- 
    plot_grid(paper_fig$norm_tumor$bottom_plots$tissue + 
                theme(legend.position = 'none', 
                      plot.subtitle = element_blank()), 
              plot_grid(paper_fig$norm_tumor$bottom_plots$gleason + 
                          theme(legend.position = 'none', 
                                plot.subtitle = element_blank()), 
                        get_legend(paper_fig$norm_tumor$bottom_plots[[1]] + 
                                     theme(legend.position = 'bottom')), 
                        nrow = 2, 
                        rel_heights = c(1, 1)), 
              ncol = 2, 
              rel_widths = c(1, 1.3), 
              labels = c('B', 'C'), 
              label_size = 10)
  
  ## the entire figure
  
  paper_fig$norm_tumor <- 
    plot_grid(paper_fig$norm_tumor$upper, 
              paper_fig$norm_tumor$bottom, 
              nrow = 2, 
              rel_heights = c(1, 1.3),
              labels = c('A', ''), 
              label_size = 10) %>% 
    as_figure(label = 'figure_1_normal_tumor', 
              ref_name = 'norm_tumor', 
              caption = paste('Differences in tissue diffusion capacity', 
                              'and expression of collagen-related genes', 
                              'between the prostate cancer and benign tissue.', 
                              'Expression of collagen-related genes in cancers', 
                              'stratified by Gleason scores.'), 
              w = 180, 
              h = 180)
  
# Figure 2: infiltration and cellular localization of collagen genes -------
  
  insert_msg('Figure 2: collagen genes and tumor microenvironment composition')
  
  ## upper panel: MCP counter, EC and CAFs
  
  paper_fig$tme$upper <- 
    ana_mcp$box_plots[c("Cancer associated fibroblast", "Endothelial cell")] %>% 
    map2(., 
         paste0(c('Fibroblasts', 'EC'), ', MCP Counter'), 
         ~.x + 
          scale_y_discrete(limits = c("gse54460", "gse70768", 
                                      "gse70769", "gse220095", 
                                      "dkfz", "tcga"), 
                           labels = globals$study_labels) + 
          theme(legend.position = 'none', 
                axis.title.y = element_text(size = 8, angle = 90), 
                axis.title.x = element_blank()) + 
          coord_flip() + 
          guides(x = guide_axis(angle = 45)) + 
          labs(x = 'Z-score') + 
          scale_x_continuous(trans = 'pseudo_log')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    plot_grid(get_legend(ana_mcp$box_plots[[1]] + 
                           scale_fill_manual(values = globals$cluster_colors, 
                                             labels = c('low', 'high'), 
                                             name = 'Collagen') + 
                           scale_color_manual(values = globals$cluster_colors, 
                                              labels = c('low', 'high'), 
                                              name = 'Collagen')),
              ncol = 2, 
              rel_widths = c(0.9, 0.1))
  
  ## bottom panel: scRNA results, placeholder
  
  paper_fig$tme$bottom <- ggdraw()
  
  ## the entire figure
  
  paper_fig$tme <- plot_grid(paper_fig$tme$upper, 
                             paper_fig$tme$bottom, 
                             nrow = 2, 
                             rel_heights = c(1, 3), 
                             labels = LETTERS, 
                             label_size = 10) %>% 
    as_figure(label = 'figure_2_tumor_microenvironment', 
              ref_name = 'tme', 
              caption = paste('Expression of collagen-related genes', 
                              'in components of the tumor microenvironment.'), 
              w = 180, 
              h = 210)
  
# Figure 3: transcriptional collagen score and RFS -------
  
  insert_msg('Figure 3: transcriptional collagen score')
  
  ## upper panel: development and performance of the transcriptional
  ## collagen score
  
  paper_fig$coll_score$upper <- 
    plot_grid(cs_plots$coef_plot + 
                theme(plot.subtitle = element_blank(), 
                      legend.position = 'none', 
                      axis.text.y = element_text(size = 7)) + 
                labs(title = 'Transcriptional Collagen Score, development'), 
              plot_grid(cs_plots$stat_plot + 
                          theme(plot.subtitle = element_blank(), 
                                legend.position = 'none') + 
                          labs(title = 'Transcriptional Collagen Score, RFS prediction'), 
                        get_legend(cs_plots$stat_plot + 
                                     theme(legend.position = 'bottom')), 
                        ggdraw(), 
                        nrow = 3, 
                        rel_heights = c(1, 0.12, 0.12)), 
              ncol = 2, 
              labels = LETTERS, 
              label_size = 10)
  
  ## bottom panel: Kaplan-Meier plots of RFS in tertiles of the Collagen Score, 
  ## the GSE70768, DKFZ and TCGA cohorts, where the score performed the best
  
  paper_fig$coll_score$bottom <- 
    cs_plots$cox_cal_plots[c("tcga", "gse70768", "dkfz")] %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    c(list(get_legend(cs_plots$cox_cal_plots[[1]] + 
                        scale_color_manual(values = c('darkolivegreen4', 
                                                      'steelblue4', 
                                                      'firebrick4'), 
                                           labels = c('low', 'int', 'high'), 
                                           name = 'Transcriptional Collagen Score') + 
                        theme(legend.position = 'right')))) %>% 
    plot_grid(plotlist = .,
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')
  
  ## the entire figure
  
  paper_fig$coll_score <- 
    plot_grid(paper_fig$coll_score$upper, 
              paper_fig$coll_score$bottom, 
              nrow = 2, 
              rel_heights = c(1, 1.2), 
              labels = c('', 'C'), 
              label_size = 10) %>% 
    as_figure(label = 'figure_3_transcriptional_collagen_score', 
              ref_name = 'coll_score', 
              caption = paste('Transcriptional Collagen Score and prediction',
                              'of biochemical relapse-free survival.'),
              w = 180, 
              h = 210)
  
# Saving the figures ------
  
  insert_msg('Saving the figures')
  
  paper_fig %>% 
    walk(pickle, 
         path = './report/manuscript figures', 
         format = 'pdf', 
         device = cairo_pdf)
  
# END ------
  
  insert_tail()
  