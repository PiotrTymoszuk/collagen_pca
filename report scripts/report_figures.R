# Main figures of the analysis report 

  insert_head()
  
# container -----
  
  report_fig <- list()
  
# Figure 1: regulation of collagen genes in the normal and tumor tissue ------

  insert_msg('Figure 1: normal - tumor gene expression')
  
  report_fig$norm_tumor$upper_panel <- 
    plot_grid(plotlist = norm_tumor$venn_plots, 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')
  
  report_fig$norm_tumor$bottom_panel <- 
    norm_tumor$volcano_plots[c('tcga', 
                               'GSE40272', 
                               'GSE70768')] %>% 
    map(~.x + 
          theme(legend.position = 'none')) %>% 
    c(list(legend = get_legend(norm_tumor$volcano_plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')
  
  report_fig$norm_tumor <- plot_grid(report_fig$norm_tumor$upper_panel, 
                                    report_fig$norm_tumor$bottom_panel, 
                                    nrow = 2, 
                                    rel_heights = c(1, 2), 
                                    labels = LETTERS, 
                                    label_size = 10) %>% 
    as_figure(label = 'figure_1_normal_tumor', 
              ref_name = 'norm_tumor', 
              caption = paste('Expression of collagen pathway genes in', 
                              'the normal prostate and prostate cancer', 
                              'tissue.'), 
              w = 180, 
              h = 210)
  
# Figure 2: collagen clusters, ribbon plots ------
  
  insert_msg('Figure 2: collagen clusters')
  
  report_fig$clusters <- coll_clust$ribbon_plots[c('tcga', 
                                                  'GSE16560', 
                                                  'GSE40272', 
                                                  'GSE70768', 
                                                  'GSE70769')] %>% 
    map2(., 
         c('TCGA, training', 
           'GSE16560, test', 
           'GSE40272, test', 
           'GSE70768, test', 
           'GSE70769, test'), 
         ~.x + 
           labs(title = .y) + 
           theme(legend.position = 'none', 
                 plot.title.position = 'plot', 
                 strip.text = element_text(size = 6), 
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
    as_figure(label = 'figure_2_clusters', 
              ref_name = 'clusters', 
              caption = paste('Clustering of prostate cancers by', 
                              'collagen pathway gene expression.'), 
              w = 180, 
              h = 225)
  
# Figure 3: collagen clusters, infiltration of CAFS and endothelial cells -----
  
  insert_msg('Figure 3: collagen clusters and infiltration')
  
  report_fig$infil$upper_panel <- 
    clust_infil$plots$mcp[c('tcga', 
                                   'GSE16560', 
                                   'GSE40272', 
                                   'GSE70768', 
                                   'GSE70769')] %>% 
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
  
  report_fig$infil$bottom_panel <- 
    clust_infil$plots$mcp[c('tcga', 
                            'GSE16560', 
                            'GSE40272', 
                            'GSE70768', 
                            'GSE70769')] %>% 
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
  
  report_fig$infil <- report_fig$infil %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 3, 
                   align = 'hv', 
                   axis = 'tblr')) %>% 
    plot_grid(plotlist = ., 
              nrow = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_3_cluster_infiltration', 
              ref_name = 'infil', 
              caption = paste('Estimates of cancer-associated fibroblast', 
                              'and endothelial cell infiltration in the ', 
                              'collagen clusters.'), 
              w = 180, 
              h = 230)
  
# Figure 4: GO enrichment, collagen clusters ------
  
  insert_msg('Figure 4: collagen clusters, GO enrichment')
  
  report_fig$gsva <- 
    plot_grid(dge_gsva$hm_plot$plot + 
                theme(legend.position = 'bottom')) %>% 
    as_figure(label = 'figure_4_clusters_gsva', 
              ref_name = 'gsva', 
              caption = paste('Common reactome pathway gene signatures', 
                              'significantly differing between the', 
                              'collagen clusters.'), 
              w = 180, 
              h = 225)

# Figure 5: collagen clusters, signaling ------
  
  insert_msg('Figure 5: collagen clusters, signaling')
  
  report_fig$spia <- dge_spia$cmm_plots %>% 
    map(~.x + 
          theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              nrow = 2, 
              rel_heights = c(1, 2.2), 
              align = 'hv', 
              labels = LETTERS, 
              label_size = 10) %>% 
    plot_grid(get_legend(dge_spia$cmm_plots[[1]] + 
                           labs(fill = 'Regulation\ntA', 
                                color = 'Regulation\ntA')), 
              ncol = 2, 
              rel_widths = c(0.8, 0.2)) %>% 
    as_figure(label = 'figure_5_signaling', 
              ref_name = 'spia', 
              caption = 'Common regulated signaling pathways in Collagen intermediate and Collagen high tumors.', 
              w = 180, 
              h = 150)
  
# Figure 6: development of the survival score ------
  
  insert_msg('Figure 6: development of the collagen score')
  
  report_fig$coll_score <- 
    plot_grid(cs_plots$c_ibs_plot + 
                theme(legend.position = 'none'), 
              plot_grid(get_legend(cs_plots$c_ibs_plot + 
                                     theme(legend.position = 'bottom', 
                                           legend.justification = c(0, 0))), 
                        get_legend(cs_plots$coef_plot + 
                                     theme(legend.position = 'bottom', 
                                           legend.justification = c(0, 0))), 
                        nrow = 4, 
                        align = 'hv', 
                        axis = 'tblr'), 
              nrow = 2) %>% 
    plot_grid(cs_plots$coef_plot + 
                expand_limits(x = 0.75) + 
                labs(subtitle = paste0(stri_replace(cs_plots$n_tags$tcga, 
                                                    fixed = '\n', 
                                                    replacement = ''), 
                                       ', \u03BB = ', 
                                       signif(coll_score$opt_lambda$lambda, 2))) + 
                theme(legend.position = 'none', 
                      plot.tag = element_blank()), 
              ., 
              ncol = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_6_collagen_score', 
              ref_name = 'coll_score', 
              caption = paste('Development of the Collagen Score for', 
                              'prediction of prostate cancer survival.'), 
              w = 180, 
              h = 160)

# Figure 7: collagen score and survival ------
  
  insert_msg('Figure 7: collagen score and survival')
  
  report_fig$coll_score_km <- 
    cs_plots$cox_cal_plots[c('tcga', 
                             'GSE16560', 
                             'GSE40272', 
                             'GSE70768', 
                             'GSE70769')] %>% 
    map2(., 
         c('TCGA, training', 
           'GSE16560, test', 
           'GSE40272, test', 
           'GSE70768, test', 
           'GSE70769, test'), 
         ~.x + 
           labs(title = .y, 
                color = 'Collagen\nScore')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_7_collagen_score_km', 
              ref_name = 'coll_score_km', 
              caption = 'Prediction of prostate cancer survival by the Collagen Score.', 
              w = 180, 
              h = 230)
  
# Saving the figures ------
  
  insert_msg('Saving the figures')
  
  report_fig %>% 
    walk(pickle, 
         path = './report/figures', 
         format = 'pdf', 
         device = cairo_pdf)
  
# END ------
  
  insert_tail()