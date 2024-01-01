# Definition of collagen clusters based on expression of the 55 collagen-related 
# genes.
#
# 1) Tuning of the clustering structure in the TCGA cohort. We're trying out 
# several common algorithms (Ward's hierarchical clustering, KMEANs and PAM), 
# as well as a novel regularized KMEANS algorithm along with several distance 
# metrics. The cluster number choice is motivated by the peak silhouette 
# statistic as well as the elbow of the within-cluster sum of squares curve. 
# The optimal clustering algorithm is chosen with silhouette width, explained 
# clustering variance, neighborhood preservation and low percentage of 
# observations with negative silhouette widths (potential of misclassification). 
# The performance stats are computed for the training data set and in 5-fold 
# cross validation. Pre-processing: normalization with mean centering.
#
# 2) Semi-supervised clustering. The collagen clusters trained in the TCGA cohort
# are projected onto the remaining cohorts with use of an inverse-weighted 
# k-nearest neighbor classifier. The quality of predictions is assessed with 
# numeric statistics listed above as well as visual and numeric analysis of cross
# distances between the clusters.
#
# 3) Comparison of levels of the cluster-defining factors between the collagen 
# clusters (two cluster solution, two-tailed T test with Cohen's d effect size 
# statisic). 
#
# 4) Permutation importance of clustering factors in the training TCGA cohort.

# tools -------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)

  library(clustTools)
  library(exda)
  library(microViz)
  
  library(furrr)
  library(soucer)

  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  set_rownames <- trafo::set_rownames
  extract <- clustTools::extract
  var <- clustTools::var
  
  c('./tools/globals.R', 
    './tools/functions.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# analysis scripts --------
  
  insert_msg('Analysis scripts')
  
  ## analysis globals
  
  source_all('./clustering scripts/globals.R', 
             message = TRUE, crash = TRUE)
  
  ## cluster tuning and semi-supervised clustering
  
  list(cache_path = c('./cache/clust_dev.RData', 
                      './cache/clust_semi.RData', 
                      './cache/clust_imp.RData'), 
       script_path = c('./clustering scripts/clust_devel.R', 
                       './clustering scripts/semi_clustering.R', 
                       './clustering scripts/importance.R'), 
       message = c('Loading cached cluster tuning results', 
                   'Cached results of semi-supervised clustering', 
                   'Cached permutation importance testing')) %>% 
    pwalk(access_cache)
  
  ## clustering features
  
  c('./clustering scripts/cluster_features.R') %>% 
    source_all(message = TRUE, crash = TRUE)

# END --------
  
  insert_tail()