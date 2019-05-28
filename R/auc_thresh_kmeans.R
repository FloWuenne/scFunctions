#' This function takes the regulon AUC values and determines thresholds for binarizatio based on kmeans clustering.
#' Less specific than the SCENIC function AUCell_exploreThresholds but much faster.
#'
#'
#' @param regulonAUC The AUC values for all regulons as calculated by SCENIC (content of file:3.4_regulonAUC.Rds).
#' @keywords SCENIC, regulons, binary activity, kmeans, thresholds
#' @export
#' @examples
#' regulon_thresholds <- auc_thresh_kmeans(regulonAUC)


auc_thresh_kmeans <- function(regulonAUC){

  require(SCENIC)
  require(svMisc)
  require(dplyr)
  require(tidyr)

  ## Iterate over each regulon in the AUC matrix
  regulons <- rownames(regulonAUC@assays$data@listData$AUC)

  kmeans_thresholds <- list()

  print("Processing regulon distributions...")

  for(regulon_no in 1:length(regulons)){

    progress(regulon_no)

    regulon <- regulons[regulon_no]

    df <- data.frame("auc" = regulonAUC@assays$data@listData$AUC[regulon,],
                     "cells"= names(regulonAUC@assays$data@listData$AUC[regulon,]),
                     "regulon" = regulon)

    ## Remove cells with 0 AUC as they interfere with kmeans clustering
    df <- df %>%
      subset(auc > 0)

    km <- kmeans(df$auc,centers=2)
    df$cluster <- as.factor(km$cluster)

    cluster1_max <- max(subset(df,cluster == 1)$auc)
    cluster2_max <- max(subset(df,cluster == 2)$auc)

    if(cluster1_max > cluster2_max){

      df <- df %>%
        mutate("cluster" = gsub(2,3,cluster)) %>%
        mutate("cluster" = gsub(1,2,cluster)) %>%
        mutate("cluster" = gsub(3,1,cluster))
    }

    df <- df %>%
      arrange(desc(auc))

    df_sub <- df %>%
      subset(cluster == 1)

    auc_thresholds <- df_sub[1,]$auc

    kmeans_thresholds[[regulon]] <- auc_thresholds

  }


  print("Done evaluating thresholds...")
  return(kmeans_thresholds)

}
