#' Calculate CSI module activity over all cell types
#'
#' @param clusters_df -
#' @param regulonAUC -
#' @param metadata -
#' @keywords SCENIC, regulons, CSI activity
#' @export
#' @examples
#'

calc_csi_module_activity <- function(clusters_df,regulonAUC,metadata,cell_type_column){

  require(tidyverse)
  require(pheatmap)
  require(viridis)
  
  metadata$cell_type <- metadata[ , cell_type_column ]
  cell_types<- unique(metadata$cell_type)
  regulons <- unique(clusters_df$regulon)

  regulonAUC_sub <- regulonAUC@assays@data@listData$AUC
  regulonAUC_sub <- regulonAUC_sub[regulons,]

  csi_activity_matrix_list <- list()
  csi_cluster_activity <- data.frame("csi_cluster" = c(),
                                     "mean_activity" = c(),
                                     "cell_type" = c())

  cell_type_counter <- 0
  regulon_counter <-
    for(ct in cell_types) {
      cell_type_counter <- cell_type_counter + 1

      cell_type_aucs <- rowMeans(regulonAUC_sub[,rownames(subset(metadata,cell_type == ct))])
      cell_type_aucs_df <- data.frame("regulon" = names(cell_type_aucs),
                                      "activtiy"= cell_type_aucs,
                                      "cell_type" = ct)
      csi_activity_matrix_list[[ct]] <- cell_type_aucs_df
    }

  for(ct in names(csi_activity_matrix_list)){
    for(cluster in unique(clusters_df$csi_cluster)){
      csi_regulon <- subset(clusters_df,csi_cluster == cluster)

      csi_regulon_activtiy <- subset(csi_activity_matrix_list[[ct]],regulon %in% csi_regulon$regulon)
      csi_activtiy_mean <- mean(csi_regulon_activtiy$activtiy)
      this_cluster_ct_activity <- data.frame("csi_cluster" = cluster,
                                             "mean_activity" = csi_activtiy_mean,
                                             "cell_type" = ct)
      csi_cluster_activity <- rbind(csi_cluster_activity,this_cluster_ct_activity)
    }
  }

  csi_cluster_activity[is.na(csi_cluster_activity)] <- 0

  csi_cluster_activity_wide <- csi_cluster_activity %>%
    spread(cell_type,mean_activity)

  rownames(csi_cluster_activity_wide) <- csi_cluster_activity_wide$csi_cluster
  csi_cluster_activity_wide <- as.matrix(csi_cluster_activity_wide[2:ncol(csi_cluster_activity_wide)])

  return(csi_cluster_activity_wide)
}
