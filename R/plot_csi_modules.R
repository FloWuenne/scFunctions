#' Plots a heatmap for the connection specificty indices for all regulons.
#'
#' @param csi_df Data frame containing CSI values for all pairwise regulons.
#' @param nclust Number of clusters to divide the heatmap into
#' @param font_size_regulons Font size for regulon names.
#' @keywords SCENIC, regulons, CSI
#' @import tidyverse
#' @import pheatmap
#' @import viridis
#' @export
#' @examples
#'

plot_csi_modules <- function(csi_df,
                             nclust = 10,
                             font_size_regulons = 6){

  ## subset csi data frame based on threshold
  csi_test_mat <- csi_df %>%
    spread(regulon_2,CSI)

  future_rownames <- csi_test_mat$regulon_1
  csi_test_mat <- as.matrix(csi_test_mat[,2:ncol(csi_test_mat)])
  rownames(csi_test_mat) <- future_rownames
  # png("./csi_heatmap.png",
  #     width = 2400,
  #     height = 1800)
  pheatmap(csi_test_mat,
           show_colnames = FALSE,
           color = viridis(n = 10),
           cutree_cols = 10,
           cutree_rows = 10,
           fontsize_row = font_size_regulons,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           treeheight_row = nclust,
           treeheight_col = nclust,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           widt = 2000,
           height = 3200)
  #dev.off()

}
