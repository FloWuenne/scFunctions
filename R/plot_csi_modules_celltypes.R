#' Plots a heatmap for the connection specificty indices for all regulons.
#'
#' @param csi_df Data frame containing CSI values for all pairwise regulons.
#' @keywords SCENIC, regulons, CSI
#' @export
#' @examples
#'

plot_csi_modules_celltypes <- function(csi_df){

  require(tidyverse)
  require(pheatmap)
  require(viridis)

  ## subset csi data frame based on threshold
  csi_test_mat <- csi_df %>%
    unique() %>%
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
           # cutree_cols = 10,
           # cutree_rows = 10,
           fontsize_row = 12,
           treeheight_row = 0,
           treeheight_col = 40,
           widt = 2000,
           height = 1800)
  #dev.off()

}
