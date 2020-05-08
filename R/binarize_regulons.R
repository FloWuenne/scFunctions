#' Uses AUC thresholds to binarize regulon activity. Binary regulons are stored as lists of dataframes.
#'
#'
#' @param regulonAUC The AUC values for all regulons as calculated by SCENIC (content of file:3.4_regulonAUC.Rds).
#' @param thresholds The AUC thresholds as determined by auc_thresh_kmeans or another method. Format has to be a named list with regulons as names and thresholds as values!
#' @keywords SCENIC, regulons, binary activity, kmeans, thresholds
#' @import SCENIC
#' @import dplyr
#' @export
#' @examples
#' \donttest{
#' binary_regulons <- binarize_regulons(regulonAUC,regulon_thresholds)
#' }


## Binarize regulons based on thresholds
binarize_regulons <- function(regulonAUC,
                              thresholds){

  binary_regulon_list <- list()
  binary_regulon_df <- data.frame("cells" = c())

  for(regulon_no in 1:length(names(thresholds))){

    progress(regulon_no)

    regulon <- names(thresholds)[regulon_no]

    auc_df <-  data.frame("auc" = regulonAUC@assays@data@listData$AUC[regulon,],
                          "cells"= names(regulonAUC@assays@data@listData$AUC[regulon,]))

    auc_df <- auc_df %>%
      mutate("regulon"= if_else(auc >= thresholds[regulon],1,0)) %>%
      select(-auc)

    colnames(auc_df) <- c("cells",regulon)

    binary_regulon_list[[regulon]] <- auc_df
  }

  return(binary_regulon_list)
}
