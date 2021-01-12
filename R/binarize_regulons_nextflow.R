#' This function is to analyze the nextflow implementation results of SCENIC!
#' Uses AUC thresholds to binarize regulon activity. Binary regulons are stored as lists of dataframes.
#'
#' @param regulonsAUC The AUC values for all regulons as calculated by SCENICprotocol.
#' @param regulonsAucThresholds Output from get_regulonThresholds(loom) from SCENICoutput.
#' @keywords SCENIC, regulons, binary activity, kmeans, thresholds
#' @import SCENIC
#' @import dplyr
#' @import svMisc
#' @export
#' @examples
#' \donttest{
#' binary_regulons <- binarize_regulons(regulonAUC,regulon_thresholds)
#' }


## Binarize regulons based on thresholds
binarize_regulons_nextflow <- function(regulonsAUC,
                                       regulonsAucThresholds){

  binary_regulon_list <- list()

  for(regulon_no in 1:length(names(regulonsAucThresholds))){

    progress(regulon_no)

    regulon <- names(regulonsAucThresholds)[regulon_no]

    auc_df <-  data.frame("auc" = regulonsAUC[regulon,],
                          "cells"= colnames(regulonsAUC))

    auc_df <- auc_df %>%
      mutate("regulon"= if_else(auc >= regulonsAucThresholds[regulon],1,0)) %>%
      select(-auc)

    colnames(auc_df) <- c("cells",regulon)

    binary_regulon_list[[regulon]] <- auc_df
  }

  return(binary_regulon_list)
}
