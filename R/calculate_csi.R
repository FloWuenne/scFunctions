#' Calculates CSI values for all regulon pairs
#'
#' @param regulonAUC The AUC values for all regulons as calculated by SCENIC (content of file:3.4_regulonAUC.Rds).
#' @keywords SCENIC, regulons, binary activity, kmeans, thresholds
#' @export
#' @examples
#' regulon_thresholds <- auc_thresh_kmeans(regulonAUC)

calculate_csi <- function(regulonAUC){

  require(tidyverse)

  regulonAUC_sub <- t(regulonAUC@assays$data@listData$AUC)

  pearson_cor <- cor(regulonAUC_sub)
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>%
    gather(regulon_2,pcc,-regulon_1) %>%
    mutate("regulon_pair" = paste(regulon_1,regulon_2,sep="_"))


  csi_regulons <- data.frame("regulon_1" = c(),
                             "regulon_2" = c(),
                             "CSI" = c())

  num_regulons <- length(unique(pearson_cor_long$regulon_1))
  for(reg in unique(pearson_cor_long$regulon_1)){
    print(reg)
    for(reg2 in unique(pearson_cor_long$regulon_2)){

      ## Don't calculate CSI for regulon with itself
      this_pair_pcc <- subset(pearson_cor_long,regulon_1 == reg & regulon_2 == reg2)$pcc

      pcc_lower1 <- pearson_cor_long %>%
        subset(regulon_1 == reg & regulon_2 != reg & regulon_2 != reg2) %>%
        summarise("pcc_lower" = matrixStats::count(pcc < this_pair_pcc))

      pcc_lower2 <- pearson_cor_long %>%
        subset(regulon_2 == reg2 & regulon_1 != reg2 & regulon_2 != reg) %>%
        summarise("pcc_lower" = matrixStats::count(pcc < this_pair_pcc))

      csi <- (pcc_lower1$pcc_lower + pcc_lower2$pcc_lower) / ((num_regulons - 2 )*2)
      this_csi <- data.frame("regulon_1" = reg,
                             "regulon_2" = reg2,
                             "CSI" = csi)

      csi_regulons <- rbind(csi_regulons,this_csi)

    }
  }
  return(csi_regulons)
}
