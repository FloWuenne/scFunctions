#' Calculates CSI values for all regulon pairs
#'
#' @param regulonAUC The AUC values for all regulons as calculated by SCENIC (content of file:3.4_regulonAUC.Rds).
#' @keywords SCENIC, regulons, binary activity, kmeans, thresholds
#' @export
#' @examples
#' regulon_thresholds <- auc_thresh_kmeans(regulonAUC)

calculate_csi <- function(regulonAUC,
                          calc_extended = FALSE,
                          verbose = FALSE){

  require(tidyverse)
  require(svMisc)

  compare_pcc <- function(vector_of_pcc,pcc){
    pcc_larger <- length(vector_of_pcc[vector_of_pcc > pcc])
    if(pcc_larger == length(vector_of_pcc)){
      return(0)
    }else{
      return(length(vector_of_pcc))
    }
  }

  calc_csi <- function(reg,reg2,pearson_cor){
    test_cor <- pearson_cor[reg,reg2]
    total_n <- ncol(pearson_cor)
    pearson_cor_sub <- subset(pearson_cor,rownames(pearson_cor) == reg | rownames(pearson_cor) == reg2)

    sums <- apply(pearson_cor_sub,MARGIN = 2, FUN = compare_pcc, pcc = test_cor)
    fraction_lower <- length(sums[sums == nrow(pearson_cor_sub)]) / total_n
    return(fraction_lower)
  }

  regulonAUC_sub <- regulonAUC@assays$data@listData$AUC

  if(calc_extended == TRUE){
    regulonAUC_sub <- subset(regulonAUC_sub,grepl("extended",rownames(regulonAUC_sub)))
  } else if (calc_extended == FALSE){
    regulonAUC_sub <- subset(regulonAUC_sub,!grepl("extended",rownames(regulonAUC_sub)))
}

  regulonAUC_sub <- t(regulonAUC_sub)

  pearson_cor <- cor(regulonAUC_sub)
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>%
    gather(regulon_2,pcc,-regulon_1) %>%
    mutate("regulon_pair" = paste(regulon_1,regulon_2,sep="_"))


  regulon_names <- unique(colnames(pearson_cor))
  num_of_calculations <- length(regulon_names)*length(regulon_names)

  csi_regulons <- data.frame(matrix(nrow=num_of_calculations,ncol = 3))

  colnames(csi_regulons) <- c("regulon_1",
                              "regulon_2",
                              "CSI")

  num_regulons <- length(regulon_names)

  f <- 0
  for(reg in regulon_names){
    ## Check if user wants to print info
    if(verbose == TRUE){
      print(reg)
      }
    for(reg2 in regulon_names){
      f <- f + 1

      fraction_lower <- calc_csi(reg,reg2,pearson_cor)

      csi_regulons[f,] <- c(reg,reg2,fraction_lower)

    }
  }
  csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
  return(csi_regulons)
}
