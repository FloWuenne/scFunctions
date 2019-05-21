#' Calculates Regulon specificity score (RSS) from binary regulon activity.
#'
#' @param metadata Dataframe containing metadata about cells. Has to create a column named cell_type that assigns groupings to cells.
#' Can be the meta.data slot from a Seurat object.
#' @param binary_regulons Data frame with binary regulons, where regulons are rows and columns are cells. Can be created from output of binarize_regulons().
#' @keywords SCENIC, regulons, binary activity, kmeans, thresholds
#' @export
#' @examples
#'


## Iterate over all cell types and perform jensen shannon divergence test using binary regulon activity and genotype

calculate_rrs <- function(metadata,
                          binary_regulons){

  require(philentropy)
  require(SCENIC)
  require(svMisc)
  require(dplyr)
  require(tidyr)

  cell_types<- unique(metadata$cell_type)
  regulons <- rownames(binary_regulons)


  jsd_matrix_ct <- data.frame("regulon" = c(),
                              "cell_type" = c(),
                              "jsd" = c())


  cell_type_counter <- 0
  for(ct in unique(cell_types)) {

    cell_type_counter <- cell_type_counter + 1
    print(paste("Processing cell type:",cell_type_counter,ct,sep=" "))

    for(regulon_no in 1:length(regulons)) {

      progress(regulon_no)

      regulon <- regulons[regulon_no]

      regulon_vec <- binary_regulons[regulon,]
      regulon_norm <- regulon_vec/sum(regulon_vec)

      genotype_vec <- metadata[colnames(binary_regulons),]
      genotype_vec <- genotype_vec %>%
        mutate("cell_class" = if_else(cell_type == ct,1,0))
      genotype_vec <- genotype_vec$cell_class
      genotype_norm <- genotype_vec/sum(genotype_vec)

      dist_df <- rbind(regulon_norm,genotype_norm)

      ## Calculate the Jensen-Shannon divergence
      jsd_divergence <- suppressMessages(philentropy::JSD(dist_df))

      ## Calculate Jensen-Shannon distance
      rss <- 1-sqrt(jsd_divergence)

      regulon_jsd <- data.frame("regulon" = regulon,
                                "cell_type" = ct,
                                "RSS" = rss[1])

      jsd_matrix_ct <- rbind(jsd_matrix_ct,regulon_jsd)
    }

  }
  jsd_matrix_ct <- jsd_matrix_ct %>%
    arrange(desc(RSS))

  return(jsd_matrix_ct)
}
