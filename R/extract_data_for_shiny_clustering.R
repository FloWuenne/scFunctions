#' Function to perform differential expression analysis for all clusters in a Seurat object.
#'
#' This function will take a precomputed Seurat object and perform differential expression analysis using one of the differential expression tests
#' included in Seurat (default= wilcox). If you want to perform DE analysis using edgeR, please check the function DE_edgeR_Seurat()!
#' All the results will be saved in a folder above the current folder location named DE_Seurat (../DE_Seurat). The output folder can easily be
#' modified using the parameter 'output_dir'.
#'
#' @param seurat_object The S4 Seurat object which contains filtered and normalized cells in the data slot.
#' @param time_points The time points that should be processed
#' @param imputed logical that indicates whether data has been imputed and there is a data frame in @imputed or not. default = TRUE
#' @param output_dir Directory where the subclustering module and all gene Rds objects will be saved. default = "."
#' @keywords Seurat, Rshiny, heart maturation, speed optimization
#' @import Seurat
#' @import shiny
#' @import DT
#' @export
#' @examples
#' \donttest{
#' extract_data_for_shiny_clustering()
#' }

## dependencies:
## Seurat : https://github.com/satijalab/seurat

extract_data_for_shiny_clustering <- function(seurat_object,
                                              time_points = c("E14.5","E16.5","E18.5","P1","P4","P7"),
                                              imputed = TRUE,
                                              output_dir = ".") {

  for(sample in time_points){
    ## Print a message
    print(paste("Now working on sample: ",sample,sep=""))

    ## Reformat seurat object to only keep information required for mapping cells

    ## Normalized expression data
    norm_exprs_sparse <- seurat_object@data
    norm_exprs_matrix <- as.matrix(seurat_object@data)
    norm_exprs_df <- as.data.frame(norm_exprs_matrix)

    ## Create a new data frame that contains tSNE embeddings and cluster identities
    tsne_mappings <- seurat_object@dr$tsne@cell.embeddings
    cell_identities <- data.frame(seurat_object@ident)
    rownames(cell_identities) <- names(seurat_object@ident)
    tsne_mappings <- merge(tsne_mappings,cell_identities,by=0,all=TRUE)
    rownames(tsne_mappings) <- tsne_mappings$Row.names
    tsne_mappings <- tsne_mappings %>%
      select(-Row.names)

    ## Metadata information
    metadata <- seurat_object@meta.data

  }


}
