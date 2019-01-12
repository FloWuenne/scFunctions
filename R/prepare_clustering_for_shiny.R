#' Function to prepare as RDS object for http://andelfingerlab.heart_maturation.genap.ca from a precomputed Seurat object.
#' The new object will have optimized size and data structure for usage in the Rshiny server.
#'
#' @param seurat_object The S4 Seurat object which contains filtered and normalized cells in the data slot.
#' @param label The sample label that will be used for naming the file.
#' @param imputed logical that indicates whether data has been imputed and there is a data frame in @imputed or not. default = FALSE
#' @param output_dir Directory in which to save the new RDS object. default = "."
#' @keywords Seurat, Rshiny, heart maturation
#' @export
#' @examples
#' extract_data_for_shiny_clustering()

## dependencies:
## Seurat : https://github.com/satijalab/seurat

prepare_clustering_for_shiny <- function(seurat_object,
                                         marker_table,
                                         label = "E14.5",
                                         imputed = FALSE,
                                         output_dir = ".",
                                         marker_table) {

  ## Load libraries
  require(Seurat)
  require(shiny)
  require(DT)
  require(dplyr)

  ## Define S4 object with the required components for the shiny server
  clustering_s4_object <- setClass("clustering_module",
                                   slots=c(name = "character",
                                           mappings="data.frame",
                                           centroids = "data.frame",
                                           metadata="data.frame",
                                           norm_exprs="data.frame",
                                           marker_list="data.frame"))

  ## Reformat seurat object to only keep information required for mapping cells

  ## Check which Seurat version was used to make the object
  seurat_version <- seurat_object@version

  ## Code for Seurat v2
  if(seurat_version == "2.3.4"){

    print(paste("Seurat version:",seurat_version,"detected!",sep=" "))

    ## Normalized expression data
    norm_exprs_sparse <- seurat_object@data
    norm_exprs_matrix <- as.matrix(seurat_object@data)
    norm_exprs_df <- as.data.frame(norm_exprs_matrix)

    ## Create a new data frame that contains  dimensional embeddings (tSNE, UMAP) and cluster identities
    tsne_mappings <- seurat_object@dr$tsne@cell.embeddings
    umap_mappings <- seurat_object@dr$umap@cell.embeddings

    cell_identities <- data.frame(seurat_object@ident)
    rownames(cell_identities) <- names(seurat_object@ident)

    tsne_mappings <- merge(tsne_mappings,cell_identities,by=0,all=TRUE)
    rownames(tsne_mappings) <- tsne_mappings$Row.names
    tsne_mappings <- tsne_mappings %>%
      select(-Row.names)

    final_mappings <- merge(tsne_mappings,umap_mappings,by=0,all=TRUE)

    rownames(final_mappings) <- final_mappings$Row.names
    final_mappings <- final_mappings %>%
      select(-Row.names)

    ## Calculate centers of each cluster for plotting
    centroids_tsne <- final_mappings %>%
      dplyr::group_by(seurat_object.ident) %>%
      summarise_at(vars(tSNE_1,tSNE_2),funs(mean(., na.rm=TRUE)))
    centroids_tsne <- as.data.frame(centroids_tsne)

    centroids_umap <- final_mappings %>%
      dplyr::group_by(seurat_object.ident) %>%
      summarise_at(vars(UMAP1,UMAP2),funs(mean(., na.rm=TRUE)))
    centroids_umap <- as.data.frame(centroids_umap)

    centroids <- merge(centroids_tsne,centroids_umap,by="seurat_object.ident")

    ## Metadata information
    metadata <- seurat_object@meta.data
  }

  ## Code for Seurat v3
  ## To come!

  ####

  ## Fill the slots of the S4 object with data
  current_object <- clustering_s4_object(name = label,
                                          mappings=tsne_mappings,
                                          centroids = centroids,
                                          norm_exprs=norm_exprs_df,
                                          metadata = metadata,
                                          marker_list = marker_table)

  ## Save new S4 object containing all the necessary data for the clustering module
  saveRDS(current_object,file=paste(output_dir,"/",label,".clustering_module.Rds",sep=""))

}
