#' Function to export data tables from Seurat object.
#'
#' This simple function will save the raw UMI matrix (seurat_object@raw.data), the normalized UMI matrix (seurat_object@data)
#' and the metadata (seurat_object@meta.data) from the Seurat object into tab separated files.
#'
#' @param seurat_object The S4 Seurat object which contains filtered and normalized cells in the data slot.
#' @param output_dir Directory where the data tables will be saved.
#' @keywords Seurat, DE, differential expression
#' @export
#' @examples
#' export_data_from_seurat()

## Simple function that will export the raw expression matrix, the normalized expression as well as the metadata table from a Seurat object
## All files will be saved as tsv files. The output directory can be specified using the output_dir parameter.
export_data_from_seurat <- function(seurat_object,
                                    output_dir = ".")
{
  ## Save raw UMI matrix with all cells
  write.table(as.data.frame(as.matrix(seurat_object@raw.data)),
              file=paste(output_dir,"/","Raw_UMI_Matrix.tsv",sep=""),
              sep="\t",
              row.names=TRUE,
              col.names=TRUE,
              quote=FALSE)

  ## Save normalized UMI matrix with all cells
  write.table(as.data.frame(as.matrix(seurat_object@data)),
              file=paste(output_dir,"/","NormalizedUMI_Matrix.tsv",sep=""),
              sep="\t",
              row.names=TRUE,
              col.names=TRUE,
              quote=FALSE)

  ## Save metadata table
  write.table(seurat_object@meta.data,
              file=paste(output_dir,"/","Metadata.tsv",sep=""),
              sep="\t",
              row.names=TRUE,
              col.names=TRUE,
              quote=FALSE)

}
