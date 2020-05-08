#' A function to perform imputation on a Seurat object using MAGIC.
#'
#' This function enables you to easily impute the normalized gene expression of a Seurat object using the MAGIC package in R.
#' More information about MAGIC can be found in the github repository.
#' @param seurat_object The S4 Seurat object which contains filtered and normalized cells in the data slot.
#' @keywords Seurat, MAGIC, imputation
#' @import Seurat
#' @import Rmagic
#' @export
#' @examples
#' \donttest{
#' impute_seurat_MAGIC()
#' }

## This function takes as input a Seurat object and uses MAGIC to perform imputation of the data.
## Imputation is performed on the filtered and normalized cells in the @data slot of the seurat object!
## The imputed expression matrix will then be saved in the @data

## dependencies:
## Seurat : https://github.com/satijalab/seurat
## MAGIC : https://github.com/KrishnaswamyLab/MAGIC

impute_seurat_MAGIC <- function(seurat_object){

  ## Print status message
  print("Starting magic imputation on your seurat object...")

  ## transpose matrix for MAGIC
  transposed_data <- t(as.matrix(seurat_object@data))

  ## Actual imputation function
  MAGIC_data <- magic(transposed_data,  genes='all_genes')

  ## Retranspose the result for proper Seurat matrix format
  imputed_data <- t(MAGIC_data$result)

  ## Replace old Seurat data with imputed matrix
  seurat_object@imputed <- as.data.frame(imputed_data)

  ## Print status message
  print("Mieschief managed, MAGIC has been run on your data! ")

  return(seurat_object)
}
