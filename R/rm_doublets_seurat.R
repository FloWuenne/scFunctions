#' Uses the doubletFinder package to identify and remove potential cell doublets from a precomputed Seurat object
#'
#' This function takes as input a precomputed Seurat object, that contains initial clustering results on the data.
#' It then uses the precomputed clusters with the doubletFinder package, which calculates in silico doublets to identify
#' cells that likely represent doublets. The function will create a QC graph and a metadatable with assignment of cells
#' as either singlet or doublet. the function returns a cleaned seurat object in which doublets have been removed.
#'
#' Currently uses v2 Seurat objects as input!
#'
#' @param seurat_object The S4 Seurat object which contains filtered and normalized cells in the data slot.
#' @param percent_doublets Estimated percentage of cell doublets. This number will be used with the number of cells in the seurat objects to estimate doublets. default = 5% (0.05)
#' @param proportion_artificial Parameter for duobletFinder that c
#' @param outfile_prefix The prefix that should be used for the output files created.
#' @param replicate_column The column in the meta.data table that contains biological replicate information.
#' @param cell_cycle_column The column that contains the cell cycle annotation for each cell. If cell cycle has not been calculated, set to FALSE.
#' @param save_plot logical of whether to savea QC plot. Default = FALSE.
#' @param output_dir Output directory of where to save the plot and corresponding table. default = ".".
#' @keywords Seurat, cell doublets, QC, doubletFinder
#' @export
#' @examples
#' rm_doublets_seurat()

## dependencies:
## Seurat : https://github.com/satijalab/seurat
## ggplot2 : https://ggplot2.tidyverse.org/
## RColorBrewer : https://cran.r-project.org/web/packages/RColorBrewer/RColorBrewer.pdf
## dplyr : https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html
## doubletFinder: https://github.com/chris-mcginnis-ucsf/DoubletFinder


summarize_cell_cycle_phases <- function(seurat_object,
                                        outfile_prefix = "doublets_",
                                        replicate_column = "Replicate",
                                        cell_cycle_column = "Phase",
                                        save_plot = TRUE,
                                        output_dir = "."){

  ## Load libraries
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(doubletFinder)

  #### Identify doublets using DoubletFinder

  print("Identifying cell doublets using DoubletFinder...")


  ## pK Identification -----------------------------------------------------------------------------
  sweep.res.list_seurat_object <- doubletFinder_ParamSweep(seurat_object)
  sweep.stats_seurat_object <- summarizeSweep(sweep.res.list_seurat_object, GT = FALSE)
  bcmvn_seurat_object <- find.pK(sweep.stats_seurat_object)

  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- small_seurat@ident
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*length(seurat_object@cell.names))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seurat_object <- doubletFinder(seurat_object, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE)


  return(seurat_object)

}
