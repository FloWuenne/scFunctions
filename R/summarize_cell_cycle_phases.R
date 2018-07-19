#'
#'
#'
#'
#' @param seurat_object The S4 Seurat object which contains filtered and normalized cells in the data slot.
#' @param de_function The function that will be used to perform differential expression analysis. See ?FindMarkers in the Seurat package for all options.
#' @param output_dir The relative directory that will be used to save results.
#' @param de_groups The two group labels to use for differential expression, supplied as a vector.
#' @param clusters_to_exclude Define a vector of clusters for which you don't want to perform DE analysis.
#' @param number_of_intersections_to_show Define the number of intersections to plot in Upset R plot. default = show all intersections.
#' @keywords Seurat, DE, differential expression
#' @export
#' @examples
#' summarize_cell_cycle_phases()

## dependencies:
## Seurat : https://github.com/satijalab/seurat
## Plotly : https://plot.ly/r/
## ggplot2 : https://ggplot2.tidyverse.org/
