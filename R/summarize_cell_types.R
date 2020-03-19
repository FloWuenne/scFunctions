#' Function to summarize and analyze the number of cells and fraction of each cell type over multiple replicates in a Seurat object.
#'
#'
#' @param seurat_object The S4 Seurat object which contains filtered and normalized cells in the data slot.
#' @param Replicate_column The column in the meta.data table that contains biological replicate information.
#' @param group_by The column that contains the groups to be compared.
#' @param grouping_var_order User supplied vector of order for grouping variable.
#' @param save_plot logical of whether to save the plot as png. Default = FALSE.
#' @param output_dir Output directory of where to save the plot and corresponding table.
#' @keywords Seurat, cell type summary, compare replicates
#' @export
#' @examples

## dependencies:
## Seurat : https://github.com/satijalab/seurat
## ggplot2 : https://ggplot2.tidyverse.org/
## RColorBrewer : https://cran.r-project.org/web/packages/RColorBrewer/RColorBrewer.pdf
## dplyr : https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html


summarize_cell_types <- function(seurat_object,
                                        Replicate_column,
                                        grouping_var,
                                        grouping_var_order = NULL,
                                        save_plot = FALSE,
                                        output_dir = "."){

  ## Load libraries
  require(ggplot2)
  require(Seurat)
  require(dplyr)
  require(RColorBrewer)

  ## print start message
  print("Starting to summarize cell cycle phases...")

  ## Add cluster identities to metadata
  seurat_object@meta.data$Clusters <- seurat_object@ident

  ## Count the number of cells per cell cycle for each replicate per grouping variable without keeping cell cycle
  cell_type_stats_per_group <-  seurat_object@meta.data %>%
    group_by(get(grouping_var),get(Replicate_column),Clusters) %>%
    summarise (n = n()) %>%
    mutate(freq = n / sum(n))

  # Set column names of summary table
  colnames(cell_type_stats_per_group) <- c("Grouping_var","Replicate","Clusters","n","freq")

  ## Set order of Genotypes in plot
  if(is.null(grouping_var_order)){
    cell_type_stats_per_group$Grouping_var <- factor(cell_type_stats_per_group$Grouping_var ,levels=sort(unique(cell_type_stats_per_group$Grouping_var)))
  } else {
    cell_type_stats_per_group$Grouping_var <- factor(cell_type_stats_per_group$Grouping_var ,levels=grouping_var_order)
  }

    ## Plot boxplot for cell cycle phases per group
    cell_type_summary_plot <- ggplot(cell_type_stats_per_group,aes(Clusters,freq,fill=Grouping_var)) +
      geom_boxplot() +
      scale_fill_brewer(name= grouping_var,
                        palette = "Set1") +
      labs(x = "Cell cluster",
           y = "% of cells")


  ## Check if user wants to save plot
  if(save_plot == TRUE){

    ggsave(cell_type_summary_plot,
           file=paste(output_dir,"/Cell_type_comparison_summary.png",sep=""),
           device = "png",
           width=10,
           height=6,
           units= "in")

    write.table(cell_type_stats_per_group,
                file = paste(output_dir,"/Cell_type_comparison_summary.txt",sep=""),
                col.names=TRUE,
                row.names=FALSE,
                quote=FALSE,
                sep="\t")

  }

  return(cell_type_summary_plot)

}
