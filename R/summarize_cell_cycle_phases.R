#' Function to summarize and analyze precomputed cell cycle phases in a Seurat object between replicates or cell types.
#'
#' This function takes as input a precomputed Seurat object that contains cell cycle phases for each cell. It then caclulates distribution
#' of cells phases by one or more user defined grouping variables and performs tests to see whether there is a statistical difference in
#' cell cycle distribution.
#'
#' @param seurat_object The S4 Seurat object which contains filtered and normalized cells in the data slot.
#' @param cell_cycle_annotation The column in the meta.data table that contains cell cycle phase information.
#' @param Replicate_column The column in the meta.data table that contains biological replicate information.
#' @param group_by The column that contains the groups to be compared.
#' @param save_plot logical of whether to save the plot as png. Default = FALSE.
#' @param output_dir Output directory of where to save the plot and corresponding table.
#' @param grouping_var_order User supplied vector of order for grouping variable.
#' @keywords Seurat, cell cycle, compare replicates
#' @export
#' @examples
#' summarize_cell_cycle_phases()
## dependencies:
## Seurat : https://github.com/satijalab/seurat
## ggplot2 : https://ggplot2.tidyverse.org/
## RColorBrewer : https://cran.r-project.org/web/packages/RColorBrewer/RColorBrewer.pdf
## dplyr : https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html


summarize_cell_cycle_phases <- function(seurat_object,
                                        cell_cycle_annotation,
                                        Replicate_column,
                                        grouping_var,
                                        grouping_var_order = NULL,
                                        save_plot = FALSE,
                                        output_dir = ".",
                                        per_cell_type=FALSE){

  ## Load libraries
  require(ggplot2)
  require(Seurat)
  require(dplyr)
  require(RColorBrewer)

  ## print start message
  print("Starting to summarize cell cycle phases...")

  ## Add cluster identities to metadata
  seurat_object@meta.data$Clusters <- seurat_object@ident

  ## Check if user wants to compute statistics per cell type or over the entire dataset
  if(per_cell_type == FALSE){

    ## Count the number of cells per cell cycle for each replicate per grouping variable without keeping cell cycle
    cell_cycle_stats_per_group <-  seurat_object@meta.data %>%
      group_by(get(grouping_var),get(Replicate_column),get(cell_cycle_annotation)) %>%
      summarise (n = n()) %>%
      mutate(freq = n / sum(n))

    # Set column names of summary table
    colnames(cell_cycle_stats_per_group) <- c("Grouping_var","Replicate","Phase","n","freq")

    ## Set order of Genotypes in plot
    if(is.null(grouping_var_order)){
      cell_cycle_staats_per_group$Grouping_var <- factor(cell_cycle_stats_per_group$Grouping_var ,levels=sort(unique(cell_cycle_stats_per_group$Grouping_var)))
    } else {
      cell_cycle_stats_per_group$Grouping_var <- factor(cell_cycle_stats_per_group$Grouping_var ,levels=grouping_var_order)
    }


    ## Plot boxplot for cell cycle phases per group
    cell_cycle_plot <- ggplot(cell_cycle_stats_per_group,aes(Phase,freq,fill=Grouping_var)) +
      geom_boxplot() +
      scale_fill_brewer(name= grouping_var,
                        palette = "Set1") +
      labs(x = "Cell Cycle Phase",
           y = "% of cells")


  }else if (per_cell_type == TRUE){

    ## Count the number of cells per cell cycle for each replicate per grouping variable without keeping cell cycle
    cell_cycle_stats_per_group <-  seurat_object@meta.data %>%
      group_by(get(grouping_var),get(Replicate_column),Clusters,get(cell_cycle_annotation)) %>%
      summarise (n = n()) %>%
      mutate(freq = n / sum(n))

    # Set column names of summary table
    colnames(cell_cycle_stats_per_group) <- c("Grouping_var","Replicate","Clusters","Phase","n","freq")

    ## Set order of Genotypes in plot
    if(is.null(grouping_var_order)){
      cell_cycle_stats_per_group$Grouping_var <- factor(cell_cycle_stats_per_group$Grouping_var ,levels=sort(unique(cell_cycle_stats_per_group$Grouping_var)))
    } else {
      cell_cycle_stats_per_group$Grouping_var <- factor(cell_cycle_stats_per_group$Grouping_var ,levels=grouping_var_order)
    }


    ## Plot boxplot for cell cycle phases per group
    cell_cycle_plot <- ggplot(cell_cycle_stats_per_group,aes(Phase,freq,fill=Grouping_var)) +
      geom_boxplot() +
      facet_wrap(~ Clusters) +
      scale_fill_brewer(name= grouping_var,
                        palette = "Set1") +
      labs(x = "Cell Cycle Phase",
           y = "% of cells")

  }

  ## Check if user wants to save plot
  if(save_plot == TRUE){

    ggsave(cell_cycle_plot,
           file=paste(output_dir,"/Cell_cycle_comparison.png",sep=""),
           device = "png",
           width=10,
           height=6,
           units= "in")

    write.table(cell_cycle_stats_per_group,
                file = paste(output_dir,"/Cell_cycle_comparison.txt",sep=""),
                col.names=TRUE,
                row.names=FALSE,
                quote=FALSE,
                sep="\t")

  }

  return(cell_cycle_plot)

}
