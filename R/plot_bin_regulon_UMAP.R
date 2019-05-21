#' Calculates Regulon specificity score (RSS) from binary regulon activity.
#'
#' @param seurat_object Processed seurat object (v3) with clustering results.
#' @param binary_regulons Binary regulon activity as named list.
#' @param regulon_name Name of the regulon to plot.
#' @param orientation Orientation of the plots. Wide = UMAP plots side-by-side, long = UMAP plots stacked.
#' @param highlight_color Color with which to highlight the active cells for a regulon
#' @keywords SCENIC, regulons, RRS, cell type classification
#' @export
#' @examples
#'

## Binary regulon activity on top of UMAP plot
plot_bin_regulon_UMAP <- function(seurat_object,
                             binary_regulons,
                             regulon_name,
                             orientation = "wide",
                             highlight_color = "red"){

  require(Seurat)
  require(tidyverse)

  ## Plot binary regulon on top of UMAP
  umap_mapping <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
  umap_mapping$cells <- rownames(umap_mapping)
  cell_types <- seurat_object@meta.data %>%
    mutate("cells" = rownames(seurat_object@meta.data)) %>%
    select(cells,cell_type)

  regulon <- binary_regulons[["Lef1_extended (286g)"]]
  colnames(regulon) <- c("cells","regulon_activity")
  regulon_umap <- full_join(umap_mapping,regulon,by="cells")
  regulon_umap <- full_join(regulon_umap,cell_types,by="cells")

  cell_type_umap <- DimPlot(seurat_object,label = TRUE) +
    theme(legend.position = "none")

  binary_umap <- ggplot(regulon_umap,aes(UMAP_1,UMAP_2,colour=as.factor(regulon_activity))) +
    geom_point(data = subset(regulon_umap,regulon_activity == 0),aes(colour = as.factor(regulon_activity)), alpha =0.75) +
    geom_point(data = subset(regulon_umap,regulon_activity == 1), aes(colour = as.factor(regulon_activity)), alpha =0.75) +
    scale_color_manual("Activity",
                       labels = c("OFF","ON"),
                       values = c("grey50", highlight_color))

  legend <- get_legend(binary_umap)


  if(orientation == "wide"){
    plot_grid(cell_type_umap,binary_umap + theme(legend.position="none"),legend,
              nrow =1,
              labels = c("A","B"),
              rel_widths = c(3,3,1))

  }else if(orientation == "long"){
    plot_grid(cell_type_umap,binary_umap,legend,
              ncol = 1,
              labels = c("A","B"),
              rel_widths = c(3,3,1))
  }

}
