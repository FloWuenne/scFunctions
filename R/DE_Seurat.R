#' Function to perform differential expression analysis for all clusters in a Seurat object.
#'
#' This function will take a precomputed Seurat object and perform differential expression analysis using one of the differential expression tests
#' included in Seurat (default= wilcox). If you want to perform DE analysis using edgeR, please check the function DE_edgeR_Seurat()!
#' All the results will be saved in a folder above the current folder location named DE_Seurat (../DE_Seurat). The output folder can easily be
#' modified using the parameter 'output_dir'.
#'
#' @param seurat_object The S4 Seurat object which contains filtered and normalized cells in the data slot.
#' @param de_function The function that will be used to perform differential expression analysis. See ?FindMarkers in the Seurat package for all options.
#' @param output_dir The relative directory that will be used to save results.
#' @param de_groups The two group labels to use for differential expression, supplied as a vector.
#' @keywords Seurat, DE, differential expression
#' @export
#' @examples
#' DE_Seurat()

## dependencies:
## Seurat : https://github.com/satijalab/seurat
## MAGIC : https://github.com/KrishnaswamyLab/MAGIC

DE_Seurat <- function(seurat_object,
                      de_function='wilcox',
                      output_dir= "../DE_Seurat",
                      grouping_var = "Genotype",
                      de_groups = c("WT","KO"),
                      min_pct = 0.1)
  {

  expression_seurat <- readRDS(seurat_object)

  ## print start message
  print("Starting differential expression analysis")

  ## Initiate empty data frames and lists for comparisons of clusters
  joined_res_table <- data.frame()
  upset_Rlist_DE_genes <- list()

  ## Set cluster numbers to keep track of how many clusters have been processed
  cluster_number <- 0

  ## Iterate over each cluster in the @ident slot
  for(this_cluster in sort(unique(seurat_objects.combined@ident))){

    cluster_number <- cluster_number + 1

    ## Print status for which cluster identity is being processed
    print(paste("Working on cluster #",cluster_number,":",this_cluster,sep=""))

    ## Subset Seurat object to only contain cells from this cluster
    cells_in_this_cluster <- SubsetData(seurat_objects.combined,
                                        ident.use=this_cluster)

    ## Get vector of names for WT and ko cells
    cells_group_1 <- rownames(subset(cells_in_this_cluster@meta.data,grouping_var == de_groups[1]))
    cells_group_2 <- rownames(subset(cells_in_this_cluster@meta.data,grouping_var == de_groups[2]))

    #### EdgeR DE expression analysis
    ## Check whether there are cells in both groups, otherwise skip this cluster
    if(length(cells_group_1) > 1 & length(cells_group_2) > 1){

      cells_in_this_cluster <- StashIdent(cells_in_this_cluster, save.name = "OldIdent")
      cells_in_this_cluster <- SetAllIdent(cells_in_this_cluster, id = grouping_var)

      ## Perform differential expression test using the Seurat FindMarkers function
      this_cluster_de_genes <- data.frame()
      this_cluster_de_genes <- FindMarkers(cells_in_this_cluster,
                                                  ident.1 = de_groups[1],
                                                  ident.2 = de_groups[2],
                                                  print.bar = TRUE,
                                                  test.use = de_function,
                                                  min.pct = min_pct)

      ## Write table for all differentially expressed genes containing testing results
      write.table(this_cluster_de_genes,
                  file=paste("../DE_Seurat/Cluster_",this_cluster,"_significant_DE_genes.",de_function,".txt",sep=""),
                  sep="\t",
                  quote=FALSE,
                  row.names=TRUE,
                  col.names=TRUE)


      ## Save DE results in a joined table
      this_cluster_de_genes_wilcox$cluster <- replicate(nrow(this_cluster_de_genes_wilcox),this_cluster)

      if(cluster_number == 1){
        joined_res_table <- rbind(joined_res_table,this_cluster_de_genes_wilcox)
      }else{
        joined_res_table <- rbind(joined_res_table,this_cluster_de_genes_wilcox[2:nrow(this_cluster_de_genes_wilcox)])
      }

      ## Calculate cell type average expressions to check correlation between the two groups
      avg.cells_in_this_cluster <- log1p(AverageExpression(cells_in_this_cluster, show.progress = FALSE))
      avg.cells_in_this_cluster$gene <- rownames(avg.cells_in_this_cluster)


      ## Make a correlation plot between the two conditions
      corr_plot <- ggplot(avg.cells_in_this_cluster, aes(de_groups[1], de_groups[2],text=gene)) +
        geom_point() +
        ggtitle(paste("Cluster_",this_cluster,sep="")) +
        geom_abline(intercept = 0, slope = 1, col="red")

      ## Save normal png version of the plot
      ggsave(corr_plot,
             file=paste("../DE_Seurat/Cluster_",this_cluster,"_corrplot.png",sep=""))

      ## Also make an interactive version using plotly
      library(plotly)

      htmlwidgets::saveWidget(as.widget(ggplotly(corr_plot)), paste("../DE_Seurat/Cluster_",this_cluster,"_corrplot.plotly.html",sep=""))

      # ## Plot volcano plot
      # volcano_plot <- ggplot(this_cluster_de_genes,aes(logFC,-log(PValue))) +
      #   geom_point(size=2) +
      #   geom_point(data = subset(res_table,(logFC < neg_log_FC_thresh) & (PValue < q_value_thresh)),col="red") +
      #   geom_point(data = subset(res_table,(logFC > pos_log_FC_thresh ) & (PValue < q_value_thresh) ),col="green") +
      #   geom_text_repel(
      #     data = subset(res_table, (logFC > pos_log_FC_thresh | logFC < neg_log_FC_thresh) & (PValue < q_value_thresh)),
      #     aes(label = subset(res_table, (logFC > pos_log_FC_thresh | logFC < neg_log_FC_thresh) & (PValue < q_value_thresh))$gene),
      #     size = 5,
      #     box.padding = unit(0.35, "lines"),
      #     point.padding = unit(0.3, "lines")
      #   ) +
      #   ggtitle(res$comparison)+
      #   theme_light()
      #
      # volcano_plot
      # ggsave(volcano_plot,file=paste("../DE_edgeR/cluster-",this_cluster,"_volcano_plot.svg",sep=""))

      print("Finished with this cluster!")

    }
    else {
      print(paste("Cluster",this_cluster," only contains cells from one group!",sep=""))
    }

  }

  ## Write table for all differentially expressed genes containing testing results
  write.table(joined_res_table,
              file=paste("../DE_Seurat/All_DE_genes.tsv",sep=""),
              sep="\t",
              quote=FALSE,
              row.names=TRUE,
              col.names=TRUE)

  }
