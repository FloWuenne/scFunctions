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
#' @param clusters_to_exclude Define a vector of clusters for which you don't want to perform DE analysis.
#' @keywords Seurat, DE, differential expression
#' @export
#' @examples
#' DE_Seurat()

## dependencies:
## Seurat : https://github.com/satijalab/seurat
## Plotly : https://plot.ly/r/
## ggplot2 : https://ggplot2.tidyverse.org/

DE_Seurat <- function(seurat_object,
                      de_function='wilcox',
                      output_dir= "../DE_Seurat",
                      grouping_var = "Genotype",
                      de_groups = c("WT","KO"),
                      min_pct = 0.1,
                      clusters_to_exclude = c())
  {

  ## Load libraries
  require(plotly)
  require(ggplot2)
  require(Seurat)
  require(UpSetR)
  require(dplyr)

  ## print start message
  print("Starting differential expression analysis")

  ## Initiate empty data frames and lists for comparisons of clusters
  joined_res_table <- data.frame()
  upset_Rlist_DE_genes <- list()

  ## Set cluster numbers to keep track of how many clusters have been processed
  cluster_number <- 0

  ## Get vector of clusters to test
  clusters_to_test <- sort(unique(seurat_object@ident))
  clusters_to_test <- setdiff(clusters_to_test,clusters_to_exclude)

  ## Iterate over each cluster in the @ident slot
  for(this_cluster in clusters_to_test){

    cluster_number <- cluster_number + 1

    ## Print status for which cluster identity is being processed
    print(paste("Working on cluster #",cluster_number,":",this_cluster,sep=""))

    ## Subset Seurat object to only contain cells from this cluster
    cells_in_this_cluster <- SubsetData(seurat_object,
                                        ident.use=this_cluster)

    ## Get vector of names for both groups of cells
    cells_group_1 <- rownames(subset(cells_in_this_cluster@meta.data,get(grouping_var) == de_groups[1]))
    cells_group_2 <- rownames(subset(cells_in_this_cluster@meta.data,get(grouping_var) == de_groups[2]))

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

      ## Filter out genes with an adjusted p value that are not significant
      this_cluster_de_genes <- subset(this_cluster_de_genes,p_val_adj < 0.05)

      ## Write table for all differentially expressed genes containing testing results
      write.table(this_cluster_de_genes,
                  file=paste(output_dir,"/",this_cluster,"_significant_DE_genes.",de_function,".txt",sep=""),
                  sep="\t",
                  quote=FALSE,
                  row.names=TRUE,
                  col.names=TRUE)


      ## Save DE results in a joined table
      this_cluster_de_genes$cluster <- replicate(nrow(this_cluster_de_genes),this_cluster)
      this_cluster_de_genes$gene <- rownames(this_cluster_de_genes)
      rownames(this_cluster_de_genes) <- c()
      joined_res_table <- rbind(joined_res_table,this_cluster_de_genes)

      ## Add DE for this cell type to the UpsetR list
      upset_Rlist_DE_genes[[this_cluster]] <- c(this_cluster_de_genes$gene)

      ## Calculate cell type average expressions to check correlation between the two groups
      avg.cells_in_this_cluster <- log1p(AverageExpression(cells_in_this_cluster, show.progress = FALSE))
      avg.cells_in_this_cluster$gene <- rownames(avg.cells_in_this_cluster)


      ## Make a correlation plot between the two conditions
      corr_plot <- ggplot(avg.cells_in_this_cluster, aes(get(de_groups[1]), get(de_groups[2]),text=gene)) +
        geom_point() +
        ggtitle(paste("Cluster : ",this_cluster,sep="")) +
        geom_abline(intercept = 0, slope = 1, col="red") +
        labs(x = de_groups[1],
             y= de_groups[2])

      ## Save normal png version of the plot
      ggsave(corr_plot,
             file=paste(output_dir,"/",this_cluster,"_corrplot.png",sep=""))

      ## Also make an interactive version using plotly
      library(plotly)

      htmlwidgets::saveWidget(as.widget(ggplotly(corr_plot)), paste(output_dir,"/",this_cluster,"_corrplot.plotly.html",sep=""))

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

  ## Check the overlap of DE genes between clusters using
  ## UpsetR: https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html


  ## The two functions fromList and get_intersect_members originate from this github post:
  ## https://github.com/hms-dbmi/UpSetR/issues/85
  ## Function to run UpsetR with a list of named vectors
  fromList <- function (input) {
    # Same as original fromList()...
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
      x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = F))
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(input)
    # ... Except now it conserves your original value names!
    row.names(data) <- elements
    return(data)
  }


  ## Plot the Upet plot
    svg(paste(output_dir,"/Overlap_DE_genes.svg",sep=""),
        width = 24,
        height=20)
    upset(fromList(upset_Rlist_DE_genes),
          order.by = "freq",
          sets = names(upset_Rlist_DE_genes),
          main.bar.color = "black",
          matrix.color="#1482a5ff",
          mainbar.y.label = "Number of DE genes",
          point.size = 6,
          line.size = 2,
          group.by = "degree",
          cutoff = 2,
          text.scale= 3)
    dev.off()

    ## Function to get the members of an intersection
    get_intersect_members <- function (x, ...){
      require(dplyr)
      require(tibble)
      x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
      n <- names(x)
      x %>% rownames_to_column() -> x
      l <- c(...)
      a <- intersect(names(x), l)
      ar <- vector('list',length(n)+1)
      ar[[1]] <- x
      i=2
      for (item in n) {
        if (item %in% a){
          if (class(x[[item]])=='integer'){
            ar[[i]] <- paste(item, '>= 1')
            i <- i + 1
          }
        } else {
          if (class(x[[item]])=='integer'){
            ar[[i]] <- paste(item, '== 0')
            i <- i + 1
          }
        }
      }
      do.call(filter_, ar) %>% column_to_rownames() -> x
      return(x)
    }

  ## Write table for all differentially expressed genes containing testing results
  write.table(joined_res_table,
              file=paste(output_dir,"/All_DE_genes.tsv",sep=""),
              sep="\t",
              quote=FALSE,
              row.names=TRUE,
              col.names=TRUE)

  ## Count number of clusters each gene is found DE in and write a table for easy inspection
  ## Count number of occurences for each gene and which cell types its in
  genes_accumulated <- joined_res_table %>%
    count(gene)

  colnames(genes_accumulated) <- c("gene","Clusters_DE")

  ## Count number of occurences for each gene and which cell types its in
  genes_clusters <- joined_res_table %>%
    group_by(gene) %>%
    summarise(cluster_names = paste(cluster, collapse=","))

  ## Calculate mean avgLogFC as well as the SD
  genes_clusters_means <- joined_res_table %>%
    group_by(gene) %>%
    summarise(mean_avgLogFC = mean(avg_logFC),
              SD_avgLogFC = sd(avg_logFC),
              mean_pvalue_adj = mean(p_val_adj),
              SD_pvalue_adj = sd(p_val_adj),
              mean_pct1 = mean(pct.1),
              SD_pct1 = sd(pct.1),
              mean_pct2 = mean(pct.2),
              SD_pct2 = sd(pct.2))

  de_genes_count_intermediate <- merge(genes_accumulated,genes_clusters,by="gene")
  de_genes_counts_final <- merge(de_genes_count_intermediate,genes_clusters_means)

  de_genes_counts_final <- de_genes_counts_final %>%
    arrange(Clusters_DE,desc(mean_avgLogFC))

  ## Write table for all differentially expressed genes containing testing results
  write.table(de_genes_counts_final,
              file=paste(output_dir,"/DE_genes_summary.tsv",sep=""),
              sep="\t",
              quote=FALSE,
              row.names=FALSE,
              col.names=TRUE)

  return(upset_Rlist_DE_genes)

  }
