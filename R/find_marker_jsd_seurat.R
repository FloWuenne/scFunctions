#' Uses Jensen-Shannon divergence to find marker genes for clusters
#'
#' @param seurat_object The S4 Seurat object which contains filtered and normalized cells in the data slot.
#' @param custom_ident If specified, will use a column from meta.data instead of the ident of the seurat object.
#' @param percent_det Percentage of cells a gene has to be expressed in to be considered for being tested as a marker. default= 0.25 (10% of all cells)
#' @keywords Seurat, marker genes, Jensen-Shannon divergence
#' @export
#' @examples
#' rm_doublets_seurat()

## dependencies:
## Seurat : https://github.com/satijalab/seurat
## tidyr :
## dplyr :
##


find_marker_jsd_seurat <- function(seurat_object,
                                   custom_ident = FALSE,
                                   percent_det = 0.1){

  ## Load libraries
  require(Seurat)
  require(tidyr)
  require(dplyr)
  require(philentropy)
  require(data.table)

  if(custom_ident != FALSE){
    if(custom_ident %in% colnames(seurat_object)){

    }else{
      print("Please select a valid column from the Seurat meta.data table!")
    }

  }else{

    #### 1) Setting up data
    print("1) Formatting data...")

    ## Get cell cluster identities
    identities <- seurat_object@ident
    identities_df <- data.frame("cell" = as.character(names(identities)),
                                "cluster" = identities)

    ## Reformat gene expression into long format
    gene_exp <- as.data.frame(as.matrix(seurat_object@data))
    gene_exp$gene <- rownames(gene_exp)
    gene_exp_long <- gene_exp %>%
      gather(cell,expression,-gene) %>%
      mutate("cell" = as.character(cell))

    ## Set levels for expression data
    gene_exp_long$cell <- factor(gene_exp_long$cell,levels = levels(identities_df$cell))

    ## merge expression data with cluster annotation
    gene_exp_long_merged <- dplyr::left_join(gene_exp_long,identities_df, by="cell")

    ## Calculate normalized expression per gene such that for each gene, expression sums up to 1
    gene_exp_long_merged_norm <- gene_exp_long_merged %>%
      group_by(gene) %>%
      mutate("norm_expression" = expression / sum(expression) ) %>%
      subset(!is.na(norm_expression))


    #### Filter data based on user preferences
    #### 2) Filtering data
    print("2) Filtering input data...")

    ## Prefilter data set based on user specified thresholds

    ## Count cells that have no expression for all genes
    count_zero <- gene_exp_long_merged_norm %>%
      group_by(gene) %>%
      subset(expression == 0) %>%
      tally()

    ## Count cells that are non-zero for all genes
    count_nonzero <- gene_exp_long_merged_norm %>%
      group_by(gene) %>%
      subset(expression != 0) %>%
      tally()

    ## Calculate fraction of cells with zero expression for all genes
    count_table <- left_join(count_zero,count_nonzero,by = "gene")
    count_table <- count_table %>%
      mutate(frac_det = 1 - (n.x / (n.y + n.x))) %>%
      subset(frac_det >= percent_det)

    ## Subset expression data to only keep genes
    filtered_data <- gene_exp_long_merged_norm %>%
      subset(gene %in% count_table$gene) %>%
      select(gene,cell,cluster,norm_expression)

    filtered_data <- data.table(filtered_data)

    #### 3) Perform marker testing
    print("3) Performing marker testing...")

    ## Create data frame to collect all stats
    nGene <- length(unique(filtered_data$gene))
    nClust <- length(unique(filtered_data$cluster))

    ## Set up final matrix of markers
    all_marker_df <- matrix(nrow = nGene * nClust, ncol=4)

    ## Create counter
    counter <- 0
    gene_counter <- 0

    ## Print status message about marker identification
    print(paste("Testing :",nGene," genes", " for ",nClust," clusters!"),sep="")

    ## Create progressBar
    ## Code for progress bar: https://ryouready.wordpress.com/2009/03/16/r-monitor-function-progress-with-a-progress-bar/
    pb <- txtProgressBar(min = 0, max = (nGene * nClust), style = 3)

    ## Iterate over all genes after filtering
    for(this_gene in unique(filtered_data$gene)) {

      gene_counter <- gene_counter + 1

      ## Iterate over each cluster for this gene
      for(this_cluster in unique(filtered_data$cluster)){

        ## Print current gene
        counter <- counter + 1

        ## print progress Bar
        setTxtProgressBar(pb, counter)

        ## Calculate cluster probablity vector
        sub_data <- filtered_data %>%
          subset(gene == this_gene) %>%
          mutate("cluster_bin" = if_else(cluster == this_cluster,1,0)) %>%
          mutate("cluster_bin" = cluster_bin / sum(cluster_bin)) %>%
          mutate("cluster_bin" = ifelse(is.na(cluster_bin), 0, cluster_bin))

        jsd_df <- rbind(c(sub_data$cluster_bin),
                        c(sub_data$norm_expression))

        ## Calculate Jensen-Shannon divergence
        jsd_divergence <- suppressMessages(philentropy::JSD(jsd_df))

        ## Calculate distance from divergence
        jsd_distance <- 1 - sqrt(jsd_divergence)

        ## Add data to output dataframe
        all_marker_df[counter,1] <- this_gene
        all_marker_df[counter,2] <- this_cluster
        all_marker_df[counter,3] <- jsd_divergence
        all_marker_df[counter,4] <- jsd_distance

      }

    }
    close(pb)

    #### 4) Writing results
    print("4) Writing results...")

    ## Remove any jsd calculations that are NA
    all_marker_df <- as.data.frame(all_marker_df)
    colnames(all_marker_df) <- c("gene","cluster","jsd_divergence","jsd_distance")

    all_marker_df <- subset(all_marker_df,!is.na(jsd_distance))

    all_marker_df$jsd_distance <- as.numeric(as.character(all_marker_df$jsd_distance))

    ## Calculate the top 50 markers based on jsd per cluster
    top_marker_per_cluster <- all_marker_df %>%
      group_by(cluster) %>%
      arrange(desc(jsd_distance)) %>%
      top_n(50) %>%
      arrange(cluster)
  }


  return(top_marker_per_cluster)

}
