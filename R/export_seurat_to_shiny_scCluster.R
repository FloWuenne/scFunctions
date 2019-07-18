#' Exports data from a Seurat object to the formats required by scCluster_genap2
#' @param seurat_object The Seurat object to be processed.
#' @param seurat_object_version Version of the seurat object. Either v2 or v3. Default = v3
#' @param embedding_used Dimensional reduction technique used.Either 'UMAP' or 'TSNE'. Default = 'UMAP'.
#' @param data_slot Data slot to extract the normalized data from. Either 'RNA' or 'SCT'. Default = 'RNA'.
#' @param export_dir Directory to store the exported files into.
#' @keywords Seurat, shiny, vusialization, processing, thresholds
#' @export
#'

export_seurat_to_shiny_scCluster <- function(seurat_object,
                                             seurat_object_version = "v3",
                                             embedding_used = "UMAP",
                                             data_slot = "RNA",
                                             export_dir = "."){

  ## Load libraries
  require("getopt")
  require("Seurat")
  require("feather")
  require("data.table")
  require("dplyr")


    if(seurat_object_version == "v2"){
      seurat_object <- UpdateSeuratObject(seurat_object)
    }else if(seurat_object_version == "v3"){
      seurat_object <- seurat_object
    }

  ## Check which embedding was used
  if(embedding_used == "UMAP"){
    cell_embeddings <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
  } else if (embedding_used == "TSNE"){
    cell_embeddings <- as.data.frame(seurat_object@reductions$tsne@cell.embeddings)
  }

  ## Extract data from Seurat object
  metadata <- seurat_object@meta.data
  metadata$cell_classification <- seurat_object@active.ident

  if(data_slot == "RNA"){
    norm_data <- t(as.data.frame(as.matrix(seurat_object@assays$RNA@data)))
  }else if(data_slot == "SCT"){
    norm_data <- t(as.data.frame(as.matrix(seurat_object@assays$SCT@data)))
  }

  cell_embeddings_with_expression <- merge(cell_embeddings,metadata,by=0)
  rownames(cell_embeddings_with_expression) <- cell_embeddings_with_expression$Row.names
  cell_embeddings_with_expression <- cell_embeddings_with_expression[2:ncol(cell_embeddings_with_expression)]
  cell_embeddings_with_expression <- merge(cell_embeddings_with_expression,norm_data,by=0)

  ## make a tibble for clustering solutiosn
  clustering_solutions <- cell_embeddings_with_expression %>%
    select(cell_classification)

  ## get gene names
  if(data_slot == "RNA"){
    gene_names <- rownames(seurat_object@assays$RNA@data)
  }else if(data_slot == "SCT"){
    gene_names <- rownames(seurat_object@assays$SCT@data)
  }
  gene_names_df <- data.frame("genes" = gene_names)

  ## Format sparse matrix for presto marker calculation

  #Write files
  cell_embeddings_with_expression_transposed_sparse <- cell_embeddings_with_expression[,gene_names]
  cell_embeddings_with_expression_transposed_sparse <- as(t(cell_embeddings_with_expression_transposed_sparse), "sparseMatrix")

  ## Write file containing embeddings and gene expression
  write_feather(cell_embeddings_with_expression,
                path = paste(export_dir,"shiny_clustering_file.feather",sep="/"))

  ## Write file containing multiple annotations for clustering solutions
  write_feather(clustering_solutions,
                path = paste(export_dir,"shiny_user_clustering.feather",sep="/"))

  ## Write gene names
  fwrite(gene_names_df,
         file = paste(export_dir,"shiny_gene_names.tsv",sep="/"))

  ## Write marker table
  saveRDS(cell_embeddings_with_expression_transposed_sparse,
         file = paste(export_dir,"shiny_clustering.sparse_presto.rds",sep="/"),
         version = "2") ## This depends on the R kernel that is running on the galaxy and on the singularity container!

  cat("\n Successfully transformed data! \n")


}

