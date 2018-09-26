#' A function that will perform fgsea analysis on all cells in a given Seurat object.
#'
#'
#' @param seurat_object The S4 Seurat object which contains filtered and normalized cells in the data slot.
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

perform_fgsea <- function(seurat_object,
                          cell_types = sort(unique(seurat_object@ident)),
                          output_dir = ".",
                          pathway_file){

  ## load libraries
  require(fgsea)
  require(Seurat)
  require(gskb)

  ## Define cell types etc

  df_total = NULL
  for (celltypename in cell_types) {
    celltype_seurat_subset <- SubsetData(seurat_object,
                                         ident.use = celltypename)


    #Plot enrichment graph for the cell type (average)
    average_cell_Type <- AverageExpression(celltype_seurat_subset, genes.use = NULL, return.seurat = FALSE,
                                           add.ident = NULL, use.scale = TRUE, use.raw = FALSE,
                                           show.progress = TRUE)

    gene_names <- rownames(average_cell_Type)

    average_cell_Type <- unlist(average_cell_Type, recursive = TRUE, use.names = TRUE)

    names(average_cell_Type) <- gene_names

    #Caculate the enrichment score for each cell and save gsea NES
    celltype_seurat_subset.scaled.data <- celltype_seurat_subset@scale.data

    colnames <- colnames(celltype_seurat_subset@data)

    gsea=NULL
    cell_count <- 0

    ## Load mouse pathways from gskb
    data(mm_pathway)

    for (cellname in colnames) {
      cell_count <- cell_count + 1
      print(cell_count)
      cell <- celltype_seurat_subset.scaled.data[,cellname]
      fgseaRes <- fgsea(pathways = mm_pathway,
                        stats = cell,
                        minSize=15,
                        maxSize=500,
                        nperm=1000,
                        nproc=6)

      fgseaRes <-data.frame(fgseaRes, cell_type=celltypename)
      gsea <- rbind(gsea,fgseaRes)
    }


    df_total <- rbind(df_total,gsea)
  }


  fwrite(df_total, file ="GSEA_GOautophagy.txt", sep="\t")

  df_total <- fread("GSEA_GOautophagy.txt", header=T)


  #Make nice volcano plot

  volcano_plot <- ggplot(df_total,aes(NES,-log(padj), color=cell_type)) + geom_point(size=2) +
    geom_point(
      data = subset(df_total,(NES > 1 ) & (padj < 0.05)), col="red"
    ) +
    ggtitle("Volcano Plot")

  ggsave(volcano_plot, filename="VolcanoPlot.png")




  #Make nice bar plot


  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                        conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
  }


  tgc <- summarySE(df_total, measurevar="NES", groupvars="cell_type")

  tgc2 <- tgc
  tgc2$cell_type <- factor(tgc2$cell_type)


  png("BarPlotNES.Retina_GOautophagy_CTL_SEM_black.png", width=3000, height=2000, res=300)
  ggplot(tgc2, aes(x=reorder(cell_type,-NES), y=NES, fill=NULL)) +
    geom_bar(position=position_dodge(), colour="black", stat="identity") +
    geom_errorbar(aes(ymin=NES-se, ymax=NES+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    xlab("Cell Type") +
    ylab("Mean of Normalized Enrichement Score (± SEM)") +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  dev.off()



  aov.dfb <- aov(NES ~ cell_type, data=df_total)

  posthoc <- TukeyHSD(x=aov.dfb, "cell_type",  conf.level=0.95)


  summary(aov.dfb)

  print(model.tables(aov.dfb,"means"),digits=3)


  posthoc.save <- as.data.frame(posthoc$cell_type)

  write.table(posthoc.save, "Anova_Tukey_Retina_GOautophagy_CTL.txt", sep="\t")


}
