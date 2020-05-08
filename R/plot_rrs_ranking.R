#' Calculates Regulon specificity score (RSS) from binary regulon activity.
#'
#' @param rrs_df Data frame containing RSS scores for all regulons over all cell types. Can be created with calculate_rrs.
#' @param cell_type Cell type for which to plot jsd ranking. Select "all" to plot a facet plot over all cell types.
#' @param ggrepel_force same as the force parameter for geom_text_repel.
#' @param ggrepel_point_padding same as the force parameter for geom_text_repel.
#' @param top_genes Number of top genes to label in plot using ggrepel.
#' @keywords SCENIC, regulons, RRS, cell type classification
#' @import ggrepel
#' @import cowplot
#' @import tidyverse
#' @export
#' @examples
#'

## Plot JSD enrichment plot for specific cell type
plot_rrs_ranking <- function(rrs_df,
                             cell_type,
                             ggrepel_force = 1,
                             ggrepel_point_padding = 0.2,
                             top_genes = 4,
                             plot_extended = FALSE){

  if(plot_extended == TRUE){
    rrs_df <- rrs_df %>%
      subset(grepl("extended",regulon))
  }else if(plot_extended == FALSE){
    rrs_df <- rrs_df %>%
      subset(!grepl("extended",regulon))
  }

  if(cell_type != "all") {
    rrs_df_sub <- rrs_df %>%
      subset(cell_type == cell_type) %>%
      arrange(desc(RSS))

    rrs_df_sub <- rrs_df_sub %>%
      mutate("rank" = as.numeric(rownames(rrs_df_sub)))

    #jsd_matrix_sub$regulon <- factor(jsd_matrix_sub$regulon,levels = unique(jsd_matrix_sub$regulon))

    rrs_ranking_plot <- ggplot(rrs_df_sub,aes(rank,RSS,label = regulon)) +
      geom_point(color = "grey20",size = 2) +
      geom_point(data = subset(rrs_df_sub,rank < top_genes),
                 color = "red",size = 2) +
      geom_text_repel(data = subset(rrs_df_sub,rank < top_genes),
                      force = ggrepel_force,point.padding = ggrepel_point_padding) +
      labs(x = "Rank",
           y = "RSS",
           title = cell_type)
  }else if(cell_type == "all"){

    rrs_df_sub <- rrs_df %>%
      group_by(cell_type) %>%
      mutate("rank" = order(order(RSS, decreasing = TRUE)))

    #jsd_matrix_sub$regulon <- factor(jsd_matrix_sub$regulon,levels = unique(jsd_matrix_sub$regulon))

    rrs_ranking_plot <- ggplot(rrs_df_sub,aes(rank,RSS,label = regulon)) +
      geom_point(color = "grey20",size = 2) +
      geom_point(data = subset(rrs_df_sub,rank < top_genes),
                 color = "red",size = 2) +
      geom_text_repel(data = subset(rrs_df_sub,rank < top_genes),
                      force = ggrepel_force,point.padding = ggrepel_point_padding) +
      labs(x = "Rank",
           y = "RSS",
           title = cell_type) +
      facet_wrap(~ cell_type)
  }


  return(rrs_ranking_plot)

}
