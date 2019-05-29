Processing and visualization of SCENIC results
================
Florian Wuennemann
March 22, 2005

# Introduction

This tutorial will describe how to use the functions implemented in this
package to further process the output from a typical run of the [SCENIC
pipeline](https://github.com/aertslab/SCENIC). This tutorial assumes
that you have processed your data up until the third step of the
pipeline. The following data will be required for completely running
this tutorial:

  - **regulonAUC** - The regulon activity scores in matrix format
    (output from file 3.4\_regulonAUC.Rds)
  - **cell classifications for your cell** - For example using a
    metadata data frame from a seurat object

We provide a small test data set as part of this package, which can be
found in ./example\_data to help test the scripts and analysis and get
familiar with the different data formats and plots. The example\_data is
from [Wuennemann et
al.](http://andelfingerlab.heart_maturation.genap.ca/) and represents
the AUC values for heart cells from E14.5 embryonic hearts.

If any of the defintions or terms in this tutorial are unclear, please
visit the SCENIC FAQ page and see whether your question is answered
there already:

[SCENIC
FAQ](https://github.com/aertslab/SCENIC/blob/master/vignettes/FAQ.md)

# Installation

You can easily install the package with the following command:

# Regulon analysis

Many of the analysis concepts and statistics that are calculated in this
tutorial have been derived from a publication by [Suo et al. (2018) in
Cell
Reports](https://www.cell.com/cell-reports/fulltext/S2211-1247\(18\)31634-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124718316346%3Fshowall%3Dtrue).
In their publication, they used a modified version of SCENIC to call
regulons for the [Mouse Cell Atlas dataset](http://bis.zju.edu.cn/MCA/).

## Determine AUC thresholds

In the original implemntation of SCENIC, the authors included a function
to determine thresholds for the AUC activity via the function:
AUCell\_exploreThresholds(). In this package, we wrote a very simple
function that determines thresholds based on k-means clustering on the
AUC distribution. This function performs comparable to the original
implementation but is much faster. Please keep in mind however, that
this function might not perform very well for setting thresholds for
non-bimodal distribution, which are quite often observed for regulons.
We still advise to manually check and adjust the thresholds prior to
binarization of regulons\!

Let us use the regulons AUC values to determine thresholds using our
k-means function.

The thresholds are saved in list format where each regulon is the name
of the list and the AUC threshold is the value of the list.

``` r
head(kmeans_thresholds)
```

    ## $`Ybx1_extended (738g)`
    ## [1] 0.2018032
    ## 
    ## $`Ybx1 (219g)`
    ## [1] 0.07257359
    ## 
    ## $`E2f1_extended (1350g)`
    ## [1] 0.1336872
    ## 
    ## $`E2f1 (1289g)`
    ## [1] 0.1264914
    ## 
    ## $`Srebf2_extended (1404g)`
    ## [1] 0.1285305
    ## 
    ## $`Gabpa_extended (2442g)`
    ## [1] 0.2076677

While we use a sample k-means clustering approache with 2 clusters here
to determine thresholds, there are obviously more sophisticated
approaches than this to determine thresholds in the AUC distribution.
Feel free to develop your own approaches to determine optimal thresholds
to binarize regulons. I would be excited to hear back if you developed
your own function to perform this task\!

## Binarize regulons using thresholds

Now that we have our thresholds, it is time to binarize the regulons
using these thresholds.

Let’s take a look at the first regulon in the binary regulon list.

``` r
head(binary_regulons$`Ybx1_extended (738g)`)
```

    ##                   cells Ybx1_extended (738g)
    ## 1 E14.5_r1_TTCTGACGATCT                    1
    ## 2 E14.5_r1_CACAAGACGTTC                    1
    ## 3 E14.5_r1_ACTGGCGATTAT                    1
    ## 4 E14.5_r1_CAATAAGTGAGN                    1
    ## 5 E14.5_r1_ACAGGTCGACGC                    1
    ## 6 E14.5_r1_GATTTCCATACC                    1

Next, we have to reformat the binary regulons into a big data frame that
contains all of the binary regulons so that we can use them to calculate
RRS scores.

Let’s check that the data table is formatted correctly before
proceeding:

``` r
binary_regulons_trans[1:4,1:3]
```

    ##                       E14.5_r1_TTCTGACGATCT E14.5_r1_CACAAGACGTTC
    ## Ybx1_extended (738g)                      1                     1
    ## Ybx1 (219g)                               0                     1
    ## E2f1_extended (1350g)                     1                     1
    ## E2f1 (1289g)                              1                     1
    ##                       E14.5_r1_ACTGGCGATTAT
    ## Ybx1_extended (738g)                      1
    ## Ybx1 (219g)                               1
    ## E2f1_extended (1350g)                     1
    ## E2f1 (1289g)                              1

## Calculate Regulon Specificity Score (RSS)

We now want to use the binary regulon activity together with the cell
assignments to see how specific each predicted regulon is for each cell
type. We can do this by calculating a regulon specificity score (RSS)
which is based on the Jensen-Shannon divergence, a measure of the
similarity between two probability distributions. Basically for the
calculation of the RSS, we will calculate the Jensen-Shannon divergence
between each vector of binary regulon activity overlaps with the
assignment of cells to a specific cell type.

First, we need to load the cell assignments as a data frame. This data
frame needs to have cell names that correspond with the binary regulon
data frame as rownames and contain a column labeled “cell\_type”, which
contains the assignments for all cells. For convenience, you can use the
metadata table from a correspondin Seurat object, just make sure that
you add a column labeled
    “cell\_type”.

``` r
head(metadata_sub)
```

    ##                       orig.ident nCount_RNA nFeature_RNA Mouse_ID   Age
    ## E14.5_r1_TTCTGACGATCT      E14.5      10497         3394       M1 E14.5
    ## E14.5_r1_CACAAGACGTTC      E14.5       9826         3342       M1 E14.5
    ## E14.5_r1_ACTGGCGATTAT      E14.5      10221         3266       M1 E14.5
    ## E14.5_r1_CAATAAGTGAGN      E14.5       9258         3450       M1 E14.5
    ## E14.5_r1_ACAGGTCGACGC      E14.5      10885         3579       M1 E14.5
    ## E14.5_r1_GATTTCCATACC      E14.5      10459         3578       M1 E14.5
    ##                       Replicate Litter Illumina_flowcell Sex
    ## E14.5_r1_TTCTGACGATCT        r1     L1         HLLVHBGXY   M
    ## E14.5_r1_CACAAGACGTTC        r1     L1         HLLVHBGXY   M
    ## E14.5_r1_ACTGGCGATTAT        r1     L1         HLLVHBGXY   M
    ## E14.5_r1_CAATAAGTGAGN        r1     L1         HLLVHBGXY   M
    ## E14.5_r1_ACAGGTCGACGC        r1     L1         HLLVHBGXY   M
    ## E14.5_r1_GATTTCCATACC        r1     L1         HLLVHBGXY   M
    ##                       doublet_scores_scrublet predicted_doublet_scrublet
    ## E14.5_r1_TTCTGACGATCT              0.02802721                      FALSE
    ## E14.5_r1_CACAAGACGTTC              0.03094556                      FALSE
    ## E14.5_r1_ACTGGCGATTAT              0.02165544                      FALSE
    ## E14.5_r1_CAATAAGTGAGN              0.04369102                      FALSE
    ## E14.5_r1_ACAGGTCGACGC              0.03033740                      FALSE
    ## E14.5_r1_GATTTCCATACC              0.02253797                      FALSE
    ##                       percent.mito    S.Score  G2M.Score Phase nCount_SCT
    ## E14.5_r1_TTCTGACGATCT   0.02438792 -0.6667862 -1.3860525    G1       2530
    ## E14.5_r1_CACAAGACGTTC   0.01465500 -0.4816747 -0.7297447    G1       2458
    ## E14.5_r1_ACTGGCGATTAT   0.03238431 -0.4856854 -1.2957825    G1       2502
    ## E14.5_r1_CAATAAGTGAGN   0.03262044  0.2346583  2.3575657   G2M       2464
    ## E14.5_r1_ACAGGTCGACGC   0.02406982 -0.7266369 -1.2373289    G1       2497
    ## E14.5_r1_GATTTCCATACC   0.01960034 -0.2892592 -1.3102109    G1       2455
    ##                       nFeature_SCT SCT_snn_res.0.8 cell_type
    ## E14.5_r1_TTCTGACGATCT         1345               7      ACM1
    ## E14.5_r1_CACAAGACGTTC         1356               7      ACM1
    ## E14.5_r1_ACTGGCGATTAT         1258               7      ACM1
    ## E14.5_r1_CAATAAGTGAGN         1429               7      ACM1
    ## E14.5_r1_ACAGGTCGACGC         1402               4      VCM1
    ## E14.5_r1_GATTTCCATACC         1352               4      VCM1
    ##                       clustering_lvl2 super_group
    ## E14.5_r1_TTCTGACGATCT            ACM1          CM
    ## E14.5_r1_CACAAGACGTTC            ACM1          CM
    ## E14.5_r1_ACTGGCGATTAT            ACM1          CM
    ## E14.5_r1_CAATAAGTGAGN            ACM1          CM
    ## E14.5_r1_ACAGGTCGACGC            VCM1          CM
    ## E14.5_r1_GATTTCCATACC            VCM1          CM

Now that we are ready to calculate the RSS for all regulons over all
cell types.

The output is a data frame with a RSS score for each regulon - cell type
combination.

``` r
head(rrs_df)
```

    ##                 regulon cell_type       RSS
    ## 1          Ikzf1 (638g)      RMPH 1.0000000
    ## 2 Ikzf1_extended (882g)      RMPH 0.8211780
    ## 3 Runx3_extended (204g)      RMPH 0.8151686
    ## 4           Foxc1 (58g)     ENDO1 0.7794728
    ## 5   Fev_extended (195g)     ENDO1 0.7532211
    ## 6   Mycn_extended (17g)      VCM1 0.7044754

We can visualize the RSS by performing ranking on the RSS scores with
the most specific regulons ranking the highest per cell type. I have
included a function (plot\_rrs\_ranking) to easily plot an RSS ranking
plot from this data frame. The function has a couple of options, most of
which are cosmetic. Importantly, you can either plot the RSS ranking for
a cell type of interest or you can set cell\_type = “all” to plot the
RSS over all cell types. plot\_extended determines whether you would
like to plot the high confidence regulons only or if you want to plot
the regulons named \_extended, which also contain genes only based on
motif prediction based on similarity.

``` r
plot_rrs_ranking(rrs_df,
                 "RMPH",
                 ggrepel_force = 1,
                 ggrepel_point_padding = 0.2,
                 top_genes = 4,
                 plot_extended = FALSE)
```

    ## Loading required package: ggrepel

    ## Loading required package: cowplot

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

![](process_SCENIC_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

We can also easily visualize all regulons over all cell types using
heatmaps. Let’s first investigate the distribution of RSS over all cell
types.

``` r
library(ggridges)
```

    ## 
    ## Attaching package: 'ggridges'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     scale_discrete_manual

``` r
rrs_df_nona <- subset(rrs_df,RSS > 0)
ggplot(rrs_df_nona,aes(RSS,cell_type, fill = cell_type)) +
  geom_density_ridges(scale = 5, alpha = 0.75) +
  geom_vline(xintercept = 0.1) +
  theme(legend.position = "none")
```

    ## Picking joint bandwidth of 0.0137

![](process_SCENIC_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

The RSS distribution clearly shows that the RSS is highly dependent upon
the cell type we are investigating. As we can see, resident macrophages
show very high and specific RSS, while other cell types for which more
similar cell types exist in the dataset, like cardiomyocytes show less
specificty for the regulons. In this small toy example, it seems that ~
0.05 - 0.1 will capture specific regulons for most cell types. For
highlighting purposes, we are gonna filter with 0.4 to be able to easily
visualize the result in a heatmap.

``` r
rrs_df_wide <- rrs_df %>%
  spread(cell_type,RSS)

rownames(rrs_df_wide) <- rrs_df_wide$regulon 
rrs_df_wide <- rrs_df_wide[,2:ncol(rrs_df_wide)]

## Subset all regulons that don't have at least an RSS of 0.7 for one cell type
rrs_df_wide_specific <- rrs_df_wide[apply(rrs_df_wide,MARGIN = 1 ,FUN =  function(x) any(x > 0.4)),]
```

We can then visualize the regulons that show an RSS over the defined
threshold in this example using heatmapply, a heatmap library using
plotly.

``` r
library(heatmaply)
```

    ## Loading required package: plotly

    ## 
    ## Attaching package: 'plotly'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     last_plot

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## The following object is masked from 'package:graphics':
    ## 
    ##     layout

    ## Loading required package: viridis

    ## Loading required package: viridisLite

    ## Registered S3 method overwritten by 'seriation':
    ##   method         from 
    ##   reorder.hclust gclus

    ## 
    ## ======================
    ## Welcome to heatmaply version 0.16.0
    ## 
    ## Type citation('heatmaply') for how to cite the package.
    ## Type ?heatmaply for the main documentation.
    ## 
    ## The github page is: https://github.com/talgalili/heatmaply/
    ## Please submit your suggestions and bug-reports at: https://github.com/talgalili/heatmaply/issues
    ## Or contact: <tal.galili@gmail.com>
    ## ======================

``` r
heatmaply(rrs_df_wide_specific)
```

![](process_SCENIC_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

This concludes the section about RSS calculations.

## Calculate conecction specificity index (CSI) for all regulons

The final statistics that we want to calculate is the connection
specificty index. The CSI is a major of connectedness between the
different regulons. Regulons that share high CSI likely are
co-regulating downstream genes and are together responsible for cell
function. You can read more about the theoretical details of CSI here:

[GAIN](http://lovelace.cs.umn.edu/similarity_index/guide.php)

We can calculate the CSI scores for all regulon pairs based on the AUCs
matrix for all regulons. Again, this function has a switch to either
select high confidence regulons or run the CSI calculation only on
\_extended regulons (calc\_extended = TRUE).

``` r
regulons_csi <- calculate_csi(regulonAUC)
```
