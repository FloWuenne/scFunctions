# scFunctions
A collection of functions for single cell data analysis by Florian Wuennemann. Mainly operates on Seurat objects. For more information on Seurat and the object structure, see https://satijalab.org/seurat/.

# Installation
Install the package using the devtools install_github function as shown below:

```
library(devtools)
install_github("FloWuenne/scFunctions")
```

# SCENIC functions

You can find a tutorial on how to use the functions in this package to further process and analyze SCENIC results in this [tutorial]):

[Processing and visualization of SCENIC results](./Tutorials/process_SCENIC.Rmd)

The following functions can be used to analyze results from the SCENIC pipeline for gene regulatory networks.

**auc_thresh_kmeans()**

**binarize_regulons()**

**calculate_rrs()**

**plot_rrs_ranking()**

**plot_bin_regulon_UMAP()**
