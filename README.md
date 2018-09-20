# runibic: UniBic biclustering algorithm for R

This package contains implementation of UniBic biclustering algorithm for gene expression data [Wang2016]
The algorithm tries to locate trend-preserving biclusters within complex and noisy data.

## Functions
This package provides the following main functions:
* `BCUnibic`/`runibic` - parallel UniBic for continuous data
* `BCUnibicD` - parallel UniBic for discrete data

The package provides some additional functions:
* `pairwiseLCS` - calculates Longest Common Subsequence (LCS) between two vectors
* `calculateLCS` - calculates LCSes between all pairs of the input dataset
* `backtrackLCS` - recovers LCS from the dynamic programming matrix
* `cluster` - main part of UniBic algorithm (biclusters seeding and expanding)
* `unisort` - returns matrix of indexes based on the increasing order in each row
* `discretize` - performs discretization using Fibonacci heap (sorting method used originally in UniBic) or standard sorting


## Installation
The package may be installed as follows:
```r
install.packages("devtools")
devtools::install_github("athril/runibic")
```

## Example
### Gene expression dataset
This example presents how to use runibic package on gene expression dataset:
```r
library(runibic)
library(biclust)
data(BicatYeast)
res <- biclust(method=BCUnibic(),BicatYeast)
drawHeatmap(BicatYeast, res, 1)
parallelCoordinates(BicatYeast,res,1)
```

### Summarized experiment
This example presents how to use runibic package on SummarizedExperiment:
```r
library(runibic)
library(biclust)
library(SummarizedExperiment)
data(airway, package="airway")
se <- airway[1:20,]
res<- runibic(se)
parallelCoordinates(assays(se)[[1]], res[[1]], 2)
```

## Tutorial
Please check [runibic tutorial](https://github.com/athril/runibic/tree/master/vignettes/runibic.Rmd)

## Citation
For the original sequential version of the UniBic please use the following citation:

**Zhenjia Wang, Guojun Li, Robert W. Robinson, Xiuzhen Huang
*UniBic: Sequential row-based biclustering algorithm for analysis of gene expression data*
Scientific Reports 6, 2016; 23466, doi: https://doi:10.1038/srep23466**


If you use in your work this package with parallel version of UniBic please use the following citation:

**Patryk Orzechowski, Artur Pańszczyk, Xiuzhen Huang Jason H. Moore:
*runibic: a Bioconductor package for parallel row-based biclustering of gene expression data*
bioRxiv, 2017; 210682, doi: [https://doi.org/10.1101/210682](https://doi.org/10.1101/210682)**

BibTex entry:
```
@article{orzechowski2018runibic,
  author = {Orzechowski, Patryk and Pańszczyk, Artur and Huang, Xiuzhen and Moore, Jason H},
  title = {runibic: a Bioconductor package for parallel row-based biclustering of gene expression data},
  journal = {Bioinformatics},
  volume = {},
  number = {},
  pages = {bty512},
  year = {2018},
  doi = {10.1093/bioinformatics/bty512},
  URL = {http://dx.doi.org/10.1093/bioinformatics/bty512},
  eprint = {/oup/backfile/content_public/journal/bioinformatics/pap/10.1093_bioinformatics_bty512/4/bty512.pdf}
}
```


## References
* [Wang2016] Wang, Zhenjia, et al. "UniBic: Sequential row-based biclustering algorithm for analysis of gene expression data." Scientific reports 6 (2016): 23466.
