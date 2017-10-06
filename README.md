# runibic: UniBic biclustering algorithm for R

This package contains implementation of UniBic biclustering algorithm for gene expression data [Wang2016]
The algorithm tries to locate trend-preserving biclusters within complex and noisy data.

## Functions
This package provides the following functions
* BCUnibic - UniBic biclustering algorithm for continuous data
* BCUnibicD - UniBic biclustering algorithm for discrete data
* pairwiseLCS - calculates Longest Common Subsequence (LCS) between two vectors
* calculateLCS - calculates multiple LCS within the dataset and prepares input for the cluster function
* backtrackLCS - recovers the LCS from the matrix obtained using dynamic programming
* cluster - seeds rows based based on the results obtained from calculateLCS
* unisort - returns matrix of indexes based on the increasing order in each row
* discretize - performs discretization using Fibonacci heap (sorting method used originally in UniBic)


## Installation
The package may be installed as follows:
```r
install.packages("devtools")
devtools::install_github("athril/runibic")
```

## Example
This example presents how to use runibic package on gene expression dataset:
```r
library(runibic)
library(biclust)
data(BicatYeast)
input=BicatYeast[1:20,1:20]
biclust(method=BCUnibic(),input)
```


## References
* [Wang2016] Wang, Zhenjia, et al. "UniBic: Sequential row-based biclustering algorithm for analysis of gene expression data." Scientific reports 6 (2016): 23466.
