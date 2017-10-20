#####################################################################################################
# This example shows the basic usage of runibic on synthetic data
#####################################################################################################

library(runibic)
# prepare a random matrix
test <- matrix(rnorm(1000), 100, 100)
test

# run UniBic biclustering algorithm 
res <- biclust::biclust(test, method = BCUnibic())

# check how many biclusters were found
res@Number

# inspect tows of the third bicluster
which(res@RowxNumber[,3])

# inspect columns of the third bicluster
which(res@NumberxCol[3,])

# draw the heatmap with the third bicluster
drawHeatmap(test, res, 3)

# show parallel cords of the fourth bicluster
parallelCoordinates(test,res,4)


#####################################################################################################
# This example shows how to use runibic on GDS dataset
#####################################################################################################
#library(GEOquery)
#library(affy)
# download dataset with Rat peripheral and brain regions from Gene Omnibus
gse <- GEOquery::getGEO("GDS589", GSEMatrix = TRUE)

#convert dataset to ExpressionSet
eset <- affy::GDS2eSet(gds)
subset <-affy::exprs(eset)[1:100,]

#perform analysis on first 100 of genes
res <- runibic(subset)

# draw the heatmap with the second bicluster
drawHeatmap(subset, res, 1)



#####################################################################################################
# The following example presents the usage of runibic with data
# from an example dataset from an RNA-Seq experiment (SummarizedExperiment)
#####################################################################################################
library(runibic)
# load airway dataset
data(airway, package="airway")
se <- airway[1:20,]

# run runibic
res<-runibic(se)

# analyze the content of the second bicluster
bicluster(assays(se)[[1]], res[[1]], 2)

# draw the heatmap with the second bicluster
drawHeatmap(assays(se)[[1]], res[[1]], 2)

# show parallel cords of the second bicluster
parallelCoordinates(assays(se)[[1]], res[[1]], 2)



#####################################################################################################
# compare the results of different biclustering algorithms
#####################################################################################################
library(runibic)
data(BicatYeast)

# load package with QUBIC biclustering algorithm for comparison
library(QUBIC)

# perform biclustering using CC, Bimax, Qubic, Plaid and Unibic
resCC <- biclust::biclust(BicatYeast, method = BCCC())
resBi <- biclust::biclust(BicatYeast, method = BCBimax())
resQub <- biclust::biclust(BicatYeast, method = BCQU())
resPlaid <- biclust::biclust(BicatYeast, method = BCPlaid())
resUni <- biclust::biclust(BicatYeast, method = BCUnibic())

# compare the results
QUBIC::showinfo(BicatYeast,c(resCC, resBi, resPlaid, resQub, resUni))



#####################################################################################################
# Examples of finding the Longest Common Subsequence between two vectors
#####################################################################################################

# prepare two vectors for comparison
A=c(1,2,1,5,4,3)
B=c(2,1,3,2,1,4)

# check which values are common for both vectors
backtrackLCS(A,B)

# calculate using dynamic programming the matrix for Longest Common Subsequence
pairwiseLCS(A,B)


#####################################################################################################
# Examples for various other functions
#####################################################################################################

#prepare input data
A = matrix(c(11,17,12,10,8,9,19,15,18,13,14,7,4,6,16,2,3,1,5,20,
             17,1,8,15,5,10,2,12,9,7,3,14,11,4,6,16,20,13,19,18,
             15,8,17,12,18,14,19,11,16,20,10,13,6,3,7,9,1,2,5,4,
             15,12,16,9,19,17,10,18,11,20,8,13,2,5,7,14,1,3,4,6,
             15,10,9,6,13,19,7,18,16,17,14,4,3,1,2,20,12,5,11,8,
             1,7,4,3,2,6,8,13,5,9,12,11,16,15,17,10,19,20,14,18,
             10,5,3,9,2,11,6,13,8,1,7,4,16,14,15,12,18,17,20,19,
             10,5,1,12,8,11,7,13,6,4,3,2,18,14,15,9,17,16,20,19,
             9,6,3,10,1,12,7,13,8,2,5,4,16,14,15,11,19,17,20,18,
             12,8,1,3,2,11,4,14,9,7,10,5,16,13,15,6,18,17,20,19), nrow=10, byrow=TRUE)


#check the matrix with sorted each row
unisort(A)

#calculate the length of LCS between each pair of rows
lcs=calculateLCS(A)

#check the length of the longest common subsequence (LCS between rows 6 and 7 is equal to 10)
list(a=lcs$a[2], b=lcs$b[2], lclen=lcs$lcslen[2])

#discretize a matrix
B=replicate(10, rnorm(10))
B

#check each row
runiDiscretize(B)

#cluster example
A = matrix(c(4,3,1,2,5,8,6,7,9,10,11,12),nrow=4,byrow=TRUE)
iA = unisort(A)
lcsResults = calculateLCS(A)
cluster(iA,A,lcsResults$lcslen,lcsResults$a, lcsResults$b,nrow(A),ncol(A))
