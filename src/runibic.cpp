#include <iostream>
#include <Rcpp.h>
#include <cstdlib>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <iterator>
#include <functional>
#include "GlobalDefs.h"


using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// Enable OpenMP (exclude macOS)
// [[Rcpp::plugins(openmp)]]




int gTFindex; // Index EOF?
int gColWidth; // TODO: check usage of this option
int gDivided;// TODO: check usage of this option
Rcpp::List fromBlocks(BicBlock ** blocks, const int numBlocks, const int nr, const int nc);




//' Discretize an input matrix
//'
//' This function discretizes the input matrix
//' TODO: Implement functionality
//'
//' @param x a numeric matrix
//' @return a discretized matrix containing integers only
//'
// [[Rcpp::export]]
Rcpp::IntegerMatrix discretize(Rcpp::NumericMatrix x) {
  IntegerMatrix y(x.nrow(),x.ncol());

  gDivided = x.ncol();// TODO: check usage of this option
  
  if(gQuantile >=0.5){
    for(auto iRow = 0; iRow < x.nrow(); iRow++){
      NumericVector rowData = x(iRow,_);
      sort(rowData.begin(), rowData.end());
  
      for(auto iCol = 0; iCol < x.ncol(); iCol++){
        double dSpace = 1.0 / gDivided;
        for(auto ind=0; ind < gDivided; ind++){
          if(x(iRow,iCol) >= calculateQuantile(rowData, x.ncol(), 1.0 - dSpace * (ind+1))){
            y(iRow,iCol) = ind+1;
            break;
          }
        }
      }
    }
  }
  else{
    for(auto iRow = 0; iRow < x.nrow(); iRow++){
      NumericVector rowData = x(iRow,_);
      sort(rowData.begin(), rowData.end());

      double partOne = calculateQuantile(rowData,x.ncol(),1-gQuantile);
      double partTwo = calculateQuantile(rowData,x.ncol(),gQuantile);
      double partThree = calculateQuantile(rowData, x.ncol(), 0.5);
      double upperLimit, lowerLimit;
      
      if((partOne-partThree) >= (partThree - partTwo)){
        upperLimit = 2*partThree - partTwo;
        lowerLimit = partTwo;
      }
      else{
        upperLimit = partOne;
        lowerLimit = 2*partThree - partOne;
      }
      NumericVector biggerPart, lowerPart;
      biggerPart = rowData[rowData > upperLimit];
      lowerPart = rowData[rowData < lowerLimit];
      for(auto iCol = 0; iCol < x.ncol(); iCol++){
        double dSpace = 1.0 / gDivided;
        for(auto ind=0; ind < gDivided; ind++){
          if(lowerPart.size() > 0 && x(iRow,iCol) <= calculateQuantile(lowerPart, lowerPart.size(), dSpace * (ind))){
            y(iRow,iCol) = -ind-1;
            break;
          }
          if(biggerPart.size() > 0 && x(iRow,iCol) >= calculateQuantile(biggerPart, biggerPart.size(), 1.0 - dSpace * (ind+1))){
            y(iRow,iCol) = ind+1;
            break;
          }
        }
      }
    }
  }
  return y;
}
//' Computing the indexes of j-th smallest values of each row
//'
//' This function sorts separately each row of a numeric matrix and returns a matrix
//' in which the value in i-th row and j-th column represent the index of the j-th smallest value of the i-th row.
//'
//' @param x a numeric matrix
//' @return a numeric matrix with indexes indicating positions of j-th smallest element in each row
//'
//' @examples
//' unisort(matrix(c(4,3,1,2,5,8,6,7),nrow=2,byrow=TRUE))
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix unisort(Rcpp::NumericMatrix x) {
  int nr = x.nrow();
  int nc = x.ncol();

  NumericMatrix y(nr,nc);
  int max=omp_get_max_threads();
  omp_set_num_threads(max);


  vector< pair<float,int> > a;
  #pragma omp parallel for private(a)
  for (auto  j=0; j<nr; j++) {
    for (auto  i=0; i<nc; i++) {
      a.push_back(std::make_pair(x(j,i),i));
    }

    sort(a.begin(), a.end());
    for (auto  i=0; i<nc; i++) {
      y(j,i)=a[i].second;
    }
    a.clear();
  }
  return y;
}



//' Calculating a matrix of Longest Common Subsequence (LCS) between a pair of numeric vectors
//'
//' This function calculates the matrix with Longest Common Subsequence (LCS)
//' between two numeric vectors.
//'
//' @param x an integer vector
//' @param y an integer vector
//' @return a matrix storing Longest Common Subsequence (LCS)
//'
//' @examples
//' pairwiseLCS(c(1,2,3,4,5),c(1,2,4))
//'
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix pairwiseLCS(Rcpp::IntegerVector x, Rcpp::IntegerVector y) {

  IntegerMatrix c(x.size()+1,y.size()+1);

  for (auto i=0; i<x.size(); i++) {
    c(i,0)=0;
  }

  for (auto j=0; j<y.size(); j++) {
    c(0,j)=0;
  }

  for(auto i=1; i<x.size()+1; i++) {
    for(auto j=1; j<y.size()+1; j++) {
      if(x(i-1) == y(j-1)) {
        c(i,j) = c(i-1,j-1) + 1;
      }
      else {
        c(i,j) = std::max(c(i,j-1),c(i-1,j));
      }
    }
  }
  return c;
}




//' Retrieving from a matrix Longest Common Subsequence (LCS) between a pair of numeric vector.
//'
//' This function retrieves the Longest Common Subsequence (LCS)
//' between two numeric vectors by backtracking the matrix obtained with dynamic programming.
//'
//' @param c an integer numeric matrix prepared using pairwiseLCS()
//' @param x an integer vector
//' @param y an integer vector
//' @return an integer with the length of Longest Common Subsequence (LCS)
//'
//' @examples
//' backtrackLCS( pairwiseLCS(c(1,2,3,4,5),c(1,2,4)), c(1,2,3,4,5),c(1,2,4))
//'
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector backtrackLCS(Rcpp::IntegerMatrix c, Rcpp::IntegerVector x, Rcpp::IntegerVector y) {
  auto index=c(c.nrow()-1,c.ncol()-1);
  auto i=x.size(), j=y.size();
  Rcpp::IntegerVector lcs(index);

  while (i > 0 && j > 0) {
    if (x(i-1) == y(j-1)) {
      lcs(index-1) = x(i-1);
      i--; j--; index--;
    }
    else if (c(i-1,j) > c(i,j-1))
      i--;
    else
      j--;
   }
   // Print the lcs
   //cout << "LCS of " << X << " and " << Y << " is " << lcs;
  return lcs;
}



//' This function calculates all pairwise LCSes within the array.
//'
//' This function computes unique pairwise Longest Common Subsequences within the matrix.
//'
//' @param discreteInput is a matrix
//' @return a list with informa
//'
//' @examples
//' calculateLCS(matrix(c(4,3,1,2,5,8,6,7),nrow=2,byrow=TRUE))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List calculateLCS(Rcpp::IntegerMatrix discreteInput) {
  int size=discreteInput.nrow()*(discreteInput.nrow()-1)/2;
  vector<int> geneA(size);
  vector<int> geneB(size);
  vector<int> lcslen(size);


// TODO: change into parallel version
// there should be 1-level for loop across all combinations of pairs of rows
//  #pragma omp parallel for private(a,b,i,j,res) schedule(dynamic)
//  for ( auto k=0; k<size; k++ ) {
//    auto i = k/discreteInput.nrow(); auto j=k%discreteInput.nrow(); 
  int k=0;
  for (auto i=0; i<discreteInput.nrow(); i++) {
    for (auto j=i+1; j<discreteInput.nrow(); j++) {
      IntegerVector a = discreteInput(i,_);
      IntegerVector b = discreteInput(j,_);
      geneA[k]=i;
      geneB[k]=j;
      IntegerMatrix res=pairwiseLCS(a,b);
      lcslen[k]=res[res.size()-1,res.size()-1];
      k++;
    }
  }


/* TODO: sort according to lcslen. The following sorting doesn't work:
  std::sort( std::begin(geneA),std::end(geneA), [&](int i1, int i2) { return lcslen[i1] > lcslen[i2]; } );
  std::sort( std::begin(geneB),std::end(geneB), [&](int i1, int i2) { return lcslen[i1] > lcslen[i2]; } );
  std::sort( std::begin(lcslen),std::end(lcslen));//, [&](int i1, int i2) { return lcslen[i1] > lcslen[i2]; } );
  potential workaround: https://stackoverflow.com/questions/37368787/c-sort-one-vector-based-on-another-one
*/
  return List::create(
           Named("a") = geneA,
           Named("b") = geneB,
           Named("lcslen") = lcslen);
}




//' Calculating biclusters from sorted list of LCS scores
//'
//' TODO: Make better parameters
//'
//' @param discreteInput an integer matrix
//' @param scores a numeric vector
//' @param geneOne a numeric vector
//' @param geneTwo a numeric vector
//' @param rowNumber a int with number of rows
//' @param colNumber a int with number of columns
//' @return a number of found clusters
//'
//' @examples
//' cluster(matrix(c(4,3,1,2,5,8,6,7,9,10,11,12),nrow=4,byrow=TRUE),c(13,12,11,7,5,3),c(0,1,2,0,0,1), c(3,2,3,2,1,3),4,3)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List cluster(Rcpp::IntegerMatrix discreteInput, Rcpp::IntegerVector scores, Rcpp::IntegerVector geneOne, Rcpp::IntegerVector geneTwo, int rowNumber, int colNumber) {
  gTFindex = rowNumber-1; // Index EOF?
  gColWidth = max(3+floor(colNumber/30),4.0); // TODO: check usage of this option
  gDivided = colNumber;// TODO: check usage of this option

  int block_id = 0;
  int cnt = 0;
  vector<int> discreteInputData = as<vector<int> >(discreteInput);
  BicBlock** arrBlocks = new BicBlock*[gSchBlock];
  for(auto ind =0; ind<gSchBlock; ind++)
    arrBlocks[ind] = NULL;
  BicBlock *currBlock; // bicluster candidate

  vector<int> vecGenes, vecScores, vecBicGenes, vecAllInCluster; // helpful vectors/stacks

  int components;

  //Memory allocation
  int *colsStat = new int[colNumber]; // column statisctics array    
  long double *pvalues = new long double[rowNumber]; // pvalues array

  bool *candidates = new bool[rowNumber];
  short *lcsLength = new short[rowNumber];
  char** lcsTags = new char*[rowNumber];

  for(auto ind = 0; ind < rowNumber; ind++)  {
    lcsTags[ind] = new char[colNumber];
  }

  //Main loop
  for(auto ind = 0; ind < scores.size(); ind++) {
    /* check if both genes already enumerated in previous blocks */
    bool flag = TRUE;
    /* speed up the program if the rows bigger than 200 */
    if (rowNumber > 250) {
      auto result1 = find(vecAllInCluster.begin(), vecAllInCluster.end(), geneOne(ind));
      auto result2 = find(vecAllInCluster.begin(), vecAllInCluster.end(), geneTwo(ind));

      if ( result1 != vecAllInCluster.end() && result2 != vecAllInCluster.end())
        flag = FALSE;
      else if (gIsTFname && (geneOne(ind) !=  gTFindex) && (geneTwo(ind)!= gTFindex))
        flag = FALSE;
    }
    else {
      flag = check_seed(scores(ind),geneOne(ind), geneTwo(ind), arrBlocks, block_id, rowNumber);
      if (gIsTFname && (geneOne(ind) !=  gTFindex) && (geneTwo(ind)!= gTFindex))
        flag = FALSE;
    }
    if (!flag)    {
      continue;
    }
   
      
    //Init Current block
    currBlock = new BicBlock();
    currBlock->score = min(2, (int)scores(ind));
    currBlock->pvalue = 1;
    vecGenes.clear();
    vecScores.clear();
    //Init vectors/stack for genes and scores
    vecGenes.push_back(geneOne(ind));
    vecGenes.push_back(geneTwo(ind));

    vecScores.push_back(1);
    vecScores.push_back(currBlock->score);
    
    /* branch-and-cut condition for seed expansion */
    int candThreshold = floor(gColWidth * gTolerance);
    if (candThreshold < 2) 
      candThreshold = 2;
    /* maintain a candidate list to avoid looping through all rows */		
    for (auto j = 0; j < rowNumber; j++) 
      candidates[j] = TRUE;
    candidates[(int)geneOne(ind)] = candidates[(int)geneTwo(ind)] = FALSE;
    components = 2;
    /* expansion step, generate a bicluster without noise */
    block_init(scores(ind), geneOne(ind), geneTwo(ind), currBlock, &vecGenes, &vecScores, candidates, candThreshold, &components, &vecAllInCluster, pvalues, rowNumber, colNumber, lcsLength, lcsTags, &discreteInputData);
    /* track back to find the genes by which we get the best score*/       
    
    int  k=0;
    for(k = 0; k < components; k++) {
      if (gIsPValue)
        if ((pvalues[k] == currBlock->pvalue) &&(k >= 2) &&(vecScores[k]!=vecScores[k+1])) 
          break;
      if ((vecScores[k] == currBlock->score)&&(vecScores[k+1]!= currBlock->score)) 
        break;
    }
    components = k + 1;
    if(components > vecGenes.size())
      components = vecGenes.size();    
    for (auto ki=0; ki < rowNumber; ki++) {
      candidates[ki] = TRUE;
    }

    vecGenes.resize(components);    
    for (auto ki=0; ki < vecGenes.size() ; ki++) {
      candidates[vecGenes[ki]] = FALSE;
    }
   // candidates[vecGenes[k]] = FALSE;
    
    bool *colcand = new bool[colNumber];
    for(auto ki = 0; ki < colNumber; ki++)
      colcand[ki] = FALSE;

    // get init block 
    int threshold = floor(components * 0.7)-1;
    if(threshold <1)
      threshold=1;
    //get the statistical results of each column produced by seed
    char *temptag = new char[colNumber];

    for(auto i=0;i<colNumber;i++) {
      colsStat[i] = 0;
      temptag[i] = 0;
    }
    for(auto i=1;i<components;i++) {
      getGenesFullLCS(&discreteInputData[vecGenes[0]*colNumber], &discreteInputData[vecGenes[i]*colNumber],temptag, NULL, colNumber);
      for(auto j=0;j<colNumber;j++)
      {      
        if(temptag[j]!=0)
          colsStat[j]++;
        temptag[j]=0;
      }
    }
    cnt = 0;
    for(auto i=0;i<colNumber;i++) {
      if (colsStat[i] >= threshold) {
        colcand[i] = TRUE;
        cnt++;
      }
    }
    delete[] temptag;
    // add some new possible genes
    int m_ct=0;
    bool colChose = TRUE;
    for(auto ki=0;ki < rowNumber;ki++) {
      colChose=TRUE;
      //if(!candidates[ki]) //if ((po->IS_list && !sublist[ki]) || !candidates[ki]) TODO Sublist check;
      //  continue;
      m_ct=0;
      for (auto i=0; i< colNumber; i++) {
        if (colcand[i] && lcsTags[ki][i]!=0) 
        m_ct++;
      }
      if (candidates[ki]&& (m_ct >= floor(cnt * gTolerance)-1)) {
        int temp;
        for(temp=0;temp<colNumber;temp++) {
          if(colcand[temp]) {
            int tmpcount = colsStat[temp];
            if(lcsTags[ki][temp]!=0)
              tmpcount++;
            if(tmpcount < floor(components * 0.1)-1) {
              colChose = FALSE;
              break;
            }
          }
        }
        if(colChose==TRUE) {
          vecGenes.push_back(ki);
          components++;
          candidates[ki] = FALSE;
          for(temp=0;temp<colNumber;temp++)
          {
            if(lcsTags[ki][temp]!=0 && colcand[ki]) {
              colsStat[temp]++;
            }
          }
        }
      }
    }
    currBlock->block_rows_pre = components;
    // add genes that negative regulated to the consensus 
    char * reveTag;
    for (auto ki = 0; ki < rowNumber; ki++) {
      colChose=TRUE;
     // if (!candidates[ki]) *always false in our case
      //  continue;
      reveTag = new char[colNumber];
      for(auto i = 0; i < colNumber; i++)
        reveTag[i] = 0;
      int commonCnt=0;
      for (auto i=0;i<colNumber;i++) {
        if (discreteInputData[vecGenes[0]*colNumber+i] * (discreteInputData[ki*colNumber+i]) != 0)
          commonCnt++;
      }
      if(commonCnt< floor(cnt * gTolerance)) {
        candidates[ki] = FALSE;
        continue;
      }
      getGenesFullLCS(&discreteInputData[vecGenes[0]*colNumber],&discreteInputData[ki*colNumber],reveTag,lcsTags[vecGenes[1]],colNumber, true);
      m_ct = 0;
      for (auto i=0; i< colNumber; i++) {
        if (colcand[i] && reveTag[i]!=0)
          m_ct++;
      }
      if (candidates[ki] && (m_ct >= floor(cnt * gTolerance)-1)) {
     
        int temp;
        for(temp=0;temp<colNumber;temp++) {
          if(colcand[temp]) {
            int tmpcount = colsStat[temp];
            if(reveTag[temp]!=0)
              tmpcount++;
            if(tmpcount < floor(components * 0.1)-1) {
              colChose = FALSE;
              break;
            }
          }
        }
        if(colChose == TRUE) {
          vecGenes.push_back(ki);
          components++;
          candidates[ki] = FALSE;
          for(temp=0;temp<colNumber;temp++) {
            if(reveTag[temp]!=0 && colcand[ki]) {
              colsStat[temp]++;
            }
          }
        }
      }
      delete[] reveTag;
    }
    // save the current cluster
    for (auto ki = 0; ki < currBlock->block_rows_pre; ki++)
      vecBicGenes.push_back(vecGenes[ki]);
    /* store gene arrays inside block */
    //currBlock->genes = dsNew(components);
    // currBlock->conds = dsNew(cols);

    for (auto j = 0; j < colNumber; j++) {
      if (colcand[j]==TRUE) {
        currBlock->conds.push_back(j);
      }
    }
    currBlock->block_cols = currBlock->conds.size();

    if (currBlock->block_cols < 4 || components < 5){
      delete currBlock;
      continue;      
    }
    currBlock->block_rows = components;
    if (gIsPValue)
      currBlock->score = -(100*log(currBlock->pvalue));
    else
      currBlock->score = currBlock->block_rows * currBlock->block_cols;

    currBlock->genes.clear();    
    for (auto ki=0; ki < components; ki++){
      currBlock->genes.push_back(vecGenes[ki]);
    }
     

    for(auto ki = 0; ki < components; ki++) {
      auto result1 = find(vecAllInCluster.begin(), vecAllInCluster.end(), vecGenes[ki]);
      if(result1==vecAllInCluster.end())
        vecAllInCluster.push_back(vecGenes[ki]);
    }
    /*save the current block b to the block list bb so that we can sort the blocks by their score*/
    arrBlocks[block_id++] = currBlock;

    /* reaching the results number limit */
    if (block_id == gSchBlock) 
      break;
    delete[] colcand;
  }

  // Sorting and postprocessing of biclusters!

  sort(arrBlocks, arrBlocks+block_id, &blockComp);
  
  int n = min(block_id, gRptBlock);
  bool flag;

  BicBlock **output = new BicBlock*[n]; // Array with filtered biclusters
  BicBlock **bb_ptr = output;
  BicBlock *b_ptr;

  double cur_rows, cur_cols;
  double inter_rows, inter_cols;

  /* the major post-processing here, filter overlapping blocks*/
  int i = 0, j = 0, k=0;
  while (i < block_id && j < n) {
    b_ptr = arrBlocks[i];
    cur_rows = b_ptr->block_rows;
    cur_cols = b_ptr->block_cols;
    
    flag = TRUE;
    k = 0;
    while (k < j) {
      inter_rows =0;
      for(auto iter = output[k]->genes.begin(); iter != output[k]->genes.end(); iter++) {
        auto result1 = find(b_ptr->genes.begin(), b_ptr->genes.end(), *iter);
        if(result1!=b_ptr->genes.end()){
          inter_rows++;
        }
      }
      inter_cols=0;
      for(auto iter = output[k]->conds.begin(); iter != output[k]->conds.end(); iter++) {
        auto result1 = find(b_ptr->conds.begin(), b_ptr->conds.end(), *iter);
        if(result1!=b_ptr->conds.end()){
          inter_cols++;
        }
      }
      if (inter_rows*inter_cols > gFilter*cur_rows*cur_cols) {
        flag = FALSE;
        break;
      }
      k++;
    }
    i++;
    if (flag) {
      // print_bc(fw, b_ptr, j++); file print
      j++;
      *bb_ptr++ = b_ptr;
    }
  }
  List outList = fromBlocks(output, j, rowNumber, colNumber);
  for(auto ind =0; ind<gSchBlock; ind++)
    if(arrBlocks!=NULL)
      delete arrBlocks[ind];
  delete[] arrBlocks;
  delete[] colsStat;
  delete[] lcsLength;
  delete[] candidates;
  delete[] pvalues;
  for(auto ind = 0; ind < rowNumber; ind++) {
    delete[] lcsTags[ind];
  }
  delete[] lcsTags;
  delete[] output;


  return outList;
}
Rcpp::List fromBlocks(BicBlock ** blocks, const int numBlocks, const int nr, const int nc) {

  auto x = LogicalMatrix(nr, numBlocks);
  auto y = LogicalMatrix(numBlocks, nc);
  for (int i = 0; i < numBlocks; i++) {
    for (auto it = blocks[i]->genes.begin(); it != blocks[i]->genes.end(); ++it) 
      x(*it, i) = true;
    for (auto it = blocks[i]->conds.begin(); it != blocks[i]->conds.end(); ++it)   
      y(i, *it) = true;
  }
  return List::create(
           Named("RowxNumber") = x,
           Named("NumberxCol") = y,
           Named("Number") = numBlocks,
           Named("info") = List::create());
}



/*
// Unoptimized version of unisort
// [[Rcpp::export]]

Rcpp::NumericMatrix unisort_not_optimal(Rcpp::NumericMatrix x){
  NumericMatrix y(x);
  int nr = x.nrow();
  int nc = x.ncol();

  vector< pair<float,int> > a;
  for (int j=0; j<nr; j++) {
    for (int i=0; i<nc; i++) {
      a.push_back(std::make_pair(x(j,i),i));
    }

    sort(a.begin(), a.end());
    for (int i=0; i<nc; i++) {
      y(j,i)=a[i].second;
    }
    a.clear();
  }
  return y;
}
*/

