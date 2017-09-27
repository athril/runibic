#include <iostream>
#include <Rcpp.h>
#include <cstdlib>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <iterator>
#include "GlobalDefs.h"


using namespace std;
using namespace Rcpp;

//' Computing the indexes of j-th smallest values of each row
//'
//' This function sorts separately each row of a numeric matrix and returns a matrix
//' in which the value in i-th row and j-th column represent the index of the j-th smallest value of the i-th row.
//'
//' @param x a numeric matrix
//' @return a numeric matrix with indexes indicating positions of j-th smallest element in each row
//'
//' @examples
//' A=matrix(c(4,3,1,2,5,8,6,7),nrow=2,byrow=TRUE)
//' unisort(A)
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


//' Calculating Longest Common Subsequence (LCS) between two numeric vectors
//'
//' This function calculates using dynamic programming the Longest Common Subsequence (LCS)
//' between two numeric vectors.
//'
//' @param m a numeric vector
//' @param n a numeric vector
//' @return a matrix storing Longest Common Subsequence (LCS)
//'
//' @examples
//' calculateLCS(c(1,2,3,4,5),c(1,2,4))
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calculateLCS(Rcpp::NumericVector m, Rcpp::NumericVector n) {
  NumericMatrix pc(m.size(), n.size());
  NumericMatrix pb(m.size(), n.size());

  for (auto j=0; j<m.size(); j++) {
    pc(j,0)=0;
    pb(j,0)=0;
  }

  for (auto j=0; j<n.size(); j++) {
    pc(0,j)=0;
    pb(0,j)=0;
  }

  // TODO: implement parallel LCS
  for(auto i=1; i<m.size(); i++) {
    for(auto j=1; j<n.size(); j++) {
      if(m(i-1) == n(j-1)) {
        pc(i,j) = pc(i-1,j-1) + 1;
        pb(i,j) = 1;
      }
      else if(pc(i-1,j) >= pc(i,j-1)) {
        pc(i,j) = pc(i-1,j);
        pb(i,j) = 2;
      }
      else {
        pc(i,j) = pc(i,j-1);
        pb(i,j) = 3;
      }
    }
  }
  return pc;
}

//' Calculating biclusters from sorted list of LCS scores
//'
//' TODO: Make better parameters
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
int cluster(Rcpp::IntegerMatrix discreteInput, Rcpp::IntegerVector scores, Rcpp::IntegerVector geneOne, Rcpp::IntegerVector geneTwo, int rowNumber, int colNumber) {
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
      /*else if (gIsList&&(!sublist[geneOne(ind)] || !sublist[geneTwo(ind)]))    TODO: Check sublist usage in file reading
          flag =FALSE;*/
    }
    else {
      flag = check_seed(scores(ind),geneOne(ind), geneTwo(ind), arrBlocks, block_id, rowNumber);
      Rcout << scores(ind) << " " << geneOne(ind) << " " << geneTwo(ind) << endl;
      if (gIsTFname && (geneOne(ind) !=  gTFindex) && (geneTwo(ind)!= gTFindex))
        flag = FALSE;
        /*if ((po->IS_list)&&(!sublist[e->gene_one] || !sublist[e->gene_two])) TODO: Check sublist usage in file reading
        flag = FALSE;*/
    }

    if (!flag)
      continue;
    //Init Current block
    currBlock = new BicBlock();
    currBlock->score = min(2, (int)scores(ind));
    currBlock->pvalue = 1;
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
    int k=0;
    for(k = 0; k < components; k++) {
      if (gIsPValue)
      if ((pvalues[k] == currBlock->pvalue) &&(k >= 2) &&(vecScores[k]!=vecScores[k+1])) 
        break;
      if ((vecScores[k] == currBlock->score)&&(vecScores[k+1]!= currBlock->score)) 
        break;
    }
    components = k + 1;

    for (auto ki=0; ki < rowNumber; ki++) {
      candidates[ki] = TRUE;
    }

    for (auto ki=0; ki < components - 1 ; ki++) {
      candidates[vecGenes.at(ki)] = FALSE;
    }

    candidates[vecGenes.at(k)] = FALSE;
    swap(vecGenes[k], vecGenes[vecGenes.size()-1]);
    //------------------------------------------------------------------------------------------------------------
    bool *colcand = new bool[colNumber];
    for(auto ki = 0; ki < colNumber; ki++)
      colcand[ki] = FALSE;

    /* get init block */
    int threshold = floor(components * 0.7)-1;
    if(threshold <1)
      threshold=1;
    /*get the statistical results of each column produced by seed*/
    char *temptag = new char[colNumber];

    for(auto i=0;i<colNumber;i++) {
      colsStat[i] = 0;
      temptag[i] = 0;
    }
    Rcout << components << "Comp" << endl;
    for(auto i=1;i<components;i++) {
      getGenesFullLCS(&discreteInputData[vecGenes[0]*colNumber], &discreteInputData[vecGenes[i]*colNumber],temptag);
      for(auto j=0;j<colNumber;j++)
      {
        if(temptag[j]!=0)
          colsStat[j]++;
        temptag[j]=0;
      }
    }
    for(auto i=0;i<colNumber;i++) {
      if (colsStat[i] >= threshold) {
        colcand[i] = TRUE;
        cnt++;
      }
    }
    delete[] temptag;


    /* add some new possible genes */
    int m_ct=0;
    bool colChose = TRUE;
    for(auto ki=0;ki < rowNumber;ki++) {
      colChose=TRUE;
      if(!candidates[ki]) //if ((po->IS_list && !sublist[ki]) || !candidates[ki]) TODO Sublist check;
        continue;

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


    /* add genes that negative regulated to the consensus */
    char * reveTag;
    for (auto ki = 0; ki < rowNumber; ki++) {
      colChose=TRUE;
      if (!candidates[ki]) 
        continue;
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
      getGenesFullLCS(&discreteInputData[vecGenes[0]*colNumber],&discreteInputData[ki*colNumber],reveTag,lcsTags[vecGenes[1]]);
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
    /* save the current cluster*/
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

    if (currBlock->block_cols < 4 || components < 5) 
      continue;
    currBlock->block_rows = components;
    if (gIsPValue)
      currBlock->score = -(100*log(currBlock->pvalue));
    else
      currBlock->score = currBlock->block_rows * currBlock->block_cols;

    currBlock->genes.clear();
    for (auto ki=0; ki < components; ki++)
      currBlock->genes.push_back(vecGenes[ki]);
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


  qsort(arrBlocks, block_id, sizeof *arrBlocks, blockComp);

  int n = min(block_id, gRptBlock);
  bool flag;

  BicBlock **output = new BicBlock*[n]; // Array with filtered biclusters
  BicBlock **bb_ptr = output;
  BicBlock *b_ptr;

  double cur_rows, cur_cols;
  double inter_rows, inter_cols;

  /* the major post-processing here, filter overlapping blocks*/
  int i = 0, j = 0;
  while (i < block_id && j < n) {
    b_ptr = arrBlocks[i];
    cur_rows = b_ptr->block_rows;
    cur_cols = b_ptr->block_cols;

    flag = TRUE;
    int k = 0;
    while (k < j) {
      for(auto iter = output[k]->genes.begin(); iter != output[k]->genes.end(); iter++) {
        auto result1 = find(b_ptr->genes.begin(), b_ptr->genes.end(), *iter);
        if(result1==b_ptr->genes.end())
          inter_rows++;
      }
      for(auto iter = output[k]->conds.begin(); iter != output[k]->conds.end(); iter++) {
        auto result1 = find(b_ptr->conds.begin(), b_ptr->conds.end(), *iter);
        if(result1==b_ptr->conds.end())
          inter_cols++;
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
  // TODO RETURN OUTPUT BICLUSTERS - 
  // Return ouuput
  //-----------------------------------------------------------------------------------------------
  //Memory clearing
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


  return j;
}
// [[Rcpp::plugins(cpp11)]]

// Enable OpenMP (exclude macOS)
// [[Rcpp::plugins(openmp)]]


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

