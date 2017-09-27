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

bool check_seed(int score, int geneOne, int geneTwo,  BicBlock** vecBlk, const int block_id, int rowNum);

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
Rcpp::NumericMatrix unisort(Rcpp::NumericMatrix x){
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
//' cluster(c(13,12,11,7,5,3),c(0,1,2,0,0,1), c(3,2,3,2,1,3),3,4)
//'
//' @export
// [[Rcpp::export]]
int cluster(Rcpp::NumericVector scores, Rcpp::NumericVector geneOne, Rcpp::NumericVector geneTwo, int rowNumber, int colNumber)
{
  int block_id = 0;
  int cnt = 0;
  
  BicBlock** arrBlocks = new BicBlock*[gSchBlock];
  for(auto ind =0; ind<gSchBlock; ind++)
    arrBlocks[ind] = NULL;
  BicBlock *b; // bicluster candidate
  
	vector<int> *vecGenes, *vecScores, *vecBicGenes, *vecAllInCluster; // helpful vectors/stacks
  vecGenes = new vector<int>();
	vecScores = new vector<int>();
  vecAllInCluster = new vector<int>();
  
	int components;
  
  //Memory allocation
	int *colsStat = new int[colNumber]; // column statisctics array    
	long double *pvalues = new long double[rowNumber]; // pvalues array

	bool *candidates = new bool[rowNumber];
  short *lcsLength = new short[rowNumber];
  char** lcsTags = new char*[rowNumber];

  for(auto ind = 0; ind < rowNumber; ind++)
  {
    lcsTags[ind] = new char[colNumber];
  }

  //Main loop
  for(auto ind = 0; ind < scores.size(); ind++)
  {
    /* check if both genes already enumerated in previous blocks */
		bool flag = TRUE;
		/* speed up the program if the rows bigger than 200 */
	  if (rowNumber > 250)
		{ 
      auto result1 = find(vecAllInCluster->begin(), vecAllInCluster->end(), geneOne(ind));
      auto result2 = find(vecAllInCluster->begin(), vecAllInCluster->end(), geneTwo(ind));

			if ( result1 != vecAllInCluster->end() && result2 != vecAllInCluster->end())
				flag = FALSE;
			else if (gIsTFname && (geneOne(ind) !=  gTFindex) && (geneTwo(ind)!= gTFindex))
				flag = FALSE;
			/*else if (gIsList&&(!sublist[geneOne(ind)] || !sublist[geneTwo(ind)]))    TODO: Check sublist usage in file reading
				flag =FALSE;*/
		}
		else   
		{
			flag = check_seed(scores(ind),geneOne(ind), geneTwo(ind), arrBlocks, block_id, rowNumber);
			if (gIsTFname && (geneOne(ind) !=  gTFindex) && (geneTwo(ind)!= gTFindex))
      flag = FALSE;
			/*if ((po->IS_list)&&(!sublist[e->gene_one] || !sublist[e->gene_two])) TODO: Check sublist usage in file reading
				flag = FALSE;*/
    }
    
		if (!flag) continue;
  }
  
  //Memory clearing
  for(auto ind =0; ind<gSchBlock; ind++)
    if(arrBlocks!=NULL)
      delete arrBlocks[ind];
  delete[] arrBlocks;
  delete[] colsStat;
  delete[] lcsLength;
  delete[] candidates;
  delete[] pvalues;
  delete vecGenes;
  delete vecScores;
  delete vecAllInCluster;
  for(auto ind = 0; ind < rowNumber; ind++)
  {
    delete[] lcsTags[ind];
  }
  delete[] lcsTags; 
  return 0;
}
bool check_seed(int score, int geneOne, int geneTwo,  BicBlock** vecBlk, const int block_id, int rowNum)
{
  int profiles[rowNum];
	int b1,b2,b3; // indexes for searching of first 
  b1 = b2 = -1;

  std::fill(profiles, profiles+rowNum,0);	
  
  for (auto ind = 0; ind < block_id; ind++)
  {
    auto result1 = find(vecBlk[ind]->genes->begin(), vecBlk[ind]->genes->end(), geneOne);
    auto result2 = find(vecBlk[ind]->genes->begin(), vecBlk[ind]->genes->end(), geneTwo);
    if ( result1 != vecBlk[ind]->genes->end()  && result2 != vecBlk[ind]->genes->end() ) 
      return FALSE;
    if (result1 != vecBlk[ind]->genes->end() && b1 == -1)
    {
      b1 = ind;
    }
    if(result2 != vecBlk[ind]->genes->end() && b2 == -1)
    {
      b2 = ind;
    }
  }

	if ( (b1 == -1)||(b2 == -1) ) 
    return TRUE;
	else
	{
		for (auto i = 0; i < vecBlk[b1]->block_rows; i++)
			profiles[vecBlk[b1]->genes->at(i)]++;
		for (auto i = 0; i < vecBlk[b2]->block_rows; i++)
			profiles[vecBlk[b2]->genes->at(i)]++;
		for (auto i = 0; i < rowNum; i++)
 			if (profiles[i] > 1) 
				return FALSE;
		b3 = max(vecBlk[b1]->block_cols, vecBlk[b2]->block_cols);
		if ( score < b3) 
			return FALSE;
		else 
			return TRUE;
	}
  return FALSE;
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

