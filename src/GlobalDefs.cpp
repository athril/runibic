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
extern Params gParameters;

int edge_cmpr(void *a, void *b)
{
	float score_a, score_b;
	score_a = ((triple *)a)->lcslen;
	score_b = ((triple *)b)->lcslen;

	if (score_a < score_b) return -1;
	if (score_a == score_b) return 0;
	return 1;
}

bool check_seed(int score, int geneOne, int geneTwo,  BicBlock** vecBlk, const int block_id, int rowNum) {
  int profiles[rowNum];
  int b1,b2,b3; // indexes for searching of first encounter
  b1 = b2 = -1;

  std::fill(profiles, profiles+rowNum,0);

  for (auto ind = 0; ind < block_id; ind++) {
    auto result1 = find(vecBlk[ind]->genes.begin(), vecBlk[ind]->genes.end(), geneOne);
    auto result2 = find(vecBlk[ind]->genes.begin(), vecBlk[ind]->genes.end(), geneTwo);
    if ( result1 != vecBlk[ind]->genes.end()  && result2 != vecBlk[ind]->genes.end() ){
      return FALSE;
    }
      
    if (result1 != vecBlk[ind]->genes.end() && b1 == -1) {
      b1 = ind;
    }
    if(result2 != vecBlk[ind]->genes.end() && b2 == -1) {
      b2 = ind;
    }
  }
  if ( (b1 == -1)||(b2 == -1) )
    return TRUE;
  else {
    for (auto i = 0; i < vecBlk[b1]->block_rows; i++)
      profiles[vecBlk[b1]->genes.at(i)]++;
    for (auto i = 0; i < vecBlk[b2]->block_rows; i++)
      profiles[vecBlk[b2]->genes.at(i)]++;
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


//lcsTags is vector<vector<int>>
void block_init(int score, int geneOne, int geneTwo, BicBlock *block, vector<int> *genes, vector<int> *scores, bool *candidates, const int cand_threshold, int *components, vector<int> *allincluster, long double *pvalues, Params* params,short *lcsLength, char** lcsTags, vector<int> *inputData) {
  int rowNum = gParameters.RowNumber;
  int colNum = gParameters.ColNumber;
  int cnt = 0, cnt_all=0, pid=0,row_all = rowNum;
  float cnt_ave=0;
  long double pvalue;
  int max_cnt, max_i;
  int t0,t1;
  int *arrRows, *arrRowsB;

  arrRows = new int[rowNum];
  arrRowsB = new int[rowNum];

  bool *colcand = new bool[colNum];
  for (auto ind=0; ind< colNum; ind++)
    colcand[ind] = false;
  int *g1, *g2;
  t0=genes->at(0);
  t1=genes->at(1);

  //PO: no need for this code?
  //PO: No need to initialize lcsTags. The vector is empty. lcsLength is 0.
  /*cut-from*/

  g1 = &((*inputData)[t0*colNum]);
  g2 = &((*inputData)[t1*colNum]);
  for(auto i=0;i<rowNum;i++) {
    lcsLength[i]=0;
    for(auto j=0;j<colNum;j++)
      lcsTags[i][j]=0;
  }
  /*cut-to*/



  /*************************calculate the lcs*********************************/


  //PO: It seems to be the same as calling 
  //PO: backtrackLCS(g1,g2)
  lcsLength[t1]=getGenesFullLCS(g1,g2,lcsTags[t1],NULL,colNum); 
  for(auto i=0;i<colNum;i++) {
    if(lcsTags[t1][i]!=0) {
      colcand[i]=TRUE;
    }
  }

  for(auto j=0;j<rowNum;j++) {
    if (j==t1 || j==t0)
      continue;


  //PO: This should be modified to:
  //PO: lcsLength=backtrackLCS(g1, ??? ).size()
    lcsLength[j]= getGenesFullLCS(g1,&(*inputData)[j*colNum],lcsTags[j],lcsTags[t1],colNum); 
  }



  while (*components < rowNum) {
    max_cnt = -1;
    max_i = -1;
    (*components)++;
    cnt_all =0;
    cnt_ave = 0;
    /******************************************************/
    /*add a function of controling the bicluster by pvalue*/
    /******************************************************/
    for (auto i=0; i< rowNum; i++) {
      if (!candidates[i]) {
        continue;
      }
      //if (gIsList && !sublist[i]) continue;// TODO Check sublist
      cnt = 0;
      for (auto j=0; j< colNum; j++) {
        if (colcand[j] && lcsTags[i][j]!=0){
          cnt++;
        }
          
      }
      cnt_all += cnt;
      if (cnt < cand_threshold)
        candidates[i] = FALSE;
      if (cnt > max_cnt) {
        max_cnt = cnt;
        max_i = i;
      }
    }
    cnt_ave = cnt_all/row_all;

    long double one = 1;
    long double poisson=one/exp(cnt_ave);
    for (auto i=0;i<max_cnt+300;i++) {
      if (i>(max_cnt-1)) 
        pvalue=pvalue+poisson;
      else 
        poisson=poisson*cnt_ave/(i+1);
    }
    if (gParameters.IsCond) {
      if (max_cnt < gParameters.ColWidth || max_i < 0|| max_cnt < block->cond_low_bound) break;
    }
    else {
      if (max_cnt < gParameters.ColWidth || max_i < 0){
        break;        
      }
    }
    int tempScore = 0;
    if (gParameters.IsArea)
      tempScore = *components*max_cnt;
    else
      tempScore = min(*components, max_cnt);
    if (tempScore > block->score)
      block->score = tempScore;
    if (pvalue < block->pvalue)
      block->pvalue = pvalue;
    genes->push_back(max_i);
    scores->push_back(tempScore);
    pvalues[pid++] = pvalue;

    for (auto i=0; i< colNum; i++)
      if (colcand[i] &&(lcsTags[max_i][i]==0))
        colcand[i] = FALSE;
    candidates[max_i] = FALSE;
  }
  delete[] colcand;
  delete[] arrRows;
  delete[] arrRowsB;
}






//int Rcpp::IntegerMatrix pairwiseLCS(Rcpp::IntegerVector x, Rcpp::IntegerVector y) <- computes LCS between 2 pairs
//Rcpp::IntegerVector backtrackLCS(Rcpp::IntegerMatrix c, Rcpp::IntegerVector x, Rcpp::IntegerVector y)  <- gets exact LCS of 2 pairs: x,y
//
// Calling:
// {lcs_tg,maxvalue} = getGenesFullLCS(s1,s2,_,_._,reverse)
// should be the same as calling: 
// backtrackLCS(s1,s2).

int getGenesFullLCS(const int *s1, const int *s2,char *lcs_tg ,char *lcs_seed,  int colNum,bool reverse) {
  vector<int> maxRecord;/*record the max value of matrix*/
  int maxvalue,rank,length1,length2;
  int *temp1,*temp2;
  short **C,**B;
  temp1 = new int[colNum];
  temp2 = new int[colNum];

  maxvalue = 0;
  length1=length2=0;
  rank = gParameters.Divided;

  for(auto i=0;i<colNum;i++) {
    temp2[i]=0;
    temp1[i]=0;
  }

  /*get the sorted sequence*/
  for(auto i=1;i<=rank;i++) {
    for(auto j=0;j<colNum;j++) {
      if(lcs_seed!=NULL) {
        if(lcs_seed[j] != 0) {
          if(s1[j] == i)
            temp1[length1++]=j+1;
          /*printf("[j:%d,s2[j]:%d] ",j,s2[j]);*/
          if(reverse){
            if(s2[j] == rank-i+1)
						  temp2[length2++]=j+1;
          }
          else{
            if(s2[j] == i)
            temp2[length2++]=j+1;
          }
          
        }
      }
      else {
        if(s1[j] == i)
          temp1[length1++]=j+1;
        /*printf("[j:%d,s2[j]:%d] ",j,s2[j]);*/
        if(s2[j] == i)
          temp2[length2++]=j+1;
      }
    }
  }
  if(gParameters.Quantile < 0.5) {
    for(auto i=rank*(-1);i<=-1;i++) {
      for(auto j=0;j<colNum;j++) {
        if(lcs_seed!=NULL) {
          if(lcs_seed[j] != 0) {
            if(s1[j] == i)
              temp1[length1++]=j+1;
            /*printf("[j:%d,s2[j]:%d] ",j,s2[j]);*/
            if(s2[j] == i)
              temp2[length2++]=j+1;
          }
        }
        else {
        if(s1[j] == i)
          temp1[length1++]=j+1;
        /*printf("[j:%d,s2[j]:%d] ",j,s2[j]);*/
        if(s2[j] == i)
          temp2[length2++]=j+1;
        }
      }
    }
  }
  /*create matrix for lcs*/
  C = new short*[length1+1];
  B = new short*[length1+1];
  for(auto i=0;i<length1+1;i++) {
    C[i] = new short[length2+1];
    B[i] = new short[length2+1];
  }
  
  /************initial the edge***************/
  for(auto i=0; i<length1+1; i++) {
    C[i][0] = 0;
    B[i][0] = 0;
  }
  for(auto j=0; j<length2+1; j++) {
    C[0][j] = 0;
    B[0][j] = 0;
  }
  /************DP*****************************/
  for(auto i=1; i<length1+1; i++) {
    for(auto j=1; j<length2+1; j++) {
      if(temp1[i-1] == temp2[j-1]) {
        C[i][j] = C[i-1][j-1] + 1;
        B[i][j] = 1;
      }
      else if(C[i-1][j] >= C[i][j-1]) {
        C[i][j] = C[i-1][j];
        B[i][j] = 2;
      }
      else {
        C[i][j] = C[i][j-1];
        B[i][j] = 3;
      }
    }
  }
  maxvalue = C[length1][length2];
  for (auto j=1;j<length2+1;j++) {
    if (C[length1][j] == maxvalue)
      maxRecord.push_back(j);
  }
  /*find all the columns of all LCSs*/
  if(maxRecord.size() > 0) {
    for (auto i=maxRecord.size()-1;i>=0;i--) {
      TrackBack(C,B, length1+1,maxRecord[i]+1);      
      break;
    }
    if(lcs_tg!=NULL)
    {
      for (auto i=1;i<length1+1;i++) {
        for (auto j=1;j<length2+1;j++) {
          if (C[i][j] == -1 && B[i][j]==1) {
              lcs_tg[temp1[i-1]-1] = 1;
          }
        } //end for j
      } // end for i  
    }
  }
  
  delete[] temp1;
  delete[] temp2;
  for(auto i=0;i<length1+1;i++) {
    delete[] C[i];
    delete[] B[i];
  }
  delete[] C;
  delete[] B;

  return maxvalue;
}



/*track back the matrix*/
void TrackBack(short** pc,short** pb,int nrow,int ncolumn) {
  int ntemp;
  if(nrow == 0 || ncolumn == 0)
    return;
  ntemp = pb[nrow-1][ncolumn-1];
  pc[nrow-1][ncolumn-1] = -1;
         
  switch(ntemp) {
    case 1:
      TrackBack(pc, pb, nrow-1, ncolumn-1);
      break;
    case 2:
      TrackBack(pc, pb, nrow-1, ncolumn);
      break;
    case 3:
      TrackBack(pc, pb, nrow, ncolumn-1);
      break;
    default:
      break;
  }
}


short* getRowData(int index) {
  return NULL;
}


bool blockComp(BicBlock*lhs, BicBlock* rhs) {
/* compare function for qsort, descending by score */ 
  return lhs->score > rhs->score;
}

double calculateQuantile(Rcpp::NumericVector vecData, int size, double qParam)
{
  double delta = (size-1)*qParam;
  int i = floor(delta);
  delta=delta-i;
  if(i < size - 1)
    return (1-delta)*vecData(i) + (delta)*vecData(i+1);
  else 
    return -1;
}