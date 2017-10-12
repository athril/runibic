/***
Copyright (c) 2017 Patryk Orzechowski, Artur Pańszczyk

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
***/

#include <iostream>
#include <Rcpp.h>
#include <cstdlib>
#include <omp.h>
#include <vector>
#include <set>
#include <algorithm>
#include <utility>
#include <iterator>
#include "GlobalDefs.h"
#include  "fib.h"


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

bool is_higher(const triple* x, const triple* y) {
  if (x->lcslen > y->lcslen)
    return true;
  if (x->lcslen == y->lcslen)
    return (x->geneA<y->geneB);
  return false;
}

bool check_seed(int score, int geneOne, int geneTwo,  std::vector<BicBlock*> const &vecBlk, const int block_id, int rowNum) {
  vector<int> profiles(rowNum,0);
  int b1,b2,b3; // indexes for searching of first encounter
  b1 = b2 = -1;

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
void block_init(int score, int geneOne, int geneTwo, BicBlock *block, std::vector<int> &genes, std::vector<int> &scores, vector<bool> &candidates, const int cand_threshold, int *components, std::vector<long double> &pvalues, Params* params, std::vector<std::vector<int>> &lcsTags, std::vector<std::vector<int>> *inputData){
 
  int rowNum = gParameters.RowNumber;
  int colNum = gParameters.ColNumber;
  int cnt = 0, cnt_all=0, pid=0,row_all = rowNum;
  float cnt_ave=0;
  long double pvalue;

  int max_cnt, max_i;
  int t0,t1;

  
  t0=genes[0];
  t1=genes[1];
  /*************************calculate the lcs*********************************/
  //PO: It seems to be the same as calling 
  //PO: backtrackLCS(g1,g2)
  lcsTags.clear();
  lcsTags.resize(rowNum);
  
  lcsTags[t1] = getGenesFullLCS((*inputData)[t0],(*inputData)[t1]);
  set<int> colcand(lcsTags[t1].begin(), lcsTags[t1].end());
  std::vector<int> g1Common;
  //lcsLength[t1]=getGenesFullLCS(g1,g2,lcsTags[t1],NULL,colNum); 
  for (auto i = 0; i < (*inputData)[t0].size() ;i++){
    auto res = find(lcsTags[t1].begin(), lcsTags[t1].end(), (*inputData)[t0][i]);
    if(res!=lcsTags[t1].end())
      g1Common.push_back((*inputData)[t0][i]);
  }
 
  for(auto j=0;j<rowNum;j++) {
    std::vector<int> gJ;
    if (j==t1 || j==t0)
      continue;
    for(auto i=0;i<colNum;i++)
    {
      auto res2 = find(lcsTags[t1].begin(), lcsTags[t1].end(), (*inputData)[j][i]);
      if(res2!=lcsTags[t1].end())
        gJ.push_back((*inputData)[j][i]);
    }
    lcsTags[j] = getGenesFullLCS(g1Common,gJ);
    //lcsLength[j]= getGenesFullLCS(g1,(*inputData)[j].data(),lcsTags[j],lcsTags[t1],colNum); 
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

      cnt= count_if(lcsTags[i].begin(), lcsTags[i].end(), [&](int k) { return colcand.find(k) != colcand.end();});

      cnt_all += cnt;
      if (cnt < cand_threshold)
        candidates[i] = false;
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
      tempScore = (*components)*max_cnt;
    else
      tempScore = min(*components, max_cnt);
    if (tempScore > block->score)
      block->score = tempScore;
    if (pvalue < block->pvalue)
      block->pvalue = pvalue;
    genes.push_back(max_i);
    scores.push_back(tempScore);
    pvalues[pid++] = pvalue;

    for (auto it = colcand.begin(); it != colcand.end();){
      auto res = find(lcsTags[max_i].begin(), lcsTags[max_i].end(), *it);
      if(res==lcsTags[max_i].end())
        it = colcand.erase(it);
      else
        it++;
    }
    candidates[max_i] = FALSE;
  }
}

std::vector<int> getGenesFullLCS(std::vector<int> const &s1, std::vector<int> const &s2){

  vector<int> maxRecord;/*record the max value of matrix*/
  vector<int> lcsTag;
  int maxvalue,rank;
  
  int **C,**B;

  /*create matrix for lcs*/
  C = new int*[s1.size()+1];
  B = new int*[s1.size()+1];
  
  for(auto i=0;i<s1.size()+1;i++) {
    C[i] = new int[s2.size()+1];
    B[i] = new int[s2.size()+1];
  }
  
  /************initial the edge***************/
  for(auto i=0; i<s1.size()+1; i++) {
    C[i][0] = 0;
    B[i][0] = 0;
  }
  for(auto j=0; j<s2.size()+1; j++) {
    C[0][j] = 0;
    B[0][j] = 0; 
  }
  /************DP*****************************/
  for(auto i=1; i<s1.size()+1; i++) {
    for(auto j=1; j<s2.size()+1; j++) {
      if(s1[i-1] == s2[j-1]) {
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
  maxvalue = C[s1.size()][s1.size()];
  
  for (auto j=1;j<s2.size()+1;j++) {
    if (C[s1.size()][j] == maxvalue)
      maxRecord.push_back(j);
  }
  /*find all the columns of all LCSs*/
  if(maxRecord.size() > 0) {

    for (auto i=maxRecord.size()-1;i>=0;i--) {
      TrackBack(C,B, s1.size()+1,maxRecord[i]+1);      
      break;
    }
      for (auto i=1;i<s1.size()+1;i++) {
        for (auto j=1;j<s2.size()+1;j++) {
          if (C[i][j] == -1 && B[i][j]==1) {
              lcsTag.push_back(s1[i-1]);
          }
        } //end for j
      } // end for i  
     // Print the lcs
     //cout << "LCS of " << X << " and " << Y << " is " << lcs;
    }
  for(auto i=0;i<s1.size()+1;i++) {
    delete[] C[i];
    delete[] B[i];
  }
  delete[] C;
  delete[] B;

  return lcsTag;
}



/*track back the matrix*/
void TrackBack(int** pc,int** pb,int nrow,int ncolumn) {
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

void internalPairwiseLCS(std::vector<int> &x, std::vector<int> &y, std::vector<std::vector<int>> &c){

  for (auto i=0; i<x.size()+1; i++) {
    c[i].resize(y.size()+1);
    c[i][0]=0;
  }

  for (auto j=0; j<y.size()+1; j++) {
    c[0][j]=0;
  }
  for(auto i=1; i<x.size()+1; i++) {
    for(auto j=1; j<y.size()+1; j++) {
      if(x[i-1] == y[j-1]) {
        c[i][j] = c[i-1][j-1] + 1;
      }
      else {
        c[i][j] = std::max(c[i][j-1],c[i-1][j]);
      }
    }
  }
}
void internalCalulateLCS(std::vector<std::vector<int>> &inputMatrix, std::vector<triple> &out, bool useFib){

  int PART = 4;
  int step = inputMatrix.size()/PART;
  int size = (PART-1)*(step*(step-1)/2);
  int rest = step+(inputMatrix.size()%PART);
  size+= rest*(rest-1)/2;
  triple** triplets = new triple*[size];
  struct fibheap *heap;
  heap = fh_makeheap();
  fh_setcmp(heap, edge_cmpr);
 
// TODO: change into parallel version
// there should be 1-level for loop across all combinations of pairs of rows
//  #pragma omp parallel for private(a,b,i,j,res) schedule(dynamic)
//  for ( auto k=0; k<size; k++ ) {
//    auto i = k/discreteInput.nrow(); auto j=k%discreteInput.nrow(); 
  //triple __cur_min = {0, 0, po->COL_WIDTH};
  triple __cur_min = {0, 0, gParameters.ColWidth};
  triple *_cur_min = &__cur_min;
  triple **cur_min = & _cur_min;
  int k=0;
  for(auto p = 0; p < PART; p++){
    auto endi = (p+1)*step;
    if(p == PART-1)
      endi = inputMatrix.size();
    for (auto i=p*step; i<endi; i++) {
      for (auto j=i+1; j<endi; j++) {
        triplets[k] = new triple;
        triplets[k]->geneA = i;
        triplets[k]->geneB = j;
        k++;
      }
    }
  }
#pragma omp parallel for shared(triplets) schedule(dynamic)
  for(auto p = 0; p < k; p++){
    vector<int> a = inputMatrix[triplets[p]->geneA];
    vector<int> b = inputMatrix[triplets[p]->geneB];
    vector< vector<int> > res(a.size()+1);
    internalPairwiseLCS(a,b,res);
    triplets[p]->lcslen= res[res.size()-1][res.size()-1];
  }

  for(auto p = 0; p < k; p++){
    fh_insert(heap, (void *)triplets[p]);
    if(useFib){
      if (size < HEAP_SIZE) 
      {
        fh_insert(heap, (void *)triplets[p]);
      }
      else
      {
        if (edge_cmpr(cur_min, triplets[p]) < 0)
        {
          /* Remove least value and renew */
          fh_extractmin(heap);
          fh_insert(heap, (void *)triplets[p]);
          /* Keep a memory of the current min */
          *cur_min = (triple *)fh_min(heap);
        }
      }
    }
  }

  if(useFib){
    for(int i=size-1; i>=0; i--){
      triple *res= static_cast<triple *>(fh_extractmin(heap));
      out.push_back(*res);
    }
    reverse(out.begin(), out.end());
  }
  else  {
    sort( triplets, triplets+size, &is_higher);
    for (int i=0; i<size; i++) {
      out.push_back(*(triplets[i]));
    }
  }
  for(int i=0;i <size; i++)
    delete triplets[i];
  delete[] triplets;
  free(heap);
}

bool blockComp(BicBlock*lhs, BicBlock* rhs) {
/* compare function for qsort, descending by score */ 
  return lhs->score > rhs->score;
}

double calculateQuantile(Rcpp::NumericVector vecData, int size, double qParam)
{
  double delta = (size-1)*qParam;
  if(delta < 0)
    delta = 0;
  int i = floor(delta);
  delta=delta-i;
  if(i < size - 1)
    return (1-delta)*vecData(i) + (delta)*vecData(i+1);
  else 
    return -1;
}
