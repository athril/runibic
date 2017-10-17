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


#ifndef GLOBALDEFS_H
#define GLOBALDEFS_H

#include <vector>
#include <algorithm>
#include <Rcpp.h>


class Params{
public:
  Params()
  : IsDiscrete(false)
  , SchBlock(200)
  , Tolerance(0.85)
  , Quantile(0.5)
  , IsCond(false)
  , IsArea(false)
  , RptBlock(100)
  , Filter(1)
  , Shuffle(0)
  , Divided(0)
  , RowNumber(0)
  , ColNumber(0){};
  int RowNumber;
  int ColNumber;
  bool IsDiscrete;
  int SchBlock;
  double Tolerance;
  double Quantile;
  bool IsCond; //the flag using the lower bound of condition number (5 persents of the gene number)
  bool IsArea; // the flag using area as the value of bicluster to determine when stop
  bool IsPValue; // the flag to enlarge current biclsuter by the pvalue constrain
  int RptBlock;
  double Filter;
  int Shuffle;
  int ColWidth;
  int Divided;
  void InitOptions(int rowNum, int colNum){
    RowNumber = rowNum;
    ColNumber = colNum;
    ColWidth = std::max(3+floor(colNum/30),4.0);
    if(Divided==0) {
      if(rowNum > 2000) {
        Divided = 15;
      }
      else {
        Divided = colNum;
      }
    }
    
  }
};
/* biclustering block */
typedef struct BicBlock {
  std::vector<int> genes;
  std::vector<int> conds;
  int score;
  int block_rows;
  int block_cols;
  int block_rows_pre;
  int cond_low_bound;
  double significance;
  long double pvalue;
} BicBlock;

struct triple {
  int geneA;
  int geneB;
  int lcslen;
};
static const int HEAP_SIZE = 20000000;

int edge_cmpr(void *a, void *b);
double calculateQuantile(Rcpp::NumericVector vecData, int size, double qParam);
bool check_seed(int score, int geneOne, int geneTwo,  std::vector<BicBlock*> const &vecBlk, const int block_id, int rowNum);
void block_init(int score, int geneOne, int geneTwo, BicBlock *block, std::vector<int> &genes, std::vector<int> &scores, std::vector<bool> &candidates, const int cand_threshold, int *components, std::vector<long double> &pvalues, Params* params, std::vector<std::vector<int>> &lcsTags, std::vector<std::vector<int>> *inputData);
std::vector<int> getGenesFullLCS(std::vector<int> const &s1, std::vector<int> const &s2);
void TrackBack(int** pc,int** pb,int nrow,int ncolumn);
short* getRowData(int index);
bool blockComp(BicBlock* lhs, BicBlock* rhs);
void internalPairwiseLCS(std::vector<int> &x, std::vector<int> &y, std::vector<std::vector<int> > &c);
void internalCalulateLCS(std::vector<std::vector<int>> &inputMatrix, std::vector<triple> &out, bool useFib);
Rcpp::List fromBlocks(BicBlock ** blocks, const int numBlocks, const int nr, const int nc);
Rcpp::IntegerVector backtrackLCS(Rcpp::IntegerVector x, Rcpp::IntegerVector y);
#endif

