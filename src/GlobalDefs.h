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
  , Divided(100)
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

double calculateQuantile(Rcpp::NumericVector vecData, int size, double qParam);
bool check_seed(int score, int geneOne, int geneTwo,  BicBlock** vecBlk, const int block_id, int rowNum);
void block_init(int score, int geneOne, int geneTwo, BicBlock *block, std::vector<int> *genes, std::vector<int> *scores, bool *candidates, const int cand_threshold, int *components, std::vector<int> *allincluster, long double *pvalues, Params* params, short *lcsLength, char** lcsTags, std::vector<int> *inputData);
int getGenesFullLCS(const int *s1, const int *s2,char *lcs_tg = NULL,char *lcs_seed = NULL, int colNum = 0, bool reverse = false);
void TrackBack(short** pc,short** pb,int nrow,int ncolumn);
short* getRowData(int index);
bool blockComp(BicBlock* lhs, BicBlock* rhs);
Rcpp::List fromBlocks(BicBlock ** blocks, const int numBlocks, const int nr, const int nc);

#endif
