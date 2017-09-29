#ifndef GLOBALDEFS_H
#define GLOBALDEFS_H

#include <vector>
#include <algorithm>
// Default parameters for program
/* TODO: Add program options init*/
const int gSchBlock = 200; // TODO: check full name of this option; Blocks to post process
const bool gIsTFname = false; //TODO: check full name and usage of this option;
const bool gIsList = false; //TODO: check full name and usage of this option;
const double gTolerance = 0.85;// TODO: check usage of this option
const int gDataMode = 0;// TODO: check usage of this option
const int gQuantile = 0.5;;// TODO: check usage of this option
const bool gIsCond = false;// TODO: check usage of this option
const bool gIsArea = false; // TODO: check usage of this option
const bool gIsPValue = false;
const int gRptBlock = 100; // blocks to output
const double gFilter = 0.25;


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


bool check_seed(int score, int geneOne, int geneTwo,  BicBlock** vecBlk, const int block_id, int rowNum);
void block_init(int score, int geneOne, int geneTwo, BicBlock *block, std::vector<int> *genes, std::vector<int> *scores, bool *candidates, const int cand_threshold, int *components, std::vector<int> *allincluster, long double *pvalues, int rowNum, int colNum,short *lcsLength, char** lcsTags, std::vector<int> *inputData);
int getGenesFullLCS(const int *s1, const int *s2,char *lcs_tg = NULL,char *lcs_seed = NULL, int colNum = 0, bool reverse = false);
void TrackBack(short** pc,short** pb,int nrow,int ncolumn);
short* getRowData(int index);
bool blockComp(BicBlock* lhs, BicBlock* rhs);

#endif
