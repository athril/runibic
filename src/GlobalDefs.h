#ifndef GLOBALDEFS_H
#define GLOBALDEFS_H
#include <vector>

// Default parameters for program
/* TODO: Add program options init*/
const int gSchBlock = 200; // TODO: check full name of this option;
const bool gIsTFname = false; //TODO: check full name and usage of this option;
int gTFindex = -1; // Index EOF?
bool gIsList = false; //TODO: check full name and usage of this option;
/* biclustering block */
typedef struct BicBlock
{
	std::vector<int> *genes;
	std::vector<int> *conds;
	int score;
	int block_rows;
	int block_cols;
	int block_rows_pre;
	int cond_low_bound;
	double significance;
	long double pvalue;
} BicBlock;

#endif
