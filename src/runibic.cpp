#include <iostream>
#include <Rcpp.h>
#include <cstdlib>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <utility>

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

