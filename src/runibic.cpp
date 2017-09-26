#include <iostream>
#include <Rcpp.h>
#include <cstdlib>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <utility>

using namespace std;
using namespace Rcpp;




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


// [[Rcpp::export]]
Rcpp::NumericMatrix calculateLCS(Rcpp::NumericVector m, Rcpp::NumericVector n) {
  NumericMatrix pc(m.size(), n.size());
  NumericMatrix pb(m.size(), n.size());

//  #pragma omp parallel for shared(pc,pb)
  for (auto j=0; j<m.size(); j++) {
    pc(j,0)=0;
    pb(j,0)=0;
  }

//  #pragma omp parallel for shared(pc,pb)
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

