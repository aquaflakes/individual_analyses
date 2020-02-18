// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <iostream>
// #include "common_op.h"
using namespace Rcpp;
using namespace std;



// [[Rcpp::export]]
NumericMatrix matRevComp(NumericMatrix seqs)
  {
    NumericMatrix rc(seqs.nrow(),seqs.ncol());  //NumericMatrix(seqs.nrow(),seqs.ncol());
    for (size_t i=0; i<seqs.size(); i++) rc[i]= seqs[i];
    std::reverse(rc.begin(),rc.end());
    // std::reverse(seqs.begin(),seqs.end());
    return rc;
  }
