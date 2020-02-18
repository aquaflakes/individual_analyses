#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix place(NumericMatrix seq, NumericVector pos,NumericVector win)
{
  //get ini pos
  int win_half_size= (win.length()-1)/2;
  int halfLen= (seq.ncol()-1)/2;
  NumericVector pos_start= pos-win_half_size+ halfLen;
  NumericVector pos_end= pos+win_half_size + halfLen;


  for(int i=0; i< seq.nrow(); i++)
  {
    if (NumericVector::is_na(pos[i])) continue;
    for(int p=pos_start[i]; p<= pos_end[i]; p++)
    {
      if (seq(i,p)>win(p-pos_start[i])) seq(i,p)=win(p-pos_start[i]); //take highest block factor
      // seq(i,p)*=win(p-pos_start[i]);
    }
  }
  return(seq);
}


// [[Rcpp::export]]
DataFrame frag_mid_and_pos(LogicalMatrix MNaseCuts)
{

  int ncol_= MNaseCuts.ncol(); int half_size= ((ncol_)-1)/2;
  int nrow_=MNaseCuts.nrow();
  vector<int> midPoints;
  vector<int> lengths;

  for(int i=0; i< nrow_; i++)
  {
    // bool old_fragment=false;
    int fragLen=0;
    for(int p=0; p< ncol_; p++)
    {
      if (MNaseCuts(i,p)==false)
      {
        fragLen+=1;
        if (p==ncol_-1) //count fragment if the final pos is false
        {
          lengths.push_back(fragLen);
          midPoints.push_back( (((p+1)-1)+((p+1)-1-(fragLen-1))) /2 - half_size);
        }
      }
      if (MNaseCuts(i,p)==true)
      {
        if(fragLen!=0)
        {
          lengths.push_back(fragLen);
          midPoints.push_back( ((p-1)+(p-1-(fragLen-1))) /2 - half_size);
          fragLen=0;
        }
      }
    }
  }
  return(Rcpp::DataFrame::create(_["fragMid"]=midPoints,_["fragLen"]=lengths));
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
cutseq <-function(seq,freeDNA_cut_prob,cut_times)
{
  matsize=seqNum*(seqLenHalf*2+1)
  cutPosTotal=matrix(0,seqNum,seqLenHalf*2+1);

  for (i in 1:cut_times)
  {
    cutPos= matrix(sample(c(0L,1L),matsize,replace = T,prob = c(1-freeDNA_cut_prob,freeDNA_cut_prob)),seqNum,seqLenHalf*2+1)
    cutRnd= matrix(runif(matsize),seqNum,seqLenHalf*2+1)
    cutRnd= cutPos*cutRnd
    cutPos= ifelse(cutRnd<=seq & cutRnd!=0,TRUE,FALSE)
    cutPosTotal= cutPosTotal | cutPos
  }
  return(cutPosTotal)
}
*/
