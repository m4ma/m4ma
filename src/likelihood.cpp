#include <Rcpp.h>
#include "utility.h"
#include "pCNL.h"

using namespace Rcpp;


NumericVector get_mum(NumericVector p) {
  CharacterVector nest_names = {"Central", "NonCentral", "acc", "const", "dec"};
  NumericVector nest_params = p[nest_names];
  NumericVector mum = 1.0 / (1.0 - nest_params);
  
  return mum;
}


//' Subject Log-likelihood
//' 
//' Calculate the log-likelihood of an observation for a single pedestrian.
//'
//' @param subject List with subject data.
//' @param p Numeric vector with subject parameters.
//' @param n Integer scalar subject index.
//' @param nests List of vectors with utility indices.
//' @param alpha List of vectors with alpha values.
//' @param cell_nest Numeric matrix with nest indices for each cell.
//' @param min_like Numeric scalar minimum likelihood return value.
//' 
//' @returns Numeric scalar subject log likelihood.
//' @export
// [[Rcpp::export]]
double like_subject(List subject, NumericVector p, int n, List nests, List alpha, NumericMatrix cell_nest, double min_like = 1e-10) {

  NumericVector u = utility(p, n, subject["v"], subject["d"], subject["BA"],
                            subject["GA"], subject["ID"], subject["FL"], subject["WB"],
                            subject["ok"], subject["group"]);
  
  int cell = subject["cell"];

  NumericVector mum = get_mum(p);

  double lprob = log(pcnl_rcpp(cell_nest(cell, _), u, mum, nests, alpha, 1.0));

  return std::max(lprob, log(min_like));
}


//' State Log-likelihood
//' 
//' Calculate the log-likelihood of observations for a state as the sum of pedestrian log-likelihoods.
//'
//' @param state List of lists with subject data.
//' @param p Numeric matrix with subject parameters for each subject.
//' @param nests List of vectors with utility indices.
//' @param alpha List of vectors with alpha values.
//' @param cell_nest Numeric matrix with nest indices for each cell.
//' @param min_like Numeric scalar minimum likelihood return value.
//' 
//' @returns Numeric vector of state subject likelihoods.
//' @export
// [[Rcpp::export]]
NumericVector like_state(List state, NumericMatrix p, List nests, List alpha, NumericMatrix cell_nest, double min_like = 1e-10) {
  int k = state.length();
  NumericVector llike_state (k);
  
  for(int i = 0;  i < k; ++i) {
    List state_i = state[i];
    NumericVector p_i = p(i, _);
    p_i.names() = colnames(p);
    llike_state[i] = like_subject(state_i, p_i, i, nests, alpha, cell_nest, min_like);
  }
  
  return llike_state;
}


//' Trace Log-likelihood
//' 
//' Calculate the log-likelihood of a trace of states as the sum over states and subject log-likelihoods.
//'
//' @param p Numeric matrix with subject parameters for each subject.
//' @param trace_rcpp List of lists of lists with subject data.
//' @param nests List of vectors with utility indices.
//' @param alpha List of vectors with alpha values.
//' @param cell_nest Numeric matrix with nest indices for each cell.
//' @param min_like Numeric scalar minimum likelihood return value.
//' @param mult Numeric scalar likelihood sum multiplier.
//' 
//' @returns Numeric scalar trace log-likelihood.
//' @export
// [[Rcpp::export]]
double msumlogLike(NumericMatrix p, List trace_rcpp, List nests, List alpha, NumericMatrix cell_nest, double min_like = 1e-10, double mult = -1.0) {
  int l = trace_rcpp.length();
  double llike_trace = 0.0;
  
  for(int i = 0;  i < l; ++i) {
    List state_i = trace_rcpp[i];
    llike_trace += sum(like_state(state_i, p, nests, alpha, cell_nest, min_like));
  }
  
  return llike_trace * mult;
}