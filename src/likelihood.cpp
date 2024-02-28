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


//' Observation Log-likelihood
//' 
//' Calculate the log-likelihood of an observation for a single subject or iteration.
//'
//' @param obs List with observation data.
//' @param p Numeric vector with subject parameters.
//' @param n Integer scalar subject index.
//' @param nests List of vectors with utility indices.
//' @param alpha List of vectors with alpha values.
//' @param cell_nest Numeric matrix with nest indices for each cell.
//' @param min_like Numeric scalar minimum likelihood return value.
//' 
//' @returns Numeric scalar observation log-likelihood.
//' @export
// [[Rcpp::export]]
double like_observation(
    List obs,
    NumericVector p,
    int n, List nests,
    List alpha,
    NumericMatrix cell_nest,
    double min_like = 1e-10
) {

  NumericVector u = utility(p, obs["v"], obs["d"], obs["BA"],
                            obs["GA"], obs["ID"],
                            obs["inGroup"], obs["IDPreUtility"],
                            obs["FL"], obs["WB"],
                            obs["okID"]);
  
  int cell = obs["cell"];
  LogicalMatrix ok = obs["ok"];

//  NumericVector mum = get_mum(p);
//  double lprob = log(pcnl_rcpp(cell_nest(cell, _), u, mum, nests, alpha, 1.0));

  double lprob = log(pmnl_rcpp(cell, u, ok));
  
  return std::max(lprob, log(min_like));
}


//' State Log-likelihood
//' 
//' Calculate the log-likelihood of observations for a state as the sum of observation log-likelihoods.
//'
//' @param state List of lists with subject data.
//' @param p Numeric matrix with subject parameters for each subject.
//' @param nests List of vectors with utility indices.
//' @param alpha List of vectors with alpha values.
//' @param cell_nest Numeric matrix with nest indices for each cell.
//' @param elements Character string indicating whether to return a trace with 
//' iterations or subjects as elements. Can be either \code{"iterations"} or 
//' \code{"subjects"}.
//' @param min_like Numeric scalar minimum likelihood return value.
//' 
//' @returns Numeric vector of state observaton log-likelihoods.
//' @export
// [[Rcpp::export]]
NumericVector like_state(
    List state,
    int ti,
    NumericMatrix p,
    List nests,
    List alpha,
    NumericMatrix cell_nest,
    std::string elements = "iterations",
    double min_like = 1e-10
) {
  int k = state.length();
  NumericVector llike_state (k);
  
  for(int i = 0;  i < k; ++i) {
    List state_i = state[i];
    
    int n;
    NumericVector p_n;
    int p_i;
    
    if (elements == "iterations") {
      n = i;
      NumericVector pn = state_i["pn"];
      p_i = pn(i) - 1; // c++ indexing
      p_n = p(p_i, _);
    } else {
      n = state_i["n"];
      n = n - 1; // c++ indexing
      p_n = p(ti, _);
    }

    p_n.names() = colnames(p);
    llike_state[i] = like_observation(state_i, p_n, n, nests, alpha, cell_nest, min_like); 
  }
  
  return llike_state;
}


//' Trace Log-likelihood
//' 
//' Calculate the log-likelihood of a trace of states as the sum over states and observation log-likelihoods.
//'
//' @param p Numeric matrix with subject parameters for each subject.
//' @param trace_rcpp List of lists of state lists with observation data.
//' @param nests List of vectors with utility indices.
//' @param alpha List of vectors with alpha values.
//' @param cell_nest Numeric matrix with nest indices for each cell.
//' @param min_like Numeric scalar minimum likelihood return value.
//' @param mult Numeric scalar likelihood sum multiplier.
//' 
//' @returns Numeric scalar trace log-likelihood.
//' @export
// [[Rcpp::export]]
double msumlogLike(
    NumericMatrix p,
    List trace_rcpp,
    List nests,
    List alpha,
    NumericMatrix cell_nest,
    double min_like = 1e-10,
    double mult = -1.0
) {
  std::string elements = trace_rcpp.attr("elements");
  int l = trace_rcpp.length();
  double llike_trace = 0.0;
  
  for(int i = 0;  i < l; ++i) {
    List state_i = trace_rcpp[i];
    llike_trace += sum(like_state(state_i, i, p, nests, alpha, cell_nest, elements, min_like));
  }
  
  return llike_trace * mult;
}