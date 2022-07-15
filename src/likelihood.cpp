#include <Rcpp.h>
#include "m4ma.h"
using namespace Rcpp;


NumericVector get_mum(NumericVector p) {
  CharacterVector nest_names = {"Central", "NonCentral", "acc", "const", "dec"};
  NumericVector nest_params = p[nest_names];
  NumericVector mum = 1.0 / (1.0 - nest_params);
  
  return mum;
}


// [[Rcpp::export]]
double like_state(List state, NumericMatrix p, int n, List nests, List alpha, NumericMatrix cell_nest, double min_like = 1e-10) {
  NumericVector p_n = p(n, _);
  p_n.names() = colnames(p);

  NumericVector u = utility(p_n, n, state["v"], state["d"], state["BA"],
                            state["GA"], state["ID"], state["FL"], state["WB"],
                            state["ok"], state["group"]);
  
  int cell = state["cell"];

  NumericVector mum = get_mum(p_n);

  double lprob = log(pcnl_rcpp(cell_nest(cell, _), u, mum, nests, alpha, 1.0));

  return std::max(lprob, log(min_like));
}


// [[Rcpp::export]]
NumericVector like_states(List states, NumericMatrix p, List nests, List alpha, NumericMatrix cell_nest, double min_like = 1e-10) {
  int k = states.length();
  NumericVector llike_states (k);
  
  for(int i = 0;  i < k; ++i) {
    List state_i = states[i];
    llike_states[i] = like_state(state_i, p, i, nests, alpha, cell_nest, min_like);
  }
  
  return llike_states;
}


// [[Rcpp::export]]
double msumlogLike(List trace, NumericMatrix p, List nests, List alpha, NumericMatrix cell_nest, double min_like = 1e-10, double mult = -1.0) {
  int l = trace.length();
  double llike_trace = 0.0;
  
  for(int i = 0;  i < l; ++i) {
    List states_i = trace[i];
    llike_trace += sum(like_states(states_i, p, nests, alpha, cell_nest, min_like));
  }
  
  return llike_trace * mult;
}