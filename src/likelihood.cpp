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
double like_state(List state, NumericMatrix p, int n, List nests, List alpha, NumericMatrix cell_nest) {
  NumericVector p_n = p(n, _);
  p_n.names() = colnames(p);

  NumericVector u = utility(p_n, n, state["v"], state["d"], state["BA"],
                            state["GA"], state["ID"], state["FL"], state["WB"],
                            state["ok"], state["group"]);
  
  Rcout << u;
  
  int cell = state["cell"];

  NumericVector mum = get_mum(p_n);

  double prob = pcnl_rcpp(cell_nest(cell + 1, _), u, mum, nests, alpha, 1.0);

  return prob;
}
