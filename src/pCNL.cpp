#include <Rcpp.h>
using namespace Rcpp;


// Probability of nest given utility and parameters alpha, mu, mum
double pnest(NumericVector u_i, NumericVector alpha_i, 
             double mum, double mu_mum) {
  
  // probability mass function
  double p_nest = pow(sum(alpha_i * exp(mum * u_i)), mu_mum);
  
  return p_nest;
}


// Probability of alternative in nests given utility and parameters alpha, mum
NumericVector painnest(NumericVector u_i, NumericVector alpha_i, double mum) {
  
  // Probability mass function
  NumericVector p_a_nest = alpha_i * exp(mum * u_i);
  NumericVector p_a_nest_reg = p_a_nest;
  LogicalVector bad = is_nan(p_a_nest) | is_infinite(p_a_nest);
  
  // If any NaN or Inf set to 1 and other to 0
  if(is_true(any(bad))) {
    p_a_nest_reg = ifelse(bad, 1.0, 0.0); 
  }
  
  // Avoid zero division if all 0
  if(is_false(all(p_a_nest_reg == 0.0))) {
    p_a_nest_reg = p_a_nest_reg / sum(p_a_nest_reg);
  }
  
  return p_a_nest_reg;
}


// Check nest likelihood for inf and zero division
NumericVector check_pnest(NumericVector p_nest) {
  
  NumericVector p_nest_reg = p_nest;
  LogicalVector bad = is_nan(p_nest) | is_infinite(p_nest);
  
  // If any NaN or Inf set to 1 and other to 0
  if(is_true(any(bad))) {
    p_nest_reg = ifelse(bad, 1.0, 0.0);
  }
  
  // Avoid zero division if all 0
  if(is_false(all(p_nest_reg == 0.0))) {
    p_nest_reg = p_nest_reg / sum(p_nest_reg);
  }
  
  return p_nest_reg;
}


//' Probability of the Conditional Nested Logit Model
//'
//' @param cell Alternative vector with nest indices.
//' @param utility Vector with utility for each alternative.
//' @param mum Vector with nests associations ranging between 0 and 1.
//' @param nests List of vectors with utility indices.
//' @param alpha List of vectors with alpha values.
//' @param mu General nest association ranging between 0 and 1.
//'
//' @return Probability of alternative \code{cell} given \code{utility}, \code{alpha}, and parameters \code{mum} and \code{mu}.
//' @export
// [[Rcpp::export]]
double pcnl_rcpp(NumericVector cell, NumericVector utility, 
                 NumericVector mum, List nests, List alpha, 
                 double mu) {
  
  double max_utility = max(utility);
  int m = nests.length();
  NumericVector mu_mum = mu / mum;
  NumericVector p_nest (m);
  List p_a_nest (m);
  
  // Loop over nests
  for(int i = 0; i < m; ++i) {
    IntegerVector nest_i = nests[i];
    NumericVector utility_i = utility[nest_i]; // get nest u
    NumericVector utility_i_reg = utility_i - max_utility; // avoid numerical instabilities
    p_nest[i] = pnest(utility_i_reg, alpha[i], mum[i], mu_mum[i]); // nest prob
    p_a_nest[i] = painnest(utility_i_reg, alpha[i], mum[i]); // nest alternative probs
  }
  
  p_nest = check_pnest(p_nest);
  
  NumericVector p_a_nest_dir = p_a_nest[cell[0]]; // get probs for dir nests
  NumericVector p_a_nest_ang = p_a_nest[cell[1]]; // get probs for ang nests
  
  // Final probability of cell in nests
  double p = p_a_nest_dir[cell[2]] * p_nest[cell[0]] + 
    p_a_nest_ang[cell[3]] * p_nest[cell[1]];
  
  return p;
}

// [[Rcpp::export]]
double pmnl_rcpp(int cell, NumericVector utility, LogicalMatrix ok) {
  
  // copy ok elements of utility into u
  int nok = ok.length()+1;
  NumericVector u(sum(ok)+1);
  int ok_count = 1;
  u(0) = utility(0);
  for(int i = 1; i < nok; ++i) {
    int ii = i - 1;
    if (ok[ii]) {
      u[ok_count] = utility[i];
      ok_count = ok_count + 1;
    }
  }
  
  // Exponential protected against overflows
  double maxu = max(u);
  NumericVector u0 = u - maxu;
  u0 = exp(u0);
  // Assumes cell is in ok
  double uc = utility(cell)-maxu;
  uc = exp(uc);
  double p = uc/sum(u0);
  
  return(p);

} 

