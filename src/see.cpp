#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector line_line_intersection_rcpp(
    NumericVector P1,
    NumericVector P2,
    NumericVector P3,
    NumericVector P4,
    bool interior_only=false) {
  
  double dx1 = P1[0] - P2[0];
  double dx2 = P3[0] - P4[0];
  double dy1 = P1[1] - P2[1];
  double dy2 = P3[1] - P4[1];
  
  double D = dx1 * dy2 - dy1 * dx2;
  
  if (D == 0.0) {
    NumericVector inf_res = {R_PosInf, R_PosInf};
    return(inf_res);
  }
  
  double D1 = P1[0] * P2[1] - P1[1] * P2[0];
  double D2 = P3[0] * P4[1] - P3[1] * P4[0];
  
  double X = (D1 * dx2 - dx1 * D2) / D;
  double Y = (D1 * dy2 - dy1 * D2) / D;
  
  if (interior_only) {
    double lambda1 = -((X - P1[0]) * dx1 + (Y - P1[1]) * dy1)/(pow(dx1, 2) + pow(dy1, 2));
    double lambda2 = -((X - P3[0]) * dx2 + (Y - P3[1]) * dy2)/(pow(dx2, 2) + pow(dy2, 2));
    
    if (!((lambda1 > 0.0) & (lambda1 < 1.0) & (lambda2 > 0.0) & (lambda2 < 1.0))) {
      NumericVector na_res = {NA_REAL, NA_REAL};
      return(na_res);
    }
  }
  
  NumericVector res = NumericVector::create(X, Y);
  
  return(res);
}


// [[Rcpp::export]]
bool seesGoal_rcpp(NumericVector p_n, NumericVector P_n, List objects) {
  int l = objects.length();
  
  NumericVector idx_P3_x = NumericVector::create(0, 0, 0, 1);
  NumericVector idx_P3_y = NumericVector::create(0, 0, 1, 0);
  NumericVector idx_P4_x = NumericVector::create(1, 0, 1, 1);
  NumericVector idx_P4_y = NumericVector::create(0, 1, 1, 1); 
  
  for (int i = 0; i < l; ++i) {
    
    List object_i = objects[i];
    NumericVector x = object_i["x"];
    NumericVector y = object_i["y"];
    
    for (int j = 0; j < 4; ++j) {
      NumericVector P3 = NumericVector::create(x[idx_P3_x[j]], y[idx_P3_y[j]]);
      NumericVector P4 = NumericVector::create(x[idx_P4_x[j]], y[idx_P4_y[j]]);
      
      NumericVector intersects = line_line_intersection_rcpp(
        p_n, P_n, P3, P4, true
      );
      
      if (is_true(any(is_finite(intersects)))) {
        return(false);
      }
    }
  }
  
  return(true);
}


NumericVector get_state_P_n(int n, List state, int offset=0) {
  List state_P = state["P"];
  
  NumericMatrix state_P_n = state_P[n - 1];
  
  int i = state_P_n.attr("i");
  
  NumericVector P_n = state_P_n(i + offset - 1, _);
  IntegerVector idx = IntegerVector::create(0, 1);
  
  return(P_n[idx]);
}


// [[Rcpp::export]]
bool seesCurrentGoal_rcpp(int n, List state, List objects, int offset=0) {
  NumericMatrix state_p = state["p"];
  NumericVector p_n = state_p(n - 1, _);
  NumericVector P_n = get_state_P_n(n, state, offset);
  
  bool sees_goal = seesGoal_rcpp(p_n, P_n, objects);
  
  return(sees_goal);
}


// [[Rcpp::export]]
LogicalVector seesMany_rcpp(NumericVector p1, NumericMatrix ps, List objects) {
  int n_rows = ps.nrow();
  
  LogicalVector sees_many(n_rows);
  
  for (int i = 0; i < n_rows; ++i) {
    sees_many[i] = seesGoal_rcpp(ps(i, _), p1, objects);
  }
  
  return(sees_many);
}


// [[Rcpp::export]]
LogicalVector seesGoalOK_rcpp(int n, List objects, List state, NumericMatrix centres, LogicalVector ok) {
  if (is_true(any(ok))) {
    NumericVector P_n = get_state_P_n(n, state);
    
    for (int i = 0; i < ok.length(); ++i) {
      if (ok[i]) {
        ok[i] = seesGoal_rcpp(centres(i, _), P_n, objects);
      }
    }
  }
  
  return(ok);
}