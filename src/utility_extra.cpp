#include <Rcpp.h>
#include <math.h>
#include "geometry.h"
#include "see.h"
#include "utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector destinationAngle_rcpp(
    double a,
    NumericMatrix p1,
    NumericMatrix P1
  ) {
  
  NumericVector a2 = angle2_rcpp(p1, P1);
  
  NumericVector dest_angle = headingAngle_rcpp(a2, a);
  
  return(dest_angle);
}


NumericMatrix omit_rows(NumericMatrix mat, IntegerVector omit) {
  int n_rows = mat.nrow();
  int n_cols = mat.ncol();
  int l = omit.length();
  
  NumericMatrix new_mat(n_rows - l, n_cols);
  CharacterVector old_row_names = rownames(mat);
  CharacterVector new_row_names(n_rows - l);
  
  int j = 0;
  
  for(int i = 0; i < n_rows; i++) {
    if (is_false(any(omit == i))) {
      new_mat(j, _) = mat(i, _);
      new_row_names[j] = old_row_names[i];
      j += 1;
    }
  }
  
  rownames(new_mat) = new_row_names;
  // colnames(new_mat) = colnames(mat);
  
  return new_mat;
}


// [[Rcpp::export]]
Nullable<NumericMatrix> predClose_rcpp(
    int n,
    NumericMatrix p1,
    double a1,
    NumericMatrix p2,
    NumericVector r,
    NumericMatrix centres,
    NumericMatrix p_pred,
    List objects
  ) {
  int p_pred_n_rows = p_pred.nrow();
  
  if (p_pred_n_rows == 1) {
    return(R_NilValue);
  }
  
  NumericMatrix ps = omit_rows(p2, n);
  
  int p2_n_rows = p2.nrow();
  
  IntegerVector occluded = seq_len(p2_n_rows) - 1;
  
  occluded.erase(n);
  
  LogicalVector is_seen = seesMany_rcpp(p1, ps, objects);
  
  occluded = occluded[!is_seen];
  
  occluded.push_front(n);
  
  p_pred = omit_rows(p_pred, occluded);
  
  if (p_pred.nrow() <= 1) {
    return(R_NilValue);
  }
  
  p2 = omit_rows(p2, occluded);
  
  LogicalVector in_front = (minAngle_rcpp(a1, angle2_rcpp(p1, p2)) < 85.0) &
    (minAngle_rcpp(a1, angle2_rcpp(p1, p_pred)) < 85.0);
  
  if (is_false(any(in_front))) {
    return(R_NilValue);
  }
  
  IntegerVector idx_in_front = seq_len(in_front.length()) - 1;
  NumericMatrix p_pred_in_front = omit_rows(p_pred, idx_in_front[!in_front]);
  
  CharacterVector names_in_front = rownames(p_pred_in_front);
  NumericVector r_in_front = r[names_in_front];
  
  NumericMatrix dist(r_in_front.length(), centres.nrow());
  
  for(int i = 0; i < centres.nrow(); i++) {
     dist(_, i) = dist1_rcpp(centres(i, _), p_pred_in_front);
  }
  
  for(int i = 0; i < r_in_front.length(); i++) {
    dist(i, _) = dist(i, _) - (r[n] + r_in_front[i]);
  } 
  
  return dist;
}


// [[Rcpp::export]]
List eObjects_rcpp(NumericMatrix p1, NumericMatrix p2, NumericVector r) {
  NumericVector d = dist1_rcpp(p1, p2);
  
  NumericVector a12 = angle2_rcpp(p1, p2);
  
  int p2_n_rows = p2.nrow();
  
  NumericVector r_rep = rep_len(r, p2_n_rows);
  
  NumericVector theta = atan(r_rep / d) * 180.0 / M_PI;
  
  NumericVector angles_ac(p2_n_rows);
  NumericVector angles_cw(p2_n_rows);
  
  for(int i = 0; i < p2_n_rows; i++) {
    angles_ac[i] = fmod(360.0 + a12[i] + theta[i], 360.0);
    angles_cw[i] = fmod(360.0 + a12[i] - theta[i], 360.0);
  }
  
  NumericMatrix ac = aTOd_rcpp(angles_ac);
  NumericMatrix cw = aTOd_rcpp(angles_cw);
  
  for(int i = 0; i < p2_n_rows; i++) {
    ac(i, _) = ac(i, _) * d[i];
    cw(i, _) = cw(i, _) * d[i];
  }
  
  NumericMatrix ac_t = transpose(ac);
  NumericMatrix cw_t = transpose(cw);
  
  for(int j = 0; j < p2.ncol(); j++) {
    ac_t(j, _) = ac_t(j, _) + p1[j];
    cw_t(j, _) = cw_t(j, _) + p1[j];
  }
  
  NumericMatrix ac_t_t = transpose(ac_t);
  NumericMatrix cw_t_t = transpose(cw_t);
  
  List ac_cw = List::create(Named("ac") = ac_t_t, Named("cw") = cw_t_t);
  
  return(ac_cw);
}


IntegerVector fix(IntegerVector x) {
  if (is_true(all(is_na(x)))) {
    IntegerVector na_res = IntegerVector::create(NA_REAL);
    return(na_res);
  }
  LogicalVector false_true = LogicalVector::create(false, true);
  LogicalVector true_false = LogicalVector::create(true, false);
  if (is_true(all(is_na(x) == false_true))) {
    x = IntegerVector::create(x[0], 11);
  } else if (is_true(all(is_na(x) == true_false))) {
    x = IntegerVector::create(1, x[1]); 
  } else if (x.length() > 1) {
    x = seq(x[0], x[1]); 
  }
  return(x);
}


List check_list_na(List x) {
  LogicalVector valid(x.length());
  for(int i = 0; i < x.length(); i++) {
    if (is_false(any(is_na(x)))) {
      valid[i] = true;
    }
  }
  
  return(x[valid]);
}


List check_list_null(List x) {
  LogicalVector valid(x.length());
  for(int i = 0; i < x.length(); i++) {
    if (x[i] != R_NilValue) {
      valid[i] = true;
    }
  }
  
  return(x[valid]);
}


List check_list_length(List x) {
  LogicalVector valid(x.length());
  for(int i = 0; i < x.length(); i++) {
    NumericVector x_i = x[i];
    if (x_i.length() > 0) {
      valid[i] = true;
    } ;
  }
  
  return(x[valid]);
}


// [[Rcpp::export]]
Nullable<NumericVector> iCones_rcpp(
    NumericMatrix p1,
    double a,
    NumericMatrix p2,
    NumericVector r,
    List objects
) {
  if (p2.nrow() == 0) {
    return(R_NilValue);
  }
  
  List ends = eObjects_rcpp(p1, p2, r);
  NumericMatrix ends_ac = ends[0];
  NumericMatrix ends_cw = ends[1];
  
  IntegerMatrix endCones(p2.nrow(), ends.length());
  rownames(endCones) = rownames(p2);
  
  for(int i = 0; i < ends.length(); i++) {
    endCones(_, i) = Iangle_rcpp(p1, a, ends[i]);
  }
  
  List cList(p2.nrow());
  cList.names() = rownames(p2);
  
  for(int i = 0; i < cList.length(); i++) {
    cList[i] = fix(endCones(i, _));
  }
  
  cList = check_list_na(cList);
  
  if (cList.length() == 0) {
    return(R_NilValue);
  }
  
  cList = check_list_null(cList);
  
  if (cList.length() == 0) {
    return(R_NilValue);
  }
  
  NumericMatrix coneLineEnds = c_vd_rcpp(
    seq_len(11) - 1, p1, rep(1.0, 11), a, get_vels(), get_angles(), 0.5
  );
  
  List cDist(cList.length());
  
  for(int i = 0; i < cList.length(); i++) {
    NumericVector cList_i = cList[i];
    if (cList_i.length() == 1) {
      bool sees_goal = seesGoal_rcpp(p1, p2(i, _), objects);
      if (!sees_goal) {
        cList[i] = NumericVector::create();
      } else {
        cDist[i] = dist1_rcpp(p1, p2(Range(i, i), _));
      }
    } else {
      for(int j = 0; j < cList_i.length(); j++) {
        NumericVector P_n_vec = line_line_intersection_rcpp(
          p1, coneLineEnds(cList_i[j], _), ends_ac(i, _), ends_cw(i, _)
        );
        
        bool sees_goal = seesGoal_rcpp(p1, P_n_vec, objects);
        if (!sees_goal) {
          cList_i[j] = NA_REAL;
        } else {
          P_n_vec.attr("dim") = Dimension(1, 2);
          NumericMatrix P_n = as<NumericMatrix>(P_n_vec);
          if (j == 0) {
            cDist[i] = dist1_rcpp(p1, P_n);
          } else {
            NumericVector d = dist1_rcpp(p1, P_n);
            for(int k = 0; k < d.length(); k++) {
              NumericVector cDist_i = cDist[i];
              cDist_i.push_back(d[k]);
              cDist[i] = cDist_i;
            }
          }
        }
      }
      cList[i] = cList_i[!is_na(cList_i)];
    }
  }
  
  cList = check_list_length(cList);
  
  if (cList.length() == 0) {
    return(R_NilValue);
  }
  
  cDist = check_list_null(cDist);
  
  for(int i = 0; i < cList.length(); i++) {
    NumericVector cDist_i = cDist[i];
    cDist_i.names() = int2char(cList[i]);
  }
  
  IntegerVector outCones = IntegerVector::create();
  NumericVector out = NumericVector::create();
  
  for(int i = 0; i < 11; i++) {
    NumericVector d;
    
    for(int j = 0; j < cDist.length(); j++) {
      NumericVector cDist_j = cDist[j];
      IntegerVector cDist_j_names_int = char2int(cDist_j.names());
      
      for(int k = 0; k < cDist_j.length(); k++) {
        if (cDist_j_names_int[k] == i) {
          d.push_back(cDist_j[k]);
        }
      }
    }
    
    if (d.length() > 0) {
      outCones.push_back(i);
      out.push_back(min(d));
    }
  }
  
  out.names() = int2char(outCones);
  
  return(out);
}


// [[Rcpp::export]]
NumericVector iCones2Cells_rcpp(NumericVector iC, double v, double tStep = 0.5) {
  NumericVector vels = NumericVector::create(1.5, 1.0, 0.5);
  
  NumericVector out = rep(iC, 3) - rep_each(v * tStep * vels, iC.length());
  
  IntegerVector rings = IntegerVector::create(0, 11, 22);
  out.names() = int2char(rep(char2int(iC.names()), 3) + 
    rep_each(rings, iC.length()));
  
  return(out);
}


// [[Rcpp::export]]
Nullable<NumericVector> blockedAngle_rcpp(int n, List state, NumericMatrix p_pred, List objects) {
  NumericMatrix p_mat = state["p"];
  NumericVector p1_vec = p_mat(n, _);
  p1_vec.attr("dim") = Dimension(1, 2);
  NumericMatrix p1 = as<NumericMatrix>(p1_vec);
  NumericVector a = state["a"];
  NumericMatrix p2 = omit_rows(p_pred, n);
  NumericVector r = state["r"];
  
  Nullable<NumericVector> _iC = iCones_rcpp(p1, a[n], p2, r, objects);
  if (_iC == R_NilValue) {
    return(R_NilValue);
  }
  NumericVector iC = as<NumericVector>(_iC);
  NumericVector v = state["v"];
  NumericVector cells = iCones2Cells_rcpp(iC, v[n]);
  
  return(cells);
}