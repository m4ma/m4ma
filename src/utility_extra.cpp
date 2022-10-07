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

// [[Rcpp::export]]
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


NumericMatrix get_p1(int n, NumericMatrix p_mat) {
  NumericVector p1_vec = p_mat(n, _);
  p1_vec.attr("dim") = Dimension(1, 2);
  NumericMatrix p1 = as<NumericMatrix>(p1_vec);
  CharacterVector p_mat_rownames = rownames(p_mat);
  rownames(p1) = as<CharacterVector>(p_mat_rownames[n]);
  
  return(p1);
}


// [[Rcpp::export]]
Nullable<NumericVector> blockedAngle_rcpp(int n, List state, NumericMatrix p_pred, List objects) {
  NumericMatrix p_mat = state["p"];
  NumericMatrix p1 = get_p1(n, p_mat);
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


// [[Rcpp::export]]
Nullable<List> getLeaders_rcpp(
  int n,
  List state,
  NumericMatrix centres,
  List objects,
  bool onlyGroup = false,
  bool preferGroup = true,
  bool pickBest = false
) {
  NumericMatrix p_mat = state["p"];
  NumericMatrix p1 = get_p1(n, p_mat);
  NumericVector a = state["a"];
  NumericVector v = state["v"];
  double a1 = a[n];
  double v1 = v[n];
  IntegerVector n_vec = IntegerVector::create(n);
  NumericMatrix ps = omit_rows(p_mat, n_vec);
  
  IntegerVector occluded = seq_along(v) - 1;
  occluded.erase(n);
  LogicalVector is_seen = seesMany_rcpp(p1, ps, objects);
  occluded = occluded[!is_seen];
  occluded.push_front(n);
  
  CharacterVector a_names = a.names();
  NumericVector a2;
  CharacterVector a2_names;
  for(int i = 0; i < a.length(); i++) {
    if(is_false(any(i == occluded))) {
      a2.push_back(a[i]);
      a2_names.push_back(a_names[i]);
    }
  }
  a2.names() = a2_names;
  
  NumericMatrix p2 = omit_rows(p_mat, occluded);
  if (p2.nrow() == 0) {
    return(R_NilValue);
  }
  
  NumericMatrix P1;
  
  if (is<List>(state["P"])) {
    List P1_list = state["P"];
    NumericMatrix P1_mat = P1_list[n];
    int attr_i = P1_mat.attr("i");
    P1 = P1_mat(Range(attr_i - 1, attr_i - 1), Range(0, 1));
  } else {
    NumericMatrix P1_mat = state["P"];
    P1 = P1_mat(Range(n, n), Range(0, 1));
  }
  
  NumericVector I_n = Iangle_rcpp(p1, a1, p2);
  I_n = I_n[!is_na(I_n)];
  if (I_n.length() == 0) {
    return(R_NilValue);
  }
  
  CharacterVector I_n_names = I_n.names();
  CharacterVector p2_rownames = rownames(p2);
  
  NumericMatrix p2_ring(I_n_names.length(), p2.ncol());
  
  for(int i = 0; i < I_n_names.length(); i++) {
    for(int j = 0; j < p2_rownames.length(); j++) {
      if (I_n_names[i] == p2_rownames[j]) {
        p2_ring(i, _) = p2(j, _);
      }
    }
  }
  
  NumericVector d_ring = dist1_rcpp(p1, p2_ring);
  
  NumericVector bins_ring = NumericVector::create(0, 0.5, 1, 5) * v1 * 0.5;
  
  NumericVector ring = 4 - bin_vector(d_ring, bins_ring);
  ring.names() = I_n_names;
  ring = ring[!is_na(ring)];
  if (ring.length() == 0) {
    return(R_NilValue);
  }
  
  CharacterVector ring_names = ring.names();
  NumericVector candidates = I_n[ring_names];
  candidates = candidates + 11 * (ring - 1);
  
  CharacterVector candidates_names = candidates.names();
  NumericVector group = state["group"];
  NumericVector group_candidates = group[candidates_names];
  LogicalVector inGroup = group_candidates == group[n];
  inGroup.names() = group_candidates.names();
  
  if (onlyGroup) {
    if (is_false(any(inGroup))) {
      return(R_NilValue);
    } else {
      candidates = candidates[inGroup];
    } 
  } else if (preferGroup && is_true(any(inGroup))) {
    candidates = candidates[inGroup];
  }
  
  NumericVector a2_candidates = a2[candidates_names];
  NumericVector d_n = Dn_rcpp(p1, P1);
  NumericVector angles(a2_candidates.length());
  
  for(int i = 0; i < a2_candidates.length(); i++) {
    NumericVector a2_canidates_min_angle = minAngle_rcpp(a2_candidates[i], d_n);
    angles[i] = a2_canidates_min_angle[0];
  }
  angles.names() = candidates_names;
  
  LogicalVector ok = angles < 90.0;
  if (is_false(any(ok))) {
    return(R_NilValue);
  }
  
  candidates = candidates[ok];
  angles = angles[ok];
  
  NumericVector leaders;
  CharacterVector leaders_names;
  if(is_false(any(duplicated(candidates)))) {
    leaders = candidates;
    leaders_names = leaders.names();
  } else {
    leaders = unique(candidates);
    leaders_names = leaders.names();
    for(int i = 0; i < leaders.length(); i++) {
      double leaders_i = leaders[i];
      NumericVector leaders_angles_i = angles[candidates == leaders_i];
      double min_angle_i = which_min(leaders_angles_i);
      CharacterVector candidates_names_i = candidates_names[candidates == leaders_i];
      leaders_names[i] = candidates_names_i[min_angle_i];
    }
    angles = angles[leaders_names];
  }
  
  NumericMatrix d(leaders.length(), 33);
  rownames(d) = leaders_names;
  for(int i = 0; i < leaders.length(); i++) {
    double leaders_i = leaders[i];
    NumericVector centres_i = centres(leaders_i - 1, _);
    d(i, _) = dist1_rcpp(centres_i, centres);
  }
  
  IntegerVector best;
  if (pickBest) {
    best = which_min(angles);
  } else {
    best = seq_along(angles) - 1;
  }
  
  NumericVector leaders_best = leaders[best];
  CharacterVector leaders_names_best = leaders_best.names();
  
  NumericMatrix dists(best.length(), d.ncol());
  for(int i = 0; i < best.length(); i++) {
    dists(i, _) = d(best[i], _);
  }
  rownames(dists) = leaders_names_best;
  colnames(dists) = int2char(seq_len(33));
  
  NumericMatrix leaders_out(3, best.length());
  leaders_out(0, _) = leaders_best;
  NumericVector leaders_angles = angles[leaders_names_best];
  leaders_out(1, _) = leaders_angles / 90.0;
  LogicalVector leaders_inGroup = inGroup[leaders_names_best];
  leaders_out(2, _) = as<NumericVector>(leaders_inGroup);
  rownames(leaders_out) = CharacterVector::create(
    "cell", "angleDisagree", "inGroup"
  );
  colnames(leaders_out) = leaders_names_best;
  
  List out_list = List::create(
    Named("dists") = dists,
    Named("leaders") = leaders_out
  );
  
  return(out_list);
}
  
  