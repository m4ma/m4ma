#include <Rcpp.h>
#include <math.h>
#include "geometry.h"
#include "see.h"
#include "utils.h"
using namespace Rcpp;

//' Angle to Destination
//'
//' Compute the angles between point `p1` and destinations `P1`.
//' 
//' The matrices `p1` and `p2` must have the same shape.
//' 
//' @param a Numeric scalar heading angle of `p1`.
//' @param p1 Numeric matrix with shape 1x2 (x and y).
//' @param P1 Numeric matrix with shape Nx2 (x and y).
//' 
//' @return Named numeric vector of length 11.
// [[Rcpp::export]]
NumericVector destinationAngle_rcpp(
    double a,
    NumericMatrix p1,
    NumericMatrix P1
  ) {
  
  // angle between p1 and goal P1
  NumericVector a2 = angle2_rcpp(p1, P1);
  
  // angular difference for 11 angles
  NumericVector dest_angle = headingAngle_rcpp(a2, a);
  
  return(dest_angle);
}


//' Predicted Distance to Close Front
//'
//' Compute the distance from edge of body of pedestrian `n` to other pedestrians
//' in front that can be seen.
//' 
//' Returns `NULL` if no other pedestrians are present, if they cannot be seen,
//' or if they are not in front of pedestrian `n`.
//' 
//' @param n Integer scalar index of current pedestrian in pedestrian matrix.
//' @param p1 Numeric matrix with shape 1x2 (x and y) indicating the position
//' of pedestrian `n`.
//' @param a1 Numeric scalar heading angle of pedestrian `n`.
//' @param p2 Numeric matrix with shape Nx2 (x and y) indicating the positions
//' of all pedestrians.
//' @param r Numeric vector with radii of pedestrian bodies.
//' @param centres Numeric matrix with x and y for each cell centre (33x2).
//' @param p_pred Numeric matrix with shape Nx2 (x and y) as predicted positions
//' of all pedestrians.
//' @param objects List containing a list for each object. An object has
//' two length-two numeric vectors of x- and y-coordinates.
//' 
//' @return Numeric matrix with rows for each other pedestrian and columns for
//' each cell (Nx33) or NULL.
//' 
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
  n = n - 1; // c++ indexing
  
  int p_pred_n_rows = p_pred.nrow();
  
  if (p_pred_n_rows == 1) {
    return(R_NilValue);
  }
  
  // remove row n from p2
  NumericMatrix ps = omit_rows(p2, IntegerVector::create(n));
  
  int p2_n_rows = p2.nrow();
  
  // init index vector of occluded pedestrians
  IntegerVector occluded = seq_len(p2_n_rows) - 1;
  
  occluded.erase(n); // remove self from occluded
  
  // remove seen peds from occluded
  LogicalVector is_seen = seesMany_rcpp(p1, ps, objects);
  
  occluded = occluded[!is_seen];
  
  // add self to start of occluded
  occluded.push_front(n);
  
  // remove occluded rows from predicted pos
  p_pred = omit_rows(p_pred, occluded);
  
  if (p_pred.nrow() <= 1) { // if all other are occluded (1 is self)
    return(R_NilValue);
  }
  
  // remove occluded rows from actual pos
  p2 = omit_rows(p2, occluded);
  
  // check which peds are in front
  LogicalVector in_front = (minAngle_rcpp(a1, angle2_rcpp(p1, p2)) < 85.0) &
    (minAngle_rcpp(a1, angle2_rcpp(p1, p_pred)) < 85.0);
  
  if (is_false(any(in_front))) { // if none are in front
    return(R_NilValue);
  }
  
  // remove not in front from predicted pos
  IntegerVector idx_in_front = seq_len(in_front.length()) - 1;
  NumericMatrix p_pred_in_front = omit_rows(p_pred, idx_in_front[!in_front]);
  
  CharacterVector names_in_front = rownames(p_pred_in_front);
  NumericVector r_in_front = r[names_in_front];
  
  NumericMatrix dist(r_in_front.length(), centres.nrow());
  rownames(dist) = names_in_front;
  
  // calc dist from each cell centre to predicted pos
  for(int i = 0; i < centres.nrow(); i++) {
     dist(_, i) = dist1_rcpp(centres(i, _), p_pred_in_front);
  }
  
  // add radii (object size) to dist
  for(int i = 0; i < r_in_front.length(); i++) {
    dist(i, _) = dist(i, _) - (r[n] + r_in_front[i]);
  } 
  
  return dist;
}


//' Egocentric Objects
//'
//' Compute the end points of line segments that are orthogonal to a line
//' between point `p1` and points `p2`. The width of the segments are equal to
//' the width `r` of `p2` from the perspective of `p1`.
//' 
//' @param p1 Numeric matrix with shape 1x2 (x and y).
//' @param p2 Numeric matrix with shape Nx2 (x and y).
//' @param r Numeric vector of length N.
//' 
//' @return List:
//'   \describe{
//'     \item{ac}{Numeric matrix with rows as \emph{anticlockwise} end points to each row in `p2`
//'       and columns as x- and y-coordinates.}
//'     \item{cw}{Numeric matrix with rows as \emph{clockwise} end points to each row in `p2`
//'       and columns as x- and y-coordinates.}
//'   }
//'   
// [[Rcpp::export]]
List eObjects_rcpp(NumericMatrix p1, NumericMatrix p2, NumericVector r) {
  // calc dist between p1 and p2
  NumericVector d = dist1_rcpp(p1, p2);
  
  // calc angle between p1 and p2
  NumericVector a12 = angle2_rcpp(p1, p2);
  
  int p2_n_rows = p2.nrow();
  
  // repeat r if single value
  NumericVector r_rep = rep_len(r, p2_n_rows);
  
  // get angles between line segments and lines between p1 and p2
  NumericVector theta = atan(r_rep / d) * 180.0 / M_PI;
  
  
  NumericVector angles_ac(p2_n_rows);
  NumericVector angles_cw(p2_n_rows);
  
  for(int i = 0; i < p2_n_rows; i++) {
    angles_ac[i] = fmod(360.0 + a12[i] + theta[i], 360.0);
    angles_cw[i] = fmod(360.0 + a12[i] - theta[i], 360.0);
  }
  
  // convert angles to sine and cosine
  NumericMatrix ac = aTOd_rcpp(angles_ac);
  NumericMatrix cw = aTOd_rcpp(angles_cw);
  
  // multiply sine and cosine with with distance to get endpoints
  for(int i = 0; i < p2_n_rows; i++) {
    ac(i, _) = ac(i, _) * d[i];
    cw(i, _) = cw(i, _) * d[i];
  }
  
  NumericMatrix ac_t = transpose(ac);
  NumericMatrix cw_t = transpose(cw);
  
  // adjust for pos of p1
  for(int j = 0; j < p2.ncol(); j++) {
    ac_t(j, _) = ac_t(j, _) + p1[j];
    cw_t(j, _) = cw_t(j, _) + p1[j];
  }
  
  // prepare output list
  List dim_names = List::create(
    Named("") = rownames(p2),
    Named("") = CharacterVector::create("x", "y")
  );
  NumericMatrix ac_t_t = transpose(ac_t);
  ac_t_t.attr("dimnames") = dim_names;
  NumericMatrix cw_t_t = transpose(cw_t);
  cw_t_t.attr("dimnames") = dim_names;
  
  List ac_cw = List::create(Named("ac") = ac_t_t, Named("cw") = cw_t_t);
  
  return(ac_cw);
}


// Helper function to transform end cones into cells
NumericVector fix(NumericVector x) {
  if (is_true(all(is_na(x)))) {
    NumericVector na_res = NumericVector::create(NA_REAL);
    return(na_res);
  }
  LogicalVector false_true = LogicalVector::create(false, true);
  LogicalVector true_false = LogicalVector::create(true, false);
  if (is_true(all(is_na(x) == false_true)) | is_true(all(is_nan(x) == false_true))) {
    x = NumericVector::create(x[0], 10);
  } else if (is_true(all(is_na(x) == true_false)) | is_true(all(is_nan(x) == true_false))) {
    x = NumericVector::create(0, x[1]); 
  }
  if (x.length() > 1) {
    x = seq(x[0], x[1]); 
  }
  return(x);
}


// Helper function to remove NAs from list
List check_list_na(List x) {
  LogicalVector valid(x.length());
  for(int i = 0; i < x.length(); i++) {
    NumericVector x_i = x[i];
    if (is_false(any(is_na(x_i))) & is_false(any(is_nan(x_i)))) {
      valid[i] = true;
    }
  }
  
  return(x[valid]);
}


// Helper function to remove NULLs from list
List check_list_null(List x) {
  LogicalVector valid(x.length());
  for(int i = 0; i < x.length(); i++) {
    if (x[i] != R_NilValue) {
      valid[i] = true;
    }
  }
  
  return(x[valid]);
}


// Helper function to remove zero length elements from list
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


//' Intersecting Cones
//'
//' Compute intersecting cones between pedestrian `p1` closest pedestrians `p2`
//' from the perspective of `p1`.
//' 
//' Returns `NULL` if no other pedestrians are present, if they cannot be seen,
//' if they are not in front of pedestrian `n`, or if their cones are not
//' intersecting.
//' 
//' @param p1 Numeric matrix with shape 1x2 (x and y).
//' @param a Numeric scalar heading angle of pedestrian `p1`.
//' @param p2 Numeric matrix with shape Nx2 (x and y) indicating the positions
//' of all \emph{other} pedestrians.
//' @param r Numeric vector with radii of pedestrian bodies.
//' @param objects List containing a list for each object. An object has
//' two length-two numeric vectors of x- and y-coordinates.
//' 
//' @return Named numeric vector with intersecting cones or NULL. Names
//' indicate the cone number and the elements the distances to the intersection.
//' 
// [[Rcpp::export]]
Nullable<NumericVector> iCones_rcpp(
    NumericMatrix p1,
    double a,
    NumericMatrix p2,
    NumericVector r,
    List objects
) {
  if (p2.nrow() == 0) { // if single pedestrian
    return(R_NilValue);
  }
  
  CharacterVector p2_row_names = rownames(p2);
  p2_row_names.names() = p2_row_names;
  
  // calc endpoints of intersecting cones
  List ends = eObjects_rcpp(p1, p2, r);
  NumericMatrix ends_ac = ends[0];
  NumericMatrix ends_cw = ends[1];
  
  NumericMatrix endCones(p2.nrow(), ends.length());
  rownames(endCones) = p2_row_names;
  
  // calc cones of endpoints
  for(int i = 0; i < ends.length(); i++) {
    endCones(_, i) = Iangle_rcpp(p1, a, ends[i]) - 1; // align with cpp indexing
  }
  
  List cList(p2.nrow());
  cList.names() = p2_row_names;
  
  // adjust extreme cones (most left and right)
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
  
  double v = 1.0;
  
  // calc cell centers for p1
  NumericMatrix coneLineEnds = c_vd_rcpp(
    seq_len(11), p1, v, a, get_vels(), get_angles(), 0.5
  );
  
  CharacterVector cList_names = cList.names();
  
  // init distance to cone list
  List cDist(cList.length());
  
  for(int i = 0; i < cList.length(); i++) { // for each other pedestrian
    NumericVector cList_i = cList[i];
    String name_i = cList_names[i];
    int idx = p2_row_names.offset(name_i); // get name pedestrian
    
    if (cList_i.length() == 1) {
      // check if p1 can see pedestrian
      bool sees_goal = seesGoal_rcpp(p1, p2(idx, _), objects);
      
      if (!sees_goal) {
        cList[i] = NumericVector::create();
      } else {
        cDist[i] = dist1_rcpp(p1, p2(Range(idx, idx), _));
      }
    } else {
      // for each intersecting cone
      for(int j = 0; j < cList_i.length(); j++) {
        // calc intersection point
        NumericVector P_n_vec = line_line_intersection_rcpp(
          p1, coneLineEnds(cList_i[j], _), ends_ac(idx, _), ends_cw(idx, _)
        );
        // check if p1 can see intersection
        bool sees_goal = seesGoal_rcpp(p1, P_n_vec, objects);
        if (!sees_goal) {
          cList_i[j] = NA_REAL;
        } else {
          P_n_vec.attr("dim") = Dimension(1, 2);
          NumericMatrix P_n = as<NumericMatrix>(P_n_vec);
          // calc dist to intersection
          if (j == 0) { // init dist vector
            cDist[i] = dist1_rcpp(p1, P_n);
          } else { // if dist vector is still emtpy
            Nullable<NumericVector> _cDist_i = cDist[i];
            if (_cDist_i == R_NilValue) {
              cDist[i] = dist1_rcpp(p1, P_n);
            } else {
              NumericVector cDist_i = as<NumericVector>(_cDist_i);
              NumericVector d = dist1_rcpp(p1, P_n);
              cDist_i.push_back(d[0]); // add dist to vector
              cDist[i] = cDist_i;
            }
          }
        }
      }
      cList[i] = cList_i[!is_na(cList_i)]; // remove NA dist entries
    }
  }
  
  cList = check_list_length(cList);
  
  if (cList.length() == 0) {
    return(R_NilValue);
  }
  
  cDist = check_list_null(cDist);
  
  // assign cone indices to dist list names
  for(int i = 0; i < cList.length(); i++) {
    NumericVector cDist_i = cDist[i];
    cDist_i.names() = int2char(cList[i]);
    cDist[i] = cDist_i;
  }
  
  IntegerVector outCones;
  NumericVector out;
  
  // get minimum dist to intersection for each cone
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
  
  out.names() = int2char(outCones + 1); // align with R indexing
  
  return(out);
}


//' Transform Intersecting Cones to Cells
//'
//' Transform intersecting cones into vector with cell names containing
//' distances to each cell ring for a pedestrian moving with speed `v`.
//' 
//' @param iC Numeric vector with intersecting cones returned by
//' \link[=iCones]{iCones}.
//' @param v Numeric scalar velocity of current pedestrian.
//' @param tStep Numeric scalar velocity scaling factor.
//' 
//' @return Numeric vector with distances.
// [[Rcpp::export]]
NumericVector iCones2Cells_rcpp(NumericVector iC, double v, double tStep = 0.5) {
  NumericVector vels = NumericVector::create(1.5, 1.0, 0.5);
  
  NumericVector out = rep(iC, 3) - rep_each(v * tStep * vels, iC.length());
  
  IntegerVector rings = IntegerVector::create(0, 11, 22);
  out.names() = int2char(rep(char2int(iC.names()), 3) + 
    rep_each(rings, iC.length()));
  
  return(out);
}


// Helper function to get row n from matrix p_mat returned as 1x2 matrix
NumericMatrix get_p1(int n, NumericMatrix p_mat) {
  NumericVector p1_vec = p_mat(n, _);
  p1_vec.attr("dim") = Dimension(1, 2);
  NumericMatrix p1 = as<NumericMatrix>(p1_vec);
  CharacterVector p_mat_rownames = rownames(p_mat);
  rownames(p1) = as<CharacterVector>(p_mat_rownames[n]);
  
  return(p1);
}


//' Blocked Angle
//'
//' Compute the distance from each cell to closest pedestrians.
//' 
//' @param p1 Numeric matrix with shape 1x2 (x and y) indicating the current pedestrian position.
//' @param a1 Double scalar angle of current pedestrian.
//' @param v1 Double scalar velocity of current pedestrian.
//' @param p2 Numeric matrix with shape Nx2 (x and y) as predicted positions
//' of all other pedestrians excluding the current pedestrian.
//' @param objects List containing a list for each object. An object has
//' two length-two numeric vectors of x- and y-coordinates.
//' 
//' @return Numeric matrix with rows for each other pedestrian and columns for
//' each cell (Nx33) or NULL.
//' 
// [[Rcpp::export]]
NumericVector blockedAngle_rcpp(NumericMatrix p1, double a1, double v1, NumericMatrix p2, NumericVector r, List objects) {
  // n = n - 1; // c++ indexing
  // NumericMatrix p_mat = state["p"];
  // NumericMatrix p1 = get_p1(n, p_mat);
  // NumericVector a = state["a"];
  // NumericMatrix p2 = omit_rows(p_pred, IntegerVector::create(n));
  // NumericVector r = state["r"];
  
  Nullable<NumericVector> _iC = iCones_rcpp(p1, a1, p2, r, objects);
  if (_iC == R_NilValue) { // if iCones is NULL return empty vector
    NumericVector out_empty;
    return(out_empty);
  }
  NumericVector iC = as<NumericVector>(_iC);
  // NumericVector v = state["v"];
  NumericVector cells = iCones2Cells_rcpp(iC, v1);
  
  return(cells);
}


//' Leaders
//'
//' Get leaders for pedestrian `n`.
//' 
//' Returns `NULL` if no other pedestrians are present, if they cannot be seen,
//' if they are not in front, or if no other pedestrian has the same group as
//' pedestrian `n` and `onlyGroup` is `TRUE`.
//' 
//' @param n Integer scalar index of current pedestrian in pedestrian matrix.
//' @param p_mat Numeric matrix with shape Nx2 (x and y) indicating the positions of all pedestrians.
//' @param a Numeric vector with angles of all pedestrians.
//' @param v Numeric vector with velocities of all pedestrians.
//' @param P1 Numeric matrix with shape 1x2 (x and y) indicating position of current goal of current pedestrian.
//' @param group Numeric vector with group membership indices of all pedestrians.
//' @param centres Numeric matrix with x and y for each cell centre (33x2).
//' @param objects List containing a list for each object. An object has
//' two length-two numeric vectors of x- and y-coordinates.
//' @param onlyGroup Logical scalar indicating if leaders must be from the same
//' group.
//' @param preferGroup Logical scalar indicating if leaders are preferred from
//' the same group.
//' @param pickBest Logical salar indicating if all or the best leader should be
//' returned.
//' 
//' @return List:
//'   \describe{
//'     \item{dists}{Numeric matrix with distances for each cell (columns) to
//'       each leader (rows).}
//'     \item{leaders}{Numeric matrix with three rows for every leader
//'       (columns). The first row is the cell of the leader, the second the
//'       difference in heading angle between pedestrian `n` and leader, the
//'       third the index of the shared group.}
//'   }
//' 
// [[Rcpp::export]]
Nullable<List> getLeaders_rcpp(
  int n,
  NumericMatrix p_mat,
  NumericVector a,
  NumericVector v,
  NumericMatrix P1,
  NumericVector group,
  NumericMatrix centres,
  List objects,
  bool onlyGroup = false,
  bool preferGroup = true,
  bool pickBest = false
) {
  n = n - 1; // c++ indexing
  // NumericMatrix p_mat = state["p"];
  NumericMatrix p1 = get_p1(n, p_mat);
  // NumericVector a = state["a"];
  // NumericVector v = state["v"];
  double a1 = a[n];
  double v1 = v[n];
  // omit self from pedestrian matrix
  NumericMatrix ps = omit_rows(p_mat, IntegerVector::create(n));
  
  // init occluded vector
  IntegerVector occluded = seq_along(v) - 1; // align with cpp indexing
  occluded.erase(n); // remove self
  LogicalVector is_seen = seesMany_rcpp(p1, ps, objects);
  occluded = occluded[!is_seen]; // remove seen pedestrians
  occluded.push_front(n); // add self to front
  
  // get angles from peds that are not occluded
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
  
  // omit occluded peds from matrix
  NumericMatrix p2 = omit_rows(p_mat, occluded);
  if (p2.nrow() == 0) {
    return(R_NilValue);
  }
  
  // get goal matrix
  // NumericMatrix P1;
  // 
  // if (is<List>(state["P"])) {
  //   List P1_list = state["P"];
  //   NumericMatrix P1_mat = P1_list[n];
  //   int attr_i = P1_mat.attr("i");
  //   P1 = P1_mat(Range(attr_i - 1, attr_i - 1), Range(0, 1));
  // } else {
  //   NumericMatrix P1_mat = state["P"];
  //   P1 = P1_mat(Range(n, n), Range(0, 1));
  // }
  
  // bin angle to other peds
  NumericVector I_n = Iangle_rcpp(p1, a1, p2);
  I_n = I_n[!is_na(I_n)];
  if (I_n.length() == 0) { // if none are in front
    return(R_NilValue);
  }
  
  CharacterVector I_n_names = I_n.names();
  CharacterVector p2_rownames = rownames(p2);
  
  // get pos of peds in front
  NumericMatrix p2_ring(I_n_names.length(), p2.ncol());
  
  for(int i = 0; i < I_n_names.length(); i++) {
    for(int j = 0; j < p2_rownames.length(); j++) {
      if (I_n_names[i] == p2_rownames[j]) {
        p2_ring(i, _) = p2(j, _);
      }
    }
  }
  
  // calc distance to peds in front
  NumericVector d_ring = dist1_rcpp(p1, p2_ring);
  
  // bin distance into rings
  NumericVector bins_ring = NumericVector::create(0, 0.5, 1, 5) * v1 * 0.5;
  
  NumericVector ring = 4 - bin_vector(d_ring, bins_ring);
  ring.names() = I_n_names;
  ring = ring[!is_na(ring)];
  if (ring.length() == 0) { // if none are in rings
    return(R_NilValue);
  }
  
  // init candidates for leaders
  CharacterVector ring_names = ring.names();
  NumericVector candidates = I_n[ring_names];
  candidates = candidates + 11 * (ring - 1);
  
  CharacterVector candidates_names = candidates.names();
  // NumericVector group = state["group"];
  NumericVector group_candidates = group[candidates_names];
  LogicalVector inGroup = group_candidates == group[n];
  inGroup.names() = group_candidates.names();
  
  if (onlyGroup) { // if leaders must be in same group
    if (is_false(any(inGroup))) {
      return(R_NilValue);
    } else {
      candidates = candidates[inGroup];
      candidates_names = candidates_names[inGroup];
    } 
  // if leaders are preferred from in-group
  } else if (preferGroup && is_true(any(inGroup))) {
    candidates = candidates[inGroup];
    candidates_names = candidates_names[inGroup];
  }
  
  // get candidate angles
  NumericVector a2_candidates = a2[candidates_names];
  
  // get angles between p1 and goal
  NumericVector d_n = Dn_rcpp(p1, P1);
  NumericVector angles(a2_candidates.length());
  
  // get minimum angular difference between candidates and p1 to goal
  for(int i = 0; i < a2_candidates.length(); i++) {
    NumericVector a2_canidates_min_angle = minAngle_rcpp(a2_candidates[i], d_n);
    angles[i] = a2_canidates_min_angle[0];
  }
  angles.names() = candidates_names;
  
  // check if difference is below 90 deg
  LogicalVector ok = angles < 90.0;
  if (is_false(any(ok))) {
    return(R_NilValue);
  }
  
  candidates = candidates[ok];
  candidates_names = candidates_names[ok];
  angles = angles[ok];
  NumericVector leaders;
  CharacterVector leaders_names;
  if(is_false(any(duplicated(candidates)))) {
    leaders = candidates;
    leaders_names = candidates_names;
  } else {
    leaders = unique(candidates);
    leaders_names = rep("", leaders.length());
    for(int i = 0; i < leaders.length(); i++) {
      double leaders_i = leaders[i];
      NumericVector leaders_angles_i = angles[candidates == leaders_i];
      double min_angle_i = which_min(leaders_angles_i);
      CharacterVector candidates_names_i = candidates_names[candidates == leaders_i];
      leaders_names[i] = candidates_names_i[min_angle_i];
    }
    angles = angles[leaders_names];
  }
  
  // get distance from cell centres to leaders
  NumericMatrix d(leaders.length(), 33);
  rownames(d) = leaders_names;
  
  for(int i = 0; i < leaders.length(); i++) {
    double leaders_i = leaders[i];
    NumericVector centres_i = centres(leaders_i - 1, _);
    d(i, _) = dist1_rcpp(centres_i, centres);
  }
  
  // pick closest leader?
  IntegerVector best;
  if (pickBest) {
    best = which_min(angles);
  } else {
    best = seq_along(angles) - 1;
  }
  
  // prepare output
  NumericVector leaders_best = leaders[best];
  CharacterVector leaders_names_best = leaders_names[best];
  
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


//' Buddies
//'
//' Compute the distance from each cell to buddies.
//' 
//' Returns `NULL` if no other pedestrians are present, if they cannot be seen,
//' if they are not in front, or if no other pedestrian has the same group as
//' pedestrian `n`.
//' 
//' @param n Integer scalar index of current pedestrian in pedestrian matrix.
//' @param p_mat Numeric matrix with shape Nx2 (x and y) indicating the positions of all pedestrians.
//' @param v Numeric vector with velocities of all pedestrians.
//' @param group Numeric vector with group membership indices of all pedestrians.
//' @param a Numeric vector with angles of all pedestrians.
//' @param p_pred Numeric matrix with shape Nx2 (x and y) as predicted positions
//' of all pedestrians.
//' @param centres Numeric matrix with x and y for each cell centre (33x2).
//' @param objects List containing a list for each object. An object has
//' two length-two numeric vectors of x- and y-coordinates.
//' @param pickBest Logical salar indicating if all or the best buddy should be
//' returned.
//' @param state List with state data.
//' 
//' @return List:
//'   \describe{
//'     \item{leaders}{Numeric matrix with two rows for every buddy
//'       (columns). The first row is the cell of the leader, the second the
//'       difference in heading angle between pedestrian `n` and buddy.}
//'     \item{dists}{Numeric matrix with distances for each cell (columns) to
//'       each buddy (rows).}
//'   }
//' 
// [[Rcpp::export]]
Nullable<List> getBuddy_rcpp(
  int n,
  NumericMatrix p_mat,
  NumericVector v,
  NumericVector group,
  NumericVector a,
  NumericMatrix p_pred,
  NumericMatrix centres,
  List objects,
  bool pickBest
) {
  n = n - 1; // c++ indexing
  // NumericMatrix p_mat = state["p"];
  NumericMatrix p1 = get_p1(n, p_mat);
  // NumericVector v = state["v"];
  // omit self from pedestrian matrix
  NumericMatrix ps = omit_rows(p_mat, IntegerVector::create(n));
  
  // init occluded vector
  IntegerVector occluded = seq_along(v) - 1; // align with cpp indexing
  occluded.erase(n); // remove self
  LogicalVector is_seen = seesMany_rcpp(p1, ps, objects);
  occluded = occluded[!is_seen]; // remove seen peds
  occluded.push_front(n); // add self to front
  
  // omit occluded from predicted pos matrix
  p_pred = omit_rows(p_pred, occluded);
  if (p_pred.nrow() == 0) {
    return(R_NilValue);
  }
  
  // get visible group members
  CharacterVector group_names = group.names();
  NumericVector group_visible;
  CharacterVector inGroup_names;
  for(int i = 0; i < group.length(); i++) {
    if (is_false(any(occluded == i))) {
      group_visible.push_back(group[i]);
      inGroup_names.push_back(group_names[i]);
    }
  }
  
  LogicalVector inGroup = group_visible == group[n];
  inGroup.names() = inGroup_names;
  IntegerVector not_in_group_idx = seq_along(inGroup) - 1;
  not_in_group_idx = not_in_group_idx[!inGroup];
  p_pred = omit_rows(p_pred, not_in_group_idx);
  int nped = p_pred.nrow();
  if (nped == 0) { // if none visible
    return(R_NilValue);
  }
  
  CharacterVector p_pred_row_names = rownames(p_pred);
  
  p_pred_row_names.names() = p_pred_row_names;
  NumericVector a_p_pred = a[p_pred_row_names];
  // angular difference between p1 and predicted pos for each cone
  NumericMatrix headingDifference_t = headingAngle_rcpp(a_p_pred, a[n]);
  NumericMatrix headingDifference = transpose(headingDifference_t);
  
  // get minimum angular difference for each cone
  IntegerVector parallelCone(headingDifference.ncol());
  for(int i = 0; i < parallelCone.length(); i++) {
    parallelCone[i] = which_min(headingDifference(_, i));
  }
  CharacterVector parallelCone_names = colnames(headingDifference);
  parallelCone.names() = parallelCone_names;
  
  // get ring of group members with minimal distance
  NumericMatrix d_rings(3, nped);
  IntegerVector ring(d_rings.ncol());
  for(int i = 0; i < nped; i++) { // for each group member
    int parallelCone_i = parallelCone[i];
    String name_i = parallelCone_names[i]; // get name of group member
    int name_idx = p_pred_row_names.offset(name_i);
    // get predicted pos
    NumericVector p_pred_parallel = p_pred(name_idx, _);
    IntegerVector centres_idx = IntegerVector::create(
      parallelCone_i, parallelCone_i + 11, parallelCone_i + 22
    );
    // get cell centres
    NumericMatrix centres_i(centres_idx.length(), centres.ncol());
    for(int j = 0; j < centres_idx.length(); j++) {
      centres_i(j, _) = centres(centres_idx[j], _);
    }
    // calc distance from each cell centres to predicted pos
    d_rings(_, i) = dist1_rcpp(p_pred_parallel, centres_i);
    
    // get ring with minimal distance
    ring[i] = which_min(d_rings(_, i));
  }
  
  IntegerVector parallelCone_11 = parallelCone + 11;
  IntegerVector parallelCone_22 = parallelCone + 22;
  IntegerMatrix cells = cbind(parallelCone, parallelCone_11, parallelCone_22);
  // get cells of ring with minimal distance
  NumericVector cell(nped);
  for(int i = 0; i < nped; i++) {
    cell[i] = cells(i, ring[i]);
  }
  
  cell.names() = parallelCone_names;
  
  // get angular difference to group members
  NumericVector angleDisagree(headingDifference.ncol());
  for(int i = 0; i < headingDifference.ncol(); i++) {
    angleDisagree[i] = headingDifference(parallelCone[i], i) / 90.0;
  }
  
  // get distance from each cell to group members
  NumericMatrix d_buddy(nped, centres.nrow());
  for(int i = 0; i < nped; i++) {
    NumericVector centres_cell_i = centres(cell[i], _);
    d_buddy(i, _) = dist1_rcpp(centres_cell_i, centres);
  }
  
  // pick closest buddy?
  IntegerVector best;
  if (pickBest) {
    best = which_min(angleDisagree);
  } else {
    best = seq_len(nped) - 1;
  }
  
  // prepare output
  NumericMatrix buddies(2, best.length());
  NumericVector cell_best = cell[best];
  buddies(0, _) = cell_best + 1;
  NumericVector angleDisagree_best = angleDisagree[best];
  buddies(1, _) = angleDisagree_best;
  rownames(buddies) = CharacterVector::create("cell", "angleDisagree");
  colnames(buddies) = parallelCone_names[best];
  
  NumericMatrix dists(best.length(), d_buddy.ncol());
  colnames(dists) = int2char(seq_len(d_buddy.ncol()));
  for(int i = 0; i < best.length(); i++) {
    dists(i, _) = d_buddy(best[i], _);
  }
  
  List out_list = List::create(
    Named("buddies") = buddies, Named("dists") = dists
  );
  
  return(out_list);
}