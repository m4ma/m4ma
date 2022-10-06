#include <Rcpp.h>
#include <math.h>
#include "m4ma.h"

using namespace Rcpp;

// IMPORTANT: The geometry C++ functions are currently directly integrated into the
// predped R code. That is, instead of the predped pp_geometry.R functions,
// the C++ functions are called. Thus, the arguments of the functions and their 
// order must not be changed for the integration to work. Moreover, the output 
// must be *exactly* the same as from the R functions, including the
// types, shapes, and names of dimensions (col and row names).


//' Matrix-to-matrix Distance
//'
//' Compute the Euclidian distance between two Nx2 matrices with the first
//' column in each matrix containing x- and the second column y-coordinates.
//' 
//' The matrices `p1` and `p2` must have the same shape.
//' 
//' @param p1 Numeric matrix with shape Nx2 (x and y).
//' @param p2 Numeric matrix with shape Nx2 (x and y).
//' 
//' @return Named numeric vector of length equal to the number of 
//' rows N in `p1`.
// [[Rcpp::export]]
NumericVector dist_rcpp(NumericMatrix p1, NumericMatrix p2) {
  
  int n_rows = p1.nrow();
  NumericMatrix temp(p1);
  NumericVector out(n_rows);
  
  // difference between the rows of p1 and p2 elevated to the power of 2
  for(int i = 0; i < n_rows; i++){
    temp(i, _) = pow(p1(i, _) - p2(i, _), 2);
  }
  
  // apply square root of the sum of the rows of x and y
  for(int i = 0; i < n_rows; i++){
    out[i] = std::sqrt(temp(i, 0) + temp(i, 1));
  }
  
  // assign row names of p1 if they exist
  if(rownames(p1) != R_NilValue) {
    out.names() = rownames(p1);
  }
  
  return out;
}


//' Vector-to-matrix Distance
//'
//' Compute the Euclidian distance between a two-element vector and a Nx2 
//' matrix with the first column/element containing x- and the second 
//' column/element y-coordinates.
//' 
//' @param p1 Numeric vector of length two (x and y).
//' @param p2 Numeric matrix with shape Nx2 (x and y).
//' 
//' @return Named numeric vector of length equal to the number of rows in `p2`.
// [[Rcpp::export]]
NumericVector dist1_rcpp(NumericVector p1, NumericMatrix p2) {
  
  // transpose matrix
  NumericMatrix p2_t = transpose(p2);
  
  int n_rows = p2_t.nrow();
  int n_cols = p2_t.ncol();
  
  NumericMatrix tmp(n_rows, n_cols);
  NumericVector out(n_cols);
  
  // difference between rows of p2_t and elements of p1 to the power 2
  for(int i = 0; i < n_rows; i++){
    tmp(i, _) = pow(p2_t(i, _) - p1[i], 2);
  }
  
  // square root of the column-wise sum 
  for(int j = 0; j < n_cols; j++){
    out[j] = std::sqrt(sum(tmp(_, j)));
  }
  
  // assign col names of p2_t if they exist
  if(colnames(p2_t) != R_NilValue) {
    out.names() = colnames(p2_t);
  }
  
  return out;
}


//' Anti-clockwise Angle
//'
//' Compute the anti-clockwise angles between two matrices with 
//' the first column in each matrix containing x- and the second column 
//' y-coordinates.
//' 
//' For `angle2s` and `angle2`, the angle is calculated from `p1` as origin to 
//' `p2`. Therefore, `p1` must be a 1x2 matrix. 
//' 
//' `Dn` calculates the angle between rows with the same index in 
//' between $>0$ and 360 degrees. `angle2s` computes the angle in between 
//' $>-180$ to 180 degrees and `angle2` in between $>0$ and 360 degrees. 
//' 
//' @param p1,p2,p_n,P_n Numeric matrices with shape Nx2 (x and y).
//' 
//' @return Named numeric vector of length equal to the number of rows in `p2`
//' or `p_n`.
// [[Rcpp::export]]
NumericVector angle2s_rcpp(NumericMatrix p1, NumericMatrix p2) {
  
  int n_rows = p2.nrow();
  NumericVector out(n_rows);
  
  // angle as arc tan from difference in x and y
  for(int i = 0; i < n_rows; i++){
    out[i] = (180 / M_PI) * std::atan2(
      (p2(i, 1) - p1[1]), (p2(i, 0) - p1[0])
    );
  }
  
  // round to ten decimals
  out = round(out, 10);
  
  // assign row names of p1 if they exist
  if(rownames(p2) != R_NilValue) {
    out.names() = rownames(p2);
  }
  
  return(out);
}


//' @rdname angle2s_rcpp
// [[Rcpp::export]]
NumericVector angle2_rcpp(NumericMatrix p1, NumericMatrix p2) {
  
  int n_rows = p2.nrow();
  NumericVector angle_centered(n_rows);
  
  // angle as arc tan from difference in x and y
  for(int i = 0; i < n_rows; i++){
    double angle = (180 / M_PI) * std::atan2(
      (p2(i, 1) - p1[1]), (p2(i, 0) - p1[0])
    );
    angle_centered[i] = fmod(360 + angle, 360);
  }
  
  // assign row names of p2 if they exist 
  if(rownames(p2) != R_NilValue) {
    angle_centered.names() = rownames(p2);
  }
  
  return angle_centered;
}


//' Sine and Cosine
//' 
//' Compute sine and cosine of angles.
//' 
//' @param a Numeric vector containing angles in degrees between 0 and 360.
//' 
//' @return Numeric matrix of x and y values (between -1 and 1) that are 
//' the signed normalized differences between the xy-points leading to
//' the angles `a`.
// [[Rcpp::export]]
NumericMatrix aTOd_rcpp(NumericVector a) {
  
  int n_rows = a.length();
  double radians = (M_PI / 180);
  NumericMatrix coordinates(n_rows,2);
  
  // sine and cosine of radians
  for(int i = 0; i < n_rows; i++){
    coordinates(i, 0) = cos(a[i] * radians);
    coordinates(i, 1) = sin(a[i] * radians);
  }
  
  // output matrix has x and y coordinates
  colnames(coordinates) = CharacterVector::create("x", "y");
  
  // assign names of a if they exist
  if(a.hasAttribute("names")) {
    CharacterVector a_names = a.names();
    rownames(coordinates) = a_names;
  }
  
  return(coordinates); 
}


// Helper function that converts angles to > -180 to 180
NumericVector tomp_rcpp(NumericVector x) {
  
  NumericVector out = ifelse(x <= -180.0, 360.0 + x,
                             ifelse(x > 180.0, x - 360.0,
                                    x));
  
  return(out);
}


// Helper function that bins angles according to borders and returns bin indices
NumericVector bin_angle(NumericVector a, NumericVector border) {
  int l = a.length();
  int b = border.length();
  NumericVector out(l);
  
  for(int i = 0; i < l; i++) {
    // set to NA if outside borders
    double idx = NA_REAL;
    
    // iterate over borders
    for(int j = 1; j < b; j++) {
      // if angle is > left border and <= right border
      if(a[i] > border[j-1] && a[i] <= border[j]) {
        idx = j;
      }
      out[i] = idx;
    }
  }
  
  return(out);
}


//' Bin Angles between Matrices
//' 
//' Compute anti-clockwise angles between two Nx2 matrices with the first
//' column in each matrix containing x- and the second column y-coordinates
//' (> 0 to 360 degrees) using \link[=angle2]{angle2} in relation to a scalar
//' angle. Bin relative angles according to a numeric vector of breaks.
//' 
//' Returns `NA` if relative angle is `NA` or outside the breaks. Excludes the
//' left and includes the right border of every bin.
//' 
//' @param p1 Numeric matrix with shape Nx2 (x and y).
//' @param a1 Numeric vector containing angles in degrees between $>0$ and 360.
//' @param p2 Numeric matrix with shape Nx2 (x and y).
//' @param border Numeric vector of angles indicating the breaks between bins
//' for different directions. Must be of length 11.
//' 
//' @return Numeric vector of bin indices.
// [[Rcpp::export]]
NumericVector Iangle_rcpp(NumericMatrix p1, double a1, NumericMatrix p2) {
  
  // NOTE: border was hand-coded in the original R function with the values provided below.
  NumericVector border = NumericVector::create(-85, -60, -40, -25, -15, -5, 5, 15, 
                                               25, 40, 60, 85);
  
  // calculate angle
  NumericVector angle = angle2_rcpp(p1, p2) - a1;
  
  // convert to -180 to 180
  NumericVector tmp = tomp_rcpp(angle);
  
  // bin angle
  NumericVector out = 12 - bin_angle(tmp, border);
  
  // assign row names of p2 if they exist
  if(rownames(p2) != R_NilValue) {
    out.names() = rownames(p2);
  }
  
  return(out);
}


//' @rdname angle2s_rcpp
// [[Rcpp::export]]
NumericVector Dn_rcpp(NumericMatrix p_n, NumericMatrix P_n){
  // NOTE: In contrast to angle2, this function computes the angle between the
  // pairwise rows of both matrices whereas angle2 treats one matrix as a vector
  // and computes the angle with the vector as the origin
  
  int p_n_rows = p_n.nrow();
  int p_n_cols = p_n.ncol();
  int P_n_cols = P_n.ncol();
  NumericVector out(p_n_rows);
  
  for(int i = 0; i < p_n_rows; i++){
    // create 1xN NumericMatrices required by angle_2
    NumericMatrix p_n_i(1, p_n_cols);
    p_n_i(0, _) = p_n(i, _);
    
    NumericMatrix P_n_i(1, P_n_cols);
    P_n_i(0, _) = P_n(i, _);
    
    // compute angle between rows of *both* matrices
    NumericVector angle = angle2_rcpp(p_n_i, P_n_i);
    
    out[i] = angle[0];
  }
  
  return(out);
  
}

//' Minimum Angle
//' 
//' Compute smallest absolute angles between a scalar angle and an angle
//' vector.
//' 
//' @param a1,a1_double Numeric scalar angle in between $>0$ and 360 degrees.
//' @param a2 Numeric vector of angles in between $>0$ and 360 degrees.
//' 
//' @return Numeric vector of length equal to `a2` with angles in between 
//' $>0$ and 360 degrees.
// [[Rcpp::export]]
NumericVector minAngle_rcpp(double a1_double, NumericVector a2){
  
  int l = a2.length();
  NumericVector out(l);
  
  // compute parallel minima and maxima
  for(int i = 0; i < l; i++){
    out[i] = std::min(
      abs(a1_double - a2[i]),
      abs(std::min(a1_double, a2[i]) + (360 - std::max(a1_double, a2[i])))
    );
  }
  
  // assign names of a2 if they exist
  if(a2.hasAttribute("names")) {
    CharacterVector a2_names = a2.names();
    out.names() = a2_names;
  }
  
  return(out);
  
}

//' Heading Angle
//' 
//' Compute absolute angular differences between 11 standard angles and a 
//' vector of angles.
//' 
//' Angle `a1` is added to the standard angles.
//' 
//' @param a2 Numeric vector of angles between >0 and 360 degrees.
//' @param a1 Numeric scalar angles between >0 and 360 degrees.
//' @param angles Numeric vector indicating the angle of different directions.
//' Must be of length 11.
//' 
//' @return Numeric matrix of shape Nx11 where N is the length of `a2`.
// [[Rcpp::export]]
NumericMatrix headingAngle_rcpp(NumericVector a2, double a1){
  
  // NOTE: angles was hand-coded in the original R function with the values 
  // provided below. Consider whether angles needs to be ever used as 
  // a parameter
  
  NumericVector angles = NumericVector::create(72.5, 50, 32.5, 20, 10, 0, 350, 
                                               340, 327.5, 310, 287.5);
  int n_rows = a2.length();
  int n_cols = angles.length();
  
  NumericMatrix output_angles(n_rows, n_cols);
  
  for(int i = 0; i < n_cols; i++){
    // add a1 to default angles
    double ang_temp = angles[i] + a1;
    // compute modulo
    ang_temp = fmod(360.0 + ang_temp, 360.0);
    // find the shortest absolute angle between ang_temp and a2
    output_angles(_, i) = minAngle_rcpp(ang_temp, a2);
  }
  
  // assign names of a2 if they exist
  if(a2.hasAttribute("names")) {
    CharacterVector a2_names = a2.names();
    rownames(output_angles) = a2_names;
  }
  
  return(output_angles);
  
}

//' Scaled Velocity
//' 
//' Scale velocities by a factor.
//' 
//' @param v Numeric vector with velocities.
//' @param tStep Numeric scalar scaling factor.
//' 
//' @return Numeric vector with scaled velocities of same length as `v`.
// [[Rcpp::export]]
NumericVector scaleVel_rcpp(NumericVector v, double tStep = 0.5) {
  
  NumericVector scaled_v = v * tStep;
  
  if(v.hasAttribute("names")) {
    CharacterVector v_names = v.names();
    scaled_v.names() = v_names;
  }
  
  return(scaled_v);
}



//' Cell Centers
//' 
//' Calculate centers of `cells` (index 1 to 33) given a xy-point `p1`
//' moving at velocity `v1` with angle `a1`.
//' 
//' @param cells Integer vector indicating cell indices between 1 and 33.
//' @param p1 Numeric vector of length two (x and y).
//' @param v1 Numeric scalar velocity.
//' @param a1 Numeric scalar angle.
//' @param vels Numeric matrix (33x3) of velocities for each cell.
//' @param angles Numeric matrix (33x3) of angles for each cell.
//' @param tStep Numeric scalar velocity scaling factor.
//' 
//' @return Numeric matrix (Nx2) of xy-coordinates for centers of each cell.
//' @export
// [[Rcpp::export]]
NumericMatrix c_vd_rcpp(IntegerVector cells, NumericVector p1, NumericVector v1,
                        double a1, NumericMatrix vels, NumericMatrix angles,
                        double tStep = 0.5) {
  
  IntegerVector cells_0 = cells - 1;
  
  int n_rows = cells.length();
  int n_cols = p1.length();
  NumericMatrix cell_centres(n_rows, n_cols);

  NumericVector cell_angles(n_rows);
  NumericVector cell_vels(n_rows);
  
  for(int i = 0; i < n_rows; i++) {
    // compute scaled velocity with default parameter 0.5
    cell_vels[i] = tStep * vels[cells_0[i]];
    // compute modulo
    cell_angles[i] = fmod(360.0 + (angles[cells_0[i]] + a1), 360.0);
  }
  
  // compute cosine and sine
  NumericMatrix cos_sin = aTOd_rcpp(cell_angles);
  
  // multiply cosine and sine with velocities
  for(int i = 0; i < n_rows; i++) {
    cell_centres(i, _) = cell_vels[i] * cos_sin(i, _);
  }
  
  // transpose 
  NumericMatrix cell_centres_t = transpose(cell_centres);
  
  // add p1 coordinates to centers
  for(int j = 0; j < n_cols; j++) {
    cell_centres_t(j, _) = cell_centres_t(j, _) + p1[j];
  }
  
  // transpose back
  NumericMatrix cell_centres_t_t = transpose(cell_centres_t);
  
  // assign names of p1 if they exist
  if(p1.hasAttribute("names")) {
    CharacterVector p1_names = p1.names();
    colnames(cell_centres_t_t) = p1_names;
  }
  
  return(cell_centres_t_t);
}


//' Cone Number
//'
//' @param k Numeric vector between 1 and 33.
//' @return Numeric vector of length equal to k with cone numbers between 
//' 1 and 11.
// [[Rcpp::export]]
NumericVector coneNum_rcpp(NumericVector k){
  int k_len = k.length();
  NumericVector cone_number(k_len);
  
  for(int i = 0; i < k_len; i++){
    cone_number(i) = 1 + fmod(k[i]-1, 11);
  }
  
  return(cone_number);
}


//' Ring Number
//' 
//' @param k Numeric vector between 1 and 33.
//' @return Numeric vector of length equal to k with ring numbers between 
//' 1 and 3.
// [[Rcpp::export]]
NumericVector ringNum_rcpp(NumericVector k){
  int k_len = k.length();
  NumericVector ring_number(k_len);

  for(int i = 0; i < k_len; i++){
    ring_number(i) = std::floor(1 + (k[i]-1) / 11);
  }

  return(ring_number);
}


NumericMatrix get_vels() {
  NumericMatrix vels(11, 3);
  
  for(int i = 0; i < vels.nrow(); i++) {
    vels(i, _) = NumericVector::create(1.5, 1.0, 0.5);
  }
  
  return(vels);
}


NumericMatrix get_angles() {
  NumericMatrix angles(11, 3);
  
  for(int i = 0; i < angles.ncol(); i++) {
    angles(_, i) = NumericVector::create(
      72.5, 50, 32.5, 20, 10, 0, 350, 340,327.5, 310, 287.5
    );
  }
  
  return(angles);
}