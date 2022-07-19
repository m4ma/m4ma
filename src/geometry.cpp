#include <Rcpp.h>
#include <math.h>
#include "m4ma.h"

using namespace Rcpp;

//' dist_rcpp
//'
//' Compute distance from p1 to p2, both xy column matrices
//' @param p1 Numeric matrix
//' @param p2 Numeric matrix
//' @return Named numeric vector of length equal to the number of rows in p1
//' @export
// [[Rcpp::export]]
NumericVector dist_rcpp(NumericMatrix p1, NumericMatrix p2) {
  
  // initialise variables
  int ncols = p1.ncol();
  int nrows = p1.nrow();
  NumericMatrix temp(p1);
  CharacterVector p1_names = rownames(p1);
  NumericVector out(nrows);
  out.names() = p1_names;
  
  // difference between the rows of p1 and p2 elevated to the power of 2
  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < ncols; j++){
      temp(i,j) = std::pow(p1(i,j) - p2(i,j),2);
    }
  }
  
  // apply square root of the sum of the rows of x and y
  for(int i = 0; i < nrows; i++){
    out(i) = std::sqrt(temp(i,0) + temp(i,1));
  }
  
  return out;
}


//' dist1_rcpp
//'
//' Compute distance from p1 to p2
//' @param p1 Numeric matrix of a single point (i.e., 1 row, 2 xy columns)
//' @param p2 Numeric matrix of multiple xy points
//' @return Named numeric vector of length equal to the number of rows in p2
//' @export
// [[Rcpp::export]]
NumericVector dist1_rcpp(NumericMatrix p1, NumericMatrix p2) {
  
  // allocate transposed matrix
  NumericMatrix t_p2 = transpose(p2);
  // allocate indices
  int ncols = t_p2.ncol();
  int nrows = t_p2.nrow();
  // allocate temporary matrix
  NumericMatrix temp_matrix(t_p2);
  // allocate output
  NumericVector out(ncols); 
  
  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < ncols; j++){
      temp_matrix(i,j) = std::pow(t_p2(i,j) - p1(0,i),2);
      out(j) = std::sqrt(sum(temp_matrix(_, j)));
    }
  }
  out.names() = rownames(p2);
  return out;
}


//' angle2s_rcpp
//'
//' Compute Shortest angle anti-clockwise from p1 as origin to p2 (> -180 to 180)
//' @param p1 Numeric matrix
//' @param p2 Numeric matrix
//' @return Numeric vector of length equal to the number of rows in p1 
//' @export
// [[Rcpp::export]]
NumericVector angle2s_rcpp(NumericMatrix p1, NumericMatrix p2) {
  
  // initialise variable
  NumericVector out(p1.nrow());
  CharacterVector p1_names = rownames(p1);
  int nrow = p1.nrow();
  
  for(int row = 0; row < nrow; row++){
    out(row) = (180 / M_PI) * std::atan2((p2(_,1) - p1(_,1))[row], (p2(_,0) - p1(_,0))[row]);
  }
  
  out = round(out, 10);
  
  // return a numeric vector whose names match those of p1's rows
  // only if p1's rows have names
  if(is_false(all(is_na(p1_names)))){
    out.names() = p1_names;
  }
  
  return(out);
}



//' angle2
//' 
//' Anti-clockwise angle from p1 as origin to p2 (x,y pairs matrices). The angle goes from 0 to 360.
//' @param p1 Numeric matrix
//' @param p2 Numeric matrix
//' @return Named numeric vector of length equal to the number of rows in p1
// [[Rcpp::export]]
NumericVector angle2_rcpp(NumericMatrix p1, NumericMatrix p2) {
  // initialise variable
  int n_rows = p1.nrow();
  NumericVector angle(n_rows);
  
  for(int row = 0; row < n_rows; ++row){
    double angle_centered = (180 / M_PI) * std::atan2((p2(_, 1) - p1(_, 1))[row], (p2(_, 0) - p1(_, 0))[row]);
    angle[row] = fmod(360 + angle_centered, 360);
  }
  
  // not clear here whether angle's names need to be p1 or p2' names
  angle.names() = rownames(p1);
  return angle;
}

//' aTOd
//' 
//' Compute sine and cosine of an angle.
//' @param a Numeric vector - angles in degrees between 0 and 360
//' @return Numeric Matrix of x and y coordinates (between -1 and 1) that are the signed (i.e., +/-) normalised difference between the xy points that generated the angle a
// [[Rcpp::export]]
NumericVector aTOd_rcpp(NumericVector a) {
  
  int rows = a.length();
  double radians = (M_PI / 180);
  NumericMatrix coordinates(rows,2);
  colnames(coordinates) = CharacterVector::create("x", "y");
  
  for(int i = 0; i < rows; i++){
    coordinates(i,0) = cos(a[i] * radians);
    coordinates(i,1) = sin(a[i] * radians);
  }
  
  return(coordinates); 
}

//' Iangle
//' 
//' Which angle cone (1..11, NA means outside of view) is p2 in relative to p1 heading at angel a1 
//' @param p1 Numeric matrix
//' @param p2 Numeric matrix
//' @param a1 Numeric vector - angles in degrees between 0 and 360
//' @return Numeric Vector of indices whereby the angle is within border's bins
// [[Rcpp::export]]

NumericVector Iangle_rcpp(NumericMatrix p1, NumericMatrix p2, NumericVector a1) {
  
  int rows = a1.length();
  NumericVector angle = angle2_rcpp(p1,p2);
  NumericVector a(rows);
  
  // NOTE: border was hand-coded in the original R function with the values provided below.
  NumericVector border = NumericVector::create(-85, -60, -40, -25, -15, -5, 5, 15, 
                                               25, 40, 60, 85);
  
  for(int i = 0; i < rows; i++){
    double a_temp = angle(i) - a1(i);
    if(a_temp < -180){
      a_temp = a_temp + 360;
      
    } else if (a_temp > 180){
      a_temp =  a_temp - 360; 
    }
    
    double a_bin = NA_REAL;
    for(int j = 0; j < border.length()-1; j++){
      // if a_temp is within border j (lower bound) and border j + 1 (upper bound)
      if(((a_temp - border(j)) * (a_temp - border(j+1))) <= 0){
        a_bin = 12 - (j+1);
      } 
    }
    
    a(i) = a_bin;
  }
  a.names() = rownames(p2);
  return(a);
}


//' Dn
//' 
//' Anti-clockwise angle to destination for all pedestrians. The angle goes from 0 to 360
//' @param p_n Numeric matrix of x and y coordinates
//' @param P_n Numeric matrix of x and y coordinates
//' @return Named numeric vector of length equal to the number of rows in p1
// [[Rcpp::export]]

NumericVector Dn_rcpp(NumericMatrix p_n, NumericMatrix P_n){
  // NOTE: this function may need more explanations?
  // the Dn function seems to do the same thing of angle2_rcpp
  // maybe the way the angle2 R function was previously written did not permit ciclying across rows?
  
  NumericVector out(p_n.nrow());
  out = angle2_rcpp(p_n,P_n);
  return(out);
  
}

//' minAngle
//' 
//' Shortest absolute angle between a1 and a2
//' @param a1_double scalar vector 
//' @param a2 Numeric vector 
//' @return Numeric vector of length equal to a1
// [[Rcpp::export]]

NumericVector minAngle_rcpp(double a1_double, NumericVector a2){
  // NOTE: it is not clear whether a1 has to be a double or a numeric vector like a2
  // I am opting for the former interpretation because of the way minAngles() was used in R within headingAngles()
  // i.e., headingAngles() in R included a for loop to provide a single element to minAngles's a1 parameter.
  
  NumericVector out(a2.length());
  
  for(int i = 0; i < a2.length(); i++){
    out(i) = std::min(abs(a1_double - a2(i)), abs(std::min(a1_double, a2(i)) + (360 - std::max(a1_double, a2(i)))));
    
  }
  return(out);
  
}

//' headingAngle
//' 
//' Absolute angular difference between 11 directions with zero cone having angle a1 and a2
//' @param a1 Numeric vector 
//' @param a2 Numeric vector 
//' @return Numeric matrix whose rows and columns are the same length of a2 and a1
// [[Rcpp::export]]

NumericMatrix headingAngle_rcpp(NumericVector a2, NumericVector a1){
  
  // NOTE: angles was hand-coded in the original R function with the values provided below.
  // consider whether angles needs to be ever used as a parameter
  
  NumericVector angles = NumericVector::create(72.5, 50, 32.5, 20, 10, 0, 350, 
                                               340, 327.5, 310, 287.5);
  int rows = a2.length();
  int cols = angles.length();
  
  // output 
  NumericMatrix output_angles(rows, cols); 
  // NOTE: not clear what is the correspondence between a1 and a2 names with the output's rows' and columns' names
  CharacterVector out_angles_rownames = a2.names();
  CharacterVector out_angles_colnames = a1.names();
  rownames(output_angles) = out_angles_rownames;
  colnames(output_angles) = out_angles_rownames;
  
  for(int i = 0; i < cols;i++){
    // sum a1 to all angles in angles
    double ang_temp = angles(i) + a1(i);
    // compute module
    ang_temp = fmod(360 + ang_temp, 360);
    // find the minAngle between ang_temp and a2 and append the output by columns in the output
    output_angles(_,i) = minAngle_rcpp(ang_temp, a2);
  }
  
  return(output_angles);

}

//' scaleVel
//' 
//' Scale velocity by time step (tStep)
//' @param v Numeric vector 
//' @param tStep double
//' @return Numeric vector of scaled velocity of same length of v
// [[Rcpp::export]]

NumericVector scaleVel_rcpp(NumericVector v, double tStep){
  
  // NOTE: not clear whether v is a vector, or a single value. 
  // this affects the choice of the output's type
  NumericVector scaled_v = v * tStep;
  
  return(scaled_v);
}


