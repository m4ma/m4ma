#include <Rcpp.h>
#include <math.h>
#include "geometry.h"
#include "see.h"
#include "utils.h"

using namespace Rcpp;

//' Rectangle to Lines
//'
//' Transforms a list of objects into a list of their constituent lines.
//' 
//' @param o List with length-two numeric vectors of x- and y-coordinates of an
//' object.
//' 
//' @return List with two matrices for the start and end points of the
//' constituent lines. Each matrix has two columns with x- and y-coordinates and 
//' a row for each line.
// [[Rcpp::export]]
List object2lines_rcpp(List o) {
  NumericVector x = o["x"];
  NumericVector y = o["y"];
  
  // re-arrange coordinates
  NumericVector vec_P1 = NumericVector::create(
    x[0], y[0], x[0], y[1], x[0], y[0], x[1], y[0]
  );
  NumericVector vec_P2 = NumericVector::create(
    x[1], y[0], x[1], y[1], x[0], y[1], x[1], y[1]
  );
  
  // convert vectors to matrices
  vec_P1.attr("dim") = Dimension(2, 4);
  vec_P2.attr("dim") = Dimension(2, 4);
  NumericMatrix mat_P1 = as<NumericMatrix>(vec_P1);
  NumericMatrix mat_P2 = as<NumericMatrix>(vec_P2);
  
  CharacterVector mat_rownames = CharacterVector::create("x", "y");
  CharacterVector mat_colnames = CharacterVector::create(
    "L1", "L2", "L3", "L4"
  );
  
  rownames(mat_P1) = mat_rownames;
  colnames(mat_P1) = mat_colnames;
  rownames(mat_P2) = mat_rownames;
  colnames(mat_P2) = mat_colnames;
  
  return(List::create(Named("P1") = mat_P1, Named("P2") = mat_P2));
}


// Helper function to check if the line between segments_P1 and segments_P2
// intersects with the line between P1 and P2
bool overlaps(
    NumericMatrix P1,
    NumericMatrix P2,
    NumericMatrix segments_P1,
    NumericMatrix segments_P2
) {
  for (int j = 0; j < segments_P1.ncol(); j++) { // for each segment
    for (int k = 0; k < P1.ncol(); k++) { // for each object line
      // check intersection
      if (is_true(all(is_finite(line_line_intersection_rcpp(
          segments_P1(_, j), segments_P2(_, j), P1(_, k), P2(_, k), true
      ))))) {
        return(true);     
      }
    }
  }
  return(false);
}


//' Body-object Overlap
//'
//' Checks if a body with radius `r` overlaps with different object lines at 
//' cell centre coordinates `okCentres`. These are coordinates of cells that 
//' have been labeled 'ok' previously. 
//' 
//' @param oL List with two matrices for the start and end points of the
//' constituent lines. Each matrix has two columns with x- and y-coordinates and 
//' a row for each line.
//' @param r Numeric scalar radius of the body.
//' @param okCentres Numeric matrix with x- and y-coordinates of cell centers 
//' that have been labeled 'ok'.
//' 
//' @return Logical vector indicating whether there is an overlap between the 
//' body at coordinates in `okCentres` and any object lines `oL`.
// [[Rcpp::export]]
LogicalVector bodyObjectOverlap_rcpp(
    List oL,
    double r,
    NumericMatrix okCentres
) {
  NumericMatrix P1 = oL["P1"];
  NumericMatrix P2 = oL["P2"];
  
  // calc right angles to object lines
  NumericVector a = angle2_rcpp(transpose(P1), transpose(P2)) + 90.0;
  for (int i = 0; i < a.length(); i++) {
    a[i] = fmod(a[i], 180.0); // a is always positive
  }
  
  a = a[!duplicated(a)];
  
  // calc dx and dy to move along lines
  NumericMatrix pos(2, a.length());
  
  NumericVector a_pi = a * M_PI / 180.0;
  
  pos(0, _) = r * sin(a_pi);
  pos(1, _) = r * cos(a_pi);
  
  LogicalVector out(okCentres.nrow());
  
  for (int i = 0; i < okCentres.nrow(); i++) { // for each center
    NumericVector p = okCentres(i, _);
    
    NumericMatrix segments_P1(pos.nrow(), pos.ncol());
    NumericMatrix segments_P2(pos.nrow(), pos.ncol());
    
    // add/sub center coordinates
    for (int j = 0; j < pos.nrow(); j++) {
      segments_P1(j, _) = p[j] - pos(j, _);
      segments_P2(j, _) = p[j] + pos(j, _);
    }
    
    // check overlap
    out[i] = overlaps(P1, P2, segments_P1, segments_P2);
  }
  
  return(out);
}


//' Check Cells
//'
//' Checks if cells are 'ok' or if the body radius at the cell coordinates
//'  overlaps with any objects.
//' 
//' @param r Numeric scalar radius of the body.
//' @param centres Numeric matrix with x- and y-coordinates of cell centers.
//' @param objects List containing a list for each object. An object has
//' two length-two numeric vectors of x- and y-coordinates.
//' @param ok Logical vector indicating which cells have been labeled 'ok'.
//' 
//' @return Logical matrix (11x3) indicating whether cells are 'ok' or whether 
//' there is an overlap between the body at cell coordinates and any objects.
// [[Rcpp::export]]
Nullable<LogicalMatrix> bodyObjectOK_rcpp(
    double r,
    NumericMatrix centres,
    List objects,
    LogicalVector ok
) {
  if (is_false(any(ok))) {
    return(R_NilValue);
  }
  
  IntegerVector which_ok = seq_along(ok) - 1;
  
  NumericMatrix okCentres = omit_rows(centres, which_ok[!ok]);
  
  List oLines = lapply(objects, object2lines_rcpp);
  
  LogicalMatrix body_overlap_mat(sum(ok), oLines.length());
  
  for (int i = 0; i < oLines.length(); i++) {
    body_overlap_mat(_, i) = bodyObjectOverlap_rcpp(
      oLines[i], r, okCentres
    );
  }
  
  LogicalVector out(33);
  
  for (int i = 0; i < body_overlap_mat.nrow(); i++) {
    out[ok] = is_false(any(body_overlap_mat(i, _)));
  }
  
  out.attr("dim") = Dimension(11, 3);
  
  return(as<LogicalMatrix>(out));
}