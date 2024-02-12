#include <Rcpp.h>
#include <math.h>
#include "geometry.h"
#include "see.h"
#include "utils.h"
#include "objects.h"

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
  int nL = P1.ncol();
  int nC = okCentres.nrow();
  LogicalVector shorter(nL);
  LogicalVector out(nC); 
  NumericVector len_sq(nL);
  NumericMatrix len(nL,nC);
  NumericMatrix d21(2,nL);
  NumericMatrix d01x(nL,nC);
  NumericMatrix d01y(nL,nC);
  NumericMatrix x(nL,nC);
  NumericMatrix y(nL,nC);
  NumericMatrix dx(nL,nC);
  NumericMatrix dy(nL,nC);
  NumericMatrix leng(nL,nC);
  NumericMatrix param(nL,nC);

  for (int i = 0; i < nL; i++)  {
    d21(_,i) = P2(_,i) - P1(_,i);
    len_sq(i) = sum(d21(_,i)*d21(_,i));
    for (int j = 0; j < nC; j++)  {
      d01x(i,j) = okCentres(j,0) - P1(0,i);
      d01y(i,j) = okCentres(j,1) - P1(1,i);
      param(i,j) = d21(0,i)*d01x(i,j)+d21(1,i)*d01y(i,j);
    }
    if (len_sq(i)>0) param(i, _) = param(i, _)/len_sq(i);
  }
  for (int j = 0; j < nC; j++) {
    for (int i = 0; i < nL; i++) {
      if (param(i,j) < 0) {
        x(i,j) = P1(0,i);
        y(i,j) = P1(1,i); 
      } else if (param(i,j) > 1) {
        x(i,j) = P2(0,i) ;
        y(i,j) = P2(1,i); 
      } else {
        x(i,j) = P1(0,i) + param(i,j)*d21(0,i);
        y(i,j) = P1(1,i) + param(i,j)*d21(1,i);
      }
      dx(i,j) = okCentres(j,0) - x(i,j);
      dy(i,j) = okCentres(j,1) - y(i,j);
      leng(i,j) = sqrt(dx(i,j)*dx(i,j)+dy(i,j)*dy(i,j));
    }
    shorter = (r-leng(_,j))>0;
    out(j) = is_true(any(shorter));
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
  // create rtree structure from objects
  rtree_t rtree = objects_to_rtree(objects);
  
  // allocate output
  LogicalVector out(ok);
  
  // for each cell centre
  for (int i = 0; i < ok.length(); i++) {
    if (ok[i]) {
      // transform centre into point
      NumericVector centres_i = centres(i, _);
      point_t p = as<point_t>(centres_i);
      
      // allocate rtree query result
      std::vector<rtree_elem_t> values;
      
      // query nearest bounding box in rtree
      rtree.query(bgi::nearest(p, 1), std::back_inserter(values));
      
      rtree_elem_t res = values[0];
      
      // get object corresponding to nearest bounding box
      multi_polygon_t res_mpoly = as<multi_polygon_t>(objects[std::get<1>(res)]);
      
      // if distance to nearest object is smaller than circle radius cell is blocked
      float d = bg::distance(p, res_mpoly);
      
      out[i] = d >= r;
    }
  }
  
  // transform into matrix
  out.attr("dim") = Dimension(11, 3);
  
  return(as<LogicalMatrix>(out));
}


// Helper function to check if there is an unblocked line from each cell center to the agent
LogicalVector blocked_cells_rcpp(
  S4 agent, List objects, NumericMatrix centers, LogicalVector check = LogicalVector(33)
) {
  LogicalVector check_copy = clone(check);
  NumericVector agent_pos = agent.slot("center");
  
  for (int i = check.length()-1; i>= 0; i--) {
    if (i + 11 >= check.length() || check_copy[i + 11]) {
      NumericVector centers_i = centers(i, _);
      check_copy[i] = seesGoal_rcpp(centers_i, agent_pos, objects);
    }
  }
  
  return(!check_copy);
}


//' Check Cell Movement 
//'
//' Checks if cells are outside the background, and if there is
//' a line of sight between the agent's position and the cell center.
//' 
//' @param agent Object of class `agent`.
//' @param background Object class `background`.
//' @param centers Numeric matrix with x- and y-coordinates of cell centers.
//' 
//' @return Logical matrix (11x3) indicating whether cells are in the background,
//' and if there is a line of sight between the agent agent's position and the 
//' cell center.
// [[Rcpp::export]]
LogicalMatrix free_cells_rcpp(
  S4 agent, S4 background, NumericMatrix centers
) {
  S4 background_object = background.slot("shape");
  List objects = background.slot("objects");
  
  multi_polygon_t background_poly = as<multi_polygon_t>(background_object);
  
  LogicalVector out(centers.nrow());
  
  for (int i = 0; i < out.length(); i++) {
    NumericVector centers_i = centers(i, _);
    point_t p = as<point_t>(centers_i);
    out[i] = bg::within(p, background_poly);
  }
  
  if (objects.length() > 0) {
    LogicalVector blocked = blocked_cells_rcpp(agent, objects, centers);
    out = out & !blocked;
  }
  
  // transform into matrix
  out.attr("dim") = Dimension(11, 3);
  
  return(as<LogicalMatrix>(out));
}
