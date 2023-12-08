// [[Rcpp::depends(BH)]]

// Enable C++14 via this plugin to compile boost
// [[Rcpp::plugins("cpp14")]]

#include <Rcpp.h>
#include <math.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>

namespace bg = boost::geometry;
namespace trans = boost::geometry::strategy::transform;

namespace Rcpp {

typedef bg::model::d2::point_xy<double> point_t;
typedef bg::model::polygon<point_t> polygon_t;
typedef bg::model::multi_polygon<polygon_t> multi_polygon_t;
typedef boost::polygon::rectangle_data<double> rect;

template <> point_t as(SEXP ptsexp) {
  Rcpp::NumericVector pt(ptsexp);
  return point_t(pt[0], pt[1]);
};

template <> SEXP wrap(const point_t &p) {
  return Rcpp::wrap(Rcpp::NumericVector::create(p.x(), p.y()));
}

// Template for rectangle objects
template <> rect as(SEXP rsexp) {
  Rcpp::S4 b(rsexp);
  
  NumericVector p1 = b.slot("lower");
  NumericVector p2 = b.slot("upper");
  
  rect r = boost::polygon::construct<rect>(p1[0], p1[1], p2[0], p2[1]);
  
  rect r_t;
  
  double angle = b.slot("orientation");
  
  trans::rotate_transformer<boost::geometry::radian, double, 2, 2> rotate(-angle);
  boost::geometry::transform(r, r_t, rotate);
  
  return r_t;
};

// Template for circle objects
template <> multi_polygon_t as(SEXP csexp) {
  Rcpp::S4 b(csexp);
  
  NumericVector center = b.slot("center");
  double radius = b.slot("radius");
  
  point_t p = as<point_t>(center);
  
  multi_polygon_t circle;
  
  const int n_points = 36;

  bg::strategy::buffer::point_circle point_strategy(n_points);
  bg::strategy::buffer::distance_symmetric<double> distance_strategy(radius);
  bg::strategy::buffer::join_round join_strategy(n_points);
  bg::strategy::buffer::end_round end_strategy(n_points);
  bg::strategy::buffer::side_straight side_strategy;

  bg::buffer(p, circle, distance_strategy, side_strategy, join_strategy, end_strategy, point_strategy);
  
  return circle;
};

}

using namespace Rcpp;

//' Two-line Intersection
//' 
//' Calculate the intersection between two lines.
//' 
//' The first line is defined by \code{P1} and \code{P1}, the second by 
//' \code{P3} \code{P4}. Returns \code{Inf} when the lines are parallel. If
//' \code{interior_only} is \code{TRUE}, returns \code{NA} if the intersection
//' is not within the plane spanning between the two lines.
//' 
//' @param P1,P2,P3,P4 Numeric vector with x- and y-coordinates.
//' @param interior.only,interior_only Logical scalar indicating whether the intersection must
//' be within the spans of the two lines.
//' 
//' @returns A numeric vector with x- and y-coordinates defining the
//' intersecting point between the two lines.
//' 
//' @source Weisstein, Eric W. "Line-Line Intersection.
//' "From MathWorld--A Wolfram Web Resource.
//' \url{http://mathworld.wolfram.com/Line-LineIntersection.html}
//' 
//' @export
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
  
  // Determinant of diff matrix
  double D = dx1 * dy2 - dy1 * dx2;
  
  // If lines are parallel
  if (D == 0.0) {
    NumericVector inf_res = {R_PosInf, R_PosInf};
    return(inf_res);
  }
  
  // Determinant of x and y matrix
  double D1 = P1[0] * P2[1] - P1[1] * P2[0];
  double D2 = P3[0] * P4[1] - P3[1] * P4[0];
  
  // Intersection
  double X = (D1 * dx2 - dx1 * D2) / D;
  double Y = (D1 * dy2 - dy1 * D2) / D;
  
  if (interior_only) {
    // Is intersection within span of each line?
    double lambda1 = -((X - P1[0]) * dx1 + (Y - P1[1]) * dy1)/(pow(dx1, 2) + pow(dy1, 2));
    double lambda2 = -((X - P3[0]) * dx2 + (Y - P3[1]) * dy2)/(pow(dx2, 2) + pow(dy2, 2));
    
    if (!(isgreater(lambda1, 0.0) & isless(lambda1, 1.0) & 
        isgreater(lambda2, 0.0) & isless(lambda2, 1.0))) {
      NumericVector na_res = {NA_REAL, NA_REAL};
      return(na_res);
    }
  }
  
  NumericVector res = NumericVector::create(X, Y);
  
  return(res);
}

//' Goal in Sight
//' 
//' Checks whether a goal \code{P} can be seen from point \code{p}, or if it is
//' occluded by \code{objects}.
//'
//' @param p_n,P_n Numeric vector with x- and y-coordinates.
//' @param objects List containing a list for each object. An object has
//' two length-two numeric vectors of x- and y-coordinates.
//' 
//' @returns \code{TRUE} if the goal is in sight, \code{FALSE} otherwise.
//' @examples
//' objects = list(
//'   list(x = c(0.5, 0.5), y = c(0.5, 0.5))
//' )
//' 
//' seesGoal_rcpp(c(0, 0), c(1, 1), objects)
//' # FALSE
//' 
//' @export
// [[Rcpp::export]]
bool seesGoal_rcpp(
    NumericVector p_n,
    NumericVector P_n,
    List objects
) {
  // define new point type
  typedef bg::model::d2::point_xy<double> point_t;
  
  // create linestring object from p_n to P_n
  bg::model::linestring<point_t> l_goal;
  
  bg::append(l_goal, as<point_t>(p_n));
  bg::append(l_goal, as<point_t>(P_n));
  
  // check for each object if intersects with line
  for (int i = 0; i < objects.length(); i++) {
    RObject objects_i = objects[i];
    
    if (objects_i.isS4() && objects_i.inherits("rectangle")) {
      rect new_obj = as<rect>(objects_i);
      
      // special case: when object has area == 0 (is point), don't return false
      if (bg::intersects(l_goal, new_obj) && bg::area(new_obj) > 0.0) {
        return(false);
      }
    } else if (objects_i.isS4() && objects_i.inherits("circle")) {
      // create circle as polygon with 36 points
      multi_polygon_t new_circle = as<multi_polygon_t>(objects_i);
      
      if (bg::intersects(l_goal, new_circle) && bg::area(new_circle) > 0.0) {
        return(false);
      }
    } else if (is<List>(objects_i)) {
      List new_list = as<List>(objects_i);
      NumericVector x = new_list["x"];
      NumericVector y = new_list["y"];
      // create object box from min_corner, max_corner
      bg::model::box<point_t> new_box(
          point_t(x[0], y[0]),
          point_t(x[1], y[1])
      );
      if (bg::intersects(l_goal, new_box) && bg::area(new_box) > 0.0) {
        return(false);
      }
    }
  }
  
  return(true);
}

// Helper function to get x and y from goal n
NumericVector get_state_P_n(int n, List state, int offset=0) {
  List state_P = state["P"];
  
  NumericMatrix state_P_n = state_P[n];
  
  int i = state_P_n.attr("i");
  
  NumericVector P_n = state_P_n(i + offset - 1, _);
  IntegerVector idx = IntegerVector::create(0, 1);
  
  return(P_n[idx]);
}

//' Current Goal in Sight
//' 
//' Checks whether the current goal in a \code{state} can be seen by 
//' subject \code{n}, or if it is occluded by \code{objects}.
//'
//' @param n Integer scalar subject index.
//' @param state List of list with state data.
//' @param objects List containing a list for each object. An object has
//' two length-two numeric vectors of x- and y-coordinates.
//' @param offset Integer scalar offset to be added to the current goal index
//' (default 0).
//' 
//' @returns \code{TRUE} if the goal is in sight, \code{FALSE} otherwise.
//' @examples
//' objects = list(
//'   list(x = c(0.5, 0.5), y = c(0.5, 0.5))
//' )
//' 
//' state = list(
//'   p = matrix(c(0, 0), 1, 2),
//'   P = list(
//'     matrix(c(1, 1), 1, 2)
//'   )
//' )
//' 
//' attr(state$P[[1]], "i") = 1
//' 
//' seesCurrentGoal_rcpp(1, state, objects)
//' # FALSE
//' 
//' @export
// [[Rcpp::export]]
bool seesCurrentGoal_rcpp(int n, List state, List objects, int offset=0) {
  NumericMatrix state_p = state["p"];
  NumericVector p_n = state_p(n - 1, _);
  NumericVector P_n = get_state_P_n(n - 1, state, offset);
  
  bool sees_goal = seesGoal_rcpp(p_n, P_n, objects);
  
  return(sees_goal);
}


//' Multiple Goals in Sight
//' 
//' Checks which of the goals \code{ps} can be seen from point \code{p1},
//' or if they occluded by \code{objects}.
//'
//' @param p1 Numeric vector with x- and y-coordinates.
//' @param ps Numeric Matrix with a row for every goal and x- and y-coordinates
//' as columns.
//' @param objects List containing a list for each object. An object has
//' two length-two numeric vectors of x- and y-coordinates.
//' 
//' @returns Logical vector with \code{TRUE} if a goal is in sight, \code{FALSE}
//' otherwise.
//' @examples
//' objects = list(
//'   list(x = c(0.5, 0.5), y = c(0.5, 0.5))
//' )
//' 
//' goals = rbind(
//'   c(1, 0),
//'   c(1, 1)
//' )
//' 
//' seesGoal_rcpp(c(0, 0), goals, objects)
//' # TRUE FALSE
//' 
//' @export
// [[Rcpp::export]]
LogicalVector seesMany_rcpp(NumericVector p1, NumericMatrix ps, List objects) {
  int n_rows = ps.nrow();
  
  LogicalVector sees_many(n_rows);
  
  for (int i = 0; i < n_rows; ++i) {
    sees_many[i] = seesGoal_rcpp(ps(i, _), p1, objects);
  }
  
  return(sees_many);
}

//' Goal in Sight from Cell Centres
//' 
//' Checks whether the current goal in a \code{state} can be seen by 
//' subject \code{n} from cell \code{centres} that are marked as \code{ok},
//'  or if it is occluded by \code{objects}.
//'
//' @param n Integer scalar subject index.
//' @param objects List containing a list for each object. An object has
//' two length-two numeric vectors of x- and y-coordinates.
//' @param state List of list with state data.
//' @param centres Numeric matrix with 33 cell centres as rows and 
//' x- and y-coordinates as columns.
//' @param ok Logical vector of length 33 indicating which cells are 
//' marked as 'ok'.
//' 
//' @returns Logical vector of length 33 indicating from which 'ok' cells the 
//' current goal can be seen.
//' @examples
//' objects = list(
//'   list(x = c(0.5, 0.5), y = c(0.5, 0.5))
//' )
//' 
//' state = list(
//'   p = matrix(c(0, 0), 1, 2),
//'   P = list(
//'     matrix(c(1, 1), 1, 2)
//'   )
//' )
//' 
//' attr(state$P[[1]], "i") = 1
//' 
//' # Random centres and ok
//' set.seed(123)
//' centres = matrix(rnorm(66), 33, 2)
//' 
//' ok = as.logical(sample(c(0, 1), 33, replace = TRUE, prob = c(0.3, 0.7)))
//' 
//' seesGoalOK_rcpp(1, objects, state, centres, ok)
//' 
//' @export
// [[Rcpp::export]]
LogicalVector seesGoalOK_rcpp(int n, List objects, List state, NumericMatrix centres, LogicalVector ok) {
  if (is_true(any(ok))) {
    // Get pos of current goal
    NumericVector P_n = get_state_P_n(n - 1, state);
    
    for (int i = 0; i < ok.length(); ++i) {
      // If cell is ok check if goal in sight
      if (ok[i]) {
        ok[i] = seesGoal_rcpp(centres(i, _), P_n, objects);
      }
    }
  }
  
  return(ok);
}