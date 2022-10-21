#include <Rcpp.h>
#include "geometry.h"
#include "see.h"
#include "utility_extra.h"
using namespace Rcpp;


// [[Rcpp::export]]
List prepUtility_rcpp(List state, int N, List objects) {
  List GAs(N);
  List IDs(N);
  List BAs(N);
  List FLs(N);
  List WBs(N);

  for (int i = 0; i < N; i++) {
    NumericVector a = state["a"];
    double a_i = a[i];

    NumericMatrix p = state["p"];
    NumericMatrix p1 = p(Range(i, i), _);

    NumericMatrix P = state["P"];
    NumericMatrix P1 = P(Range(i, i), _);

    NumericVector r = state["r"];

    List centres = state["centres"];
    NumericMatrix centres_i = centres[i];

    NumericMatrix p_pred = state["p_pred"];

    NumericVector group = state["group"];

    GAs[i] = List::create(destinationAngle_rcpp(a_i, p1, P1) / 90.0);
    IDs[i] = List::create(
      predClose_rcpp(i, p1, a_i, p, r, centres_i, p_pred, objects)
    );
    BAs[i] = List::create(
      blockedAngle_rcpp(i, state, p_pred, objects)
    );
    FLs[i] = List::create(
      getLeaders_rcpp(i, state, centres_i, objects)
    );
    WBs[i] = List::create(
      getBuddy_rcpp(i, group, a, p_pred, centres_i, objects, false, state)
    );
  }
  
  state["GA"] = GAs;
  state["ID"] = IDs;
  state["BA"] = BAs;
  state["FL"] = FLs;
  state["WB"] = WBs;
  
  return(state);
}



NumericMatrix get_current_goals(List P) {
  NumericMatrix out(P.length(), 2);

  for(int i = 0; i < P.length(); i++) {
    NumericMatrix P_i = P[i];
    int row_idx = P_i.attr("i");
    CharacterVector P_i_col_names = colnames(P_i);
    P_i_col_names.names() = P_i_col_names;
    int x_col_idx = P_i_col_names.offset("x");
    int y_col_idx = P_i_col_names.offset("y");

    out(i, 0) = P_i(row_idx, x_col_idx);
    out(i, 1) = P_i(row_idx, y_col_idx);
  }

  CharacterVector P_names = P.names();
  rownames(out) = P_names;
  colnames(out) = CharacterVector::create("x", "y");

  return(out);
}



NumericVector get_current_dist(List P) {
  NumericVector out(P.length());

  for(int i = 0; i < P.length(); i++) {
    NumericMatrix P_i = P[i];
    int row_idx = P_i.attr("i");
    CharacterVector P_i_col_names = colnames(P_i);
    P_i_col_names.names() = P_i_col_names;
    int dist_col_idx = P_i_col_names.offset("dist");

    out[i] = P_i(row_idx, dist_col_idx);
  }

  out.names() = P.names();

  return(out);
}


List get_centres(int N, NumericMatrix p, NumericVector v, NumericVector a) {
  List centres(N);
  IntegerVector cells = seq_len(33);
  
  for(int i = 0; i < N; i++) {
    NumericVector p_i = p(i, _);
    centres[i] = c_vd_rcpp(cells, p_i, NumericVector::create(v[i]), a[i],
                           get_vels(), get_angles(), 0.5);
  }
  
  return(centres);
}


bool inObject_rcpp(
    NumericVector xy,
    NumericVector xlim,
    NumericVector ylim,
    bool outside=true
) {
  bool ok = (xy[0] > xlim[0]) & (xy[0] < xlim[1]) & (xy[1] > ylim[0]) & 
    (xy[1] < ylim[1]);
  
  return(ok);
}


LogicalVector isBlocked_rcpp(
    NumericMatrix centres,
    NumericVector p_n,
    List objects,
    LogicalVector ok = LogicalVector(33)
) {
  IntegerVector remain = seq(11, 21);
  
  for(int i = 0; i < 11; i++) {
    NumericVector centres_i = centres(i, _);
    ok[i+22] = seesGoal_rcpp(centres_i, p_n, objects);
  }
  remain = remain[ok[seq(22, 32)]];
  
  if (remain.length() > 0) {
    for(int i = 0; i < remain.length(); i++) {
      NumericVector centres_i = centres(remain[i], _);
      ok[remain[i]] = seesGoal_rcpp(centres_i, p_n, objects);
    }
  }
  
  remain = seq_len(11);
  remain = remain[ok[seq(11, 21)]];
  
  if (remain.length() > 0) {
    for(int i = 0; i < remain.length(); i++) {
      NumericVector centres_i = centres(remain[i], _);
      ok[remain[i]] = seesGoal_rcpp(centres_i, p_n, objects);
    }
  }
  
  return(!ok);
}


NumericMatrix okObject_rcpp(
    int n,
    List objects,
    List state,
    NumericMatrix centres
) {
  LogicalVector ok(centres.nrow());
  
  List object_first = objects[1];
  
  NumericVector xlim = object_first["x"];
  NumericVector ylim = object_first["y"];
  
  for(int i = 0; i < centres.nrow(); i++) {
    NumericVector centres_i = centres(i, _);
    
    ok[i] = inObject_rcpp(centres_i, xlim, ylim, false);
  }
  
  NumericMatrix p = state["p"];
  NumericVector p_n = p(n, _);
  
  if (objects.length() > 1) {
    objects.erase(0);
    LogicalVector blocked = isBlocked_rcpp(centres, p_n, objects);
    ok = ok & !blocked;
  }
  
  ok.attr("dim") = Dimension(11, 3);
  NumericMatrix ok_mat = as<NumericMatrix>(ok);
  
  return(ok_mat);
}


List get_oks(int N, List objects, List state, List centres) {
  List out(N);
  
  for(int i = 0; i < N; i++) {
    NumericMatrix centres_i = centres[i];
    out[i] = okObject_rcpp(i, objects, state, centres_i);
  }
  
  return(out);
}


NumericMatrix predictPred_rcpp(
    NumericMatrix p,
    NumericVector v,
    NumericVector a,
    Nullable<IntegerVector> cell,
    bool stayStop = true
) {
  NumericMatrix out = clone(p);
  
  NumericVector v_scaled = scaleVel_rcpp(v, 0.5);
  NumericMatrix sine_cosine = aTOd_rcpp(a);
  
  for(int i = 0; i < p.nrow(); i++) {
    out(i, _) = p(i, _) + v_scaled[i] * sine_cosine(i, _);
  }
  
  if (stayStop && cell != R_NilValue) {
    IntegerVector cell_int = as<IntegerVector>(cell);
    
    for(int i = 0; i < p.nrow(); i++) {
      if (cell_int[i] < 1) {
        out(i, _) = p(i, _);
      }
    }
  }
  
  return(out);
}


// [[Rcpp::export]]
List prepSimTrace_rcpp(List sim_trace, Nullable<NumericVector> constant = R_NilValue) {
  List space = sim_trace.attr("space");
  List objects = space["objects"];
  
  List state_first = sim_trace[0];
  NumericMatrix pMat_first = state_first["pMat"];
  
  NumericMatrix pMat_t(pMat_first.ncol(), 0);
  CharacterVector pMat_row_names;
  
  List trace(sim_trace.length());
  
  for(int i = 0; i < sim_trace.length(); i++) {
    List state_i = sim_trace[i];
    
    IntegerVector cell = state_i["cell"];
    cell[cell == -1] = 0;
    state_i["cell"] = cell;

    NumericMatrix p = state_i["p"];
    
    int N = p.nrow();

    List P = state_i["P"];
  
    NumericVector v = state_i["v"];
    NumericVector a = state_i["a"];
    
    NumericMatrix current_goals = get_current_goals(P);
    
    NumericVector current_dist = get_current_dist(P);
    
    List centres = get_centres(N, p, v, a);
    
    List oks = get_oks(N, objects, state_i, centres);
    
    NumericMatrix p_pred = predictPred_rcpp(p, v, a, cell);
    
    state_i["P"] = current_goals;
    state_i["d"] = current_dist;
    state_i["centres"] = centres;
    state_i["p_pred"] = p_pred;
    state_i["ok"] = oks;
    
    state_i = prepUtility_rcpp(state_i, N, objects);
    
    state_i["bads"] = R_NilValue;
    state_i["a"] = R_NilValue;
    state_i["P"] = R_NilValue;
    state_i["centres"] = R_NilValue;
    state_i["p_pred"] = R_NilValue;
    state_i["r"] = R_NilValue;
    
    trace[i] = state_i;
    
    // NOTE: This is a really hacky solution; should look for better one!
    NumericMatrix pMat_i = state_i["pMat"];
    CharacterVector pMat_i_row_names = rownames(pMat_i);
    
    for(int j = 0; j < pMat_i.nrow(); j++) {
      String name_j = pMat_i_row_names[j];
      if (!pMat_row_names.containsElementNamed(pMat_i_row_names[j])) {
        NumericMatrix pMat_i_j = pMat_i(Range(j, j), _);
        pMat_t = cbind(pMat_t, transpose(pMat_i_j));
        pMat_row_names.push_back(name_j);
        pMat_row_names.names() = pMat_row_names;
      }
    }
  }
  
  NumericMatrix pMat = transpose(pMat_t);
  
  trace.attr("alpha") = sim_trace.attr("alpha");
  trace.attr("nests") = sim_trace.attr("nests");
  
  if (constant != R_NilValue) {
    NumericVector constant(constant);
    
    int n_const = constant.length();
    
    NumericVector const_vec = rep_each(constant, pMat.nrow());
    const_vec.attr("dim") = Dimension(pMat.nrow(), n_const);
    
    NumericMatrix const_mat = as<NumericMatrix>(const_vec);
    
    pMat = cbind(pMat, const_mat);
    CharacterVector constant_names = constant.names();
    colnames(pMat) = constant_names;
  }
  
  rownames(pMat) = pMat_row_names;
  colnames(pMat) = colnames(pMat_first);
  
  trace.attr("pMat") = pMat;
  trace.attr("subject_names") = pMat_row_names;
  trace.attr("param_names") = colnames(pMat_first);
  
  return(trace);
}
