// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// object2lines_rcpp
List object2lines_rcpp(List o);
RcppExport SEXP _m4ma_object2lines_rcpp(SEXP oSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type o(oSEXP);
    rcpp_result_gen = Rcpp::wrap(object2lines_rcpp(o));
    return rcpp_result_gen;
END_RCPP
}
// bodyObjectOverlap_rcpp
LogicalVector bodyObjectOverlap_rcpp(List oL, double r, NumericMatrix okCentres);
RcppExport SEXP _m4ma_bodyObjectOverlap_rcpp(SEXP oLSEXP, SEXP rSEXP, SEXP okCentresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type oL(oLSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type okCentres(okCentresSEXP);
    rcpp_result_gen = Rcpp::wrap(bodyObjectOverlap_rcpp(oL, r, okCentres));
    return rcpp_result_gen;
END_RCPP
}
// bodyObjectOK_rcpp
Nullable<LogicalMatrix> bodyObjectOK_rcpp(double r, NumericMatrix centres, List objects, LogicalVector ok);
RcppExport SEXP _m4ma_bodyObjectOK_rcpp(SEXP rSEXP, SEXP centresSEXP, SEXP objectsSEXP, SEXP okSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type centres(centresSEXP);
    Rcpp::traits::input_parameter< List >::type objects(objectsSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type ok(okSEXP);
    rcpp_result_gen = Rcpp::wrap(bodyObjectOK_rcpp(r, centres, objects, ok));
    return rcpp_result_gen;
END_RCPP
}
// free_cells_rcpp
LogicalMatrix free_cells_rcpp(S4 agent, S4 background, NumericMatrix centers);
RcppExport SEXP _m4ma_free_cells_rcpp(SEXP agentSEXP, SEXP backgroundSEXP, SEXP centersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type agent(agentSEXP);
    Rcpp::traits::input_parameter< S4 >::type background(backgroundSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type centers(centersSEXP);
    rcpp_result_gen = Rcpp::wrap(free_cells_rcpp(agent, background, centers));
    return rcpp_result_gen;
END_RCPP
}
// dist_rcpp
NumericVector dist_rcpp(NumericMatrix p1, NumericMatrix p2);
RcppExport SEXP _m4ma_dist_rcpp(SEXP p1SEXP, SEXP p2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p2(p2SEXP);
    rcpp_result_gen = Rcpp::wrap(dist_rcpp(p1, p2));
    return rcpp_result_gen;
END_RCPP
}
// dist1_rcpp
NumericVector dist1_rcpp(NumericVector p1, NumericMatrix p2);
RcppExport SEXP _m4ma_dist1_rcpp(SEXP p1SEXP, SEXP p2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p2(p2SEXP);
    rcpp_result_gen = Rcpp::wrap(dist1_rcpp(p1, p2));
    return rcpp_result_gen;
END_RCPP
}
// angle2s_rcpp
NumericVector angle2s_rcpp(NumericMatrix p1, NumericMatrix p2);
RcppExport SEXP _m4ma_angle2s_rcpp(SEXP p1SEXP, SEXP p2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p2(p2SEXP);
    rcpp_result_gen = Rcpp::wrap(angle2s_rcpp(p1, p2));
    return rcpp_result_gen;
END_RCPP
}
// angle2_rcpp
NumericVector angle2_rcpp(NumericMatrix p1, NumericMatrix p2);
RcppExport SEXP _m4ma_angle2_rcpp(SEXP p1SEXP, SEXP p2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p2(p2SEXP);
    rcpp_result_gen = Rcpp::wrap(angle2_rcpp(p1, p2));
    return rcpp_result_gen;
END_RCPP
}
// aTOd_rcpp
NumericMatrix aTOd_rcpp(NumericVector a);
RcppExport SEXP _m4ma_aTOd_rcpp(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(aTOd_rcpp(a));
    return rcpp_result_gen;
END_RCPP
}
// Iangle_rcpp
NumericVector Iangle_rcpp(NumericMatrix p1, double a1, NumericMatrix p2);
RcppExport SEXP _m4ma_Iangle_rcpp(SEXP p1SEXP, SEXP a1SEXP, SEXP p2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p2(p2SEXP);
    rcpp_result_gen = Rcpp::wrap(Iangle_rcpp(p1, a1, p2));
    return rcpp_result_gen;
END_RCPP
}
// Dn_rcpp
NumericVector Dn_rcpp(NumericMatrix p_n, NumericMatrix P_n);
RcppExport SEXP _m4ma_Dn_rcpp(SEXP p_nSEXP, SEXP P_nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p_n(p_nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P_n(P_nSEXP);
    rcpp_result_gen = Rcpp::wrap(Dn_rcpp(p_n, P_n));
    return rcpp_result_gen;
END_RCPP
}
// minAngle_rcpp
NumericVector minAngle_rcpp(double a1_double, NumericVector a2);
RcppExport SEXP _m4ma_minAngle_rcpp(SEXP a1_doubleSEXP, SEXP a2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a1_double(a1_doubleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a2(a2SEXP);
    rcpp_result_gen = Rcpp::wrap(minAngle_rcpp(a1_double, a2));
    return rcpp_result_gen;
END_RCPP
}
// headingAngle_rcpp
NumericMatrix headingAngle_rcpp(NumericVector a2, double a1);
RcppExport SEXP _m4ma_headingAngle_rcpp(SEXP a2SEXP, SEXP a1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    rcpp_result_gen = Rcpp::wrap(headingAngle_rcpp(a2, a1));
    return rcpp_result_gen;
END_RCPP
}
// scaleVel_rcpp
NumericVector scaleVel_rcpp(NumericVector v, double tStep);
RcppExport SEXP _m4ma_scaleVel_rcpp(SEXP vSEXP, SEXP tStepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type tStep(tStepSEXP);
    rcpp_result_gen = Rcpp::wrap(scaleVel_rcpp(v, tStep));
    return rcpp_result_gen;
END_RCPP
}
// c_vd_rcpp
NumericMatrix c_vd_rcpp(IntegerVector cells, NumericVector p1, double v1, double a1, NumericMatrix vels, NumericMatrix angles, double tStep);
RcppExport SEXP _m4ma_c_vd_rcpp(SEXP cellsSEXP, SEXP p1SEXP, SEXP v1SEXP, SEXP a1SEXP, SEXP velsSEXP, SEXP anglesSEXP, SEXP tStepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type cells(cellsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< double >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type vels(velsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type angles(anglesSEXP);
    Rcpp::traits::input_parameter< double >::type tStep(tStepSEXP);
    rcpp_result_gen = Rcpp::wrap(c_vd_rcpp(cells, p1, v1, a1, vels, angles, tStep));
    return rcpp_result_gen;
END_RCPP
}
// coneNum_rcpp
NumericVector coneNum_rcpp(NumericVector k);
RcppExport SEXP _m4ma_coneNum_rcpp(SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(coneNum_rcpp(k));
    return rcpp_result_gen;
END_RCPP
}
// ringNum_rcpp
NumericVector ringNum_rcpp(NumericVector k);
RcppExport SEXP _m4ma_ringNum_rcpp(SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(ringNum_rcpp(k));
    return rcpp_result_gen;
END_RCPP
}
// like_observation
double like_observation(List obs, NumericVector p, int n, List nests, List alpha, NumericMatrix cell_nest, double min_like);
RcppExport SEXP _m4ma_like_observation(SEXP obsSEXP, SEXP pSEXP, SEXP nSEXP, SEXP nestsSEXP, SEXP alphaSEXP, SEXP cell_nestSEXP, SEXP min_likeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< List >::type nests(nestsSEXP);
    Rcpp::traits::input_parameter< List >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cell_nest(cell_nestSEXP);
    Rcpp::traits::input_parameter< double >::type min_like(min_likeSEXP);
    rcpp_result_gen = Rcpp::wrap(like_observation(obs, p, n, nests, alpha, cell_nest, min_like));
    return rcpp_result_gen;
END_RCPP
}
// like_state
NumericVector like_state(List state, int ti, NumericMatrix p, List nests, List alpha, NumericMatrix cell_nest, std::string elements, double min_like);
RcppExport SEXP _m4ma_like_state(SEXP stateSEXP, SEXP tiSEXP, SEXP pSEXP, SEXP nestsSEXP, SEXP alphaSEXP, SEXP cell_nestSEXP, SEXP elementsSEXP, SEXP min_likeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< int >::type ti(tiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP);
    Rcpp::traits::input_parameter< List >::type nests(nestsSEXP);
    Rcpp::traits::input_parameter< List >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cell_nest(cell_nestSEXP);
    Rcpp::traits::input_parameter< std::string >::type elements(elementsSEXP);
    Rcpp::traits::input_parameter< double >::type min_like(min_likeSEXP);
    rcpp_result_gen = Rcpp::wrap(like_state(state, ti, p, nests, alpha, cell_nest, elements, min_like));
    return rcpp_result_gen;
END_RCPP
}
// msumlogLike
double msumlogLike(NumericMatrix p, List trace_rcpp, List nests, List alpha, NumericMatrix cell_nest, double min_like, double mult);
RcppExport SEXP _m4ma_msumlogLike(SEXP pSEXP, SEXP trace_rcppSEXP, SEXP nestsSEXP, SEXP alphaSEXP, SEXP cell_nestSEXP, SEXP min_likeSEXP, SEXP multSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p(pSEXP);
    Rcpp::traits::input_parameter< List >::type trace_rcpp(trace_rcppSEXP);
    Rcpp::traits::input_parameter< List >::type nests(nestsSEXP);
    Rcpp::traits::input_parameter< List >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cell_nest(cell_nestSEXP);
    Rcpp::traits::input_parameter< double >::type min_like(min_likeSEXP);
    Rcpp::traits::input_parameter< double >::type mult(multSEXP);
    rcpp_result_gen = Rcpp::wrap(msumlogLike(p, trace_rcpp, nests, alpha, cell_nest, min_like, mult));
    return rcpp_result_gen;
END_RCPP
}
// pcnl_rcpp
double pcnl_rcpp(NumericVector cell, NumericVector utility, NumericVector mum, List nests, List alpha, double mu);
RcppExport SEXP _m4ma_pcnl_rcpp(SEXP cellSEXP, SEXP utilitySEXP, SEXP mumSEXP, SEXP nestsSEXP, SEXP alphaSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type cell(cellSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type utility(utilitySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mum(mumSEXP);
    Rcpp::traits::input_parameter< List >::type nests(nestsSEXP);
    Rcpp::traits::input_parameter< List >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(pcnl_rcpp(cell, utility, mum, nests, alpha, mu));
    return rcpp_result_gen;
END_RCPP
}
// pmnl_rcpp
double pmnl_rcpp(int cell, NumericVector utility, LogicalMatrix ok);
RcppExport SEXP _m4ma_pmnl_rcpp(SEXP cellSEXP, SEXP utilitySEXP, SEXP okSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type cell(cellSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type utility(utilitySEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type ok(okSEXP);
    rcpp_result_gen = Rcpp::wrap(pmnl_rcpp(cell, utility, ok));
    return rcpp_result_gen;
END_RCPP
}
// line_line_intersection_rcpp
NumericVector line_line_intersection_rcpp(NumericVector P1, NumericVector P2, NumericVector P3, NumericVector P4, bool interior_only);
RcppExport SEXP _m4ma_line_line_intersection_rcpp(SEXP P1SEXP, SEXP P2SEXP, SEXP P3SEXP, SEXP P4SEXP, SEXP interior_onlySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P2(P2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P3(P3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P4(P4SEXP);
    Rcpp::traits::input_parameter< bool >::type interior_only(interior_onlySEXP);
    rcpp_result_gen = Rcpp::wrap(line_line_intersection_rcpp(P1, P2, P3, P4, interior_only));
    return rcpp_result_gen;
END_RCPP
}
// seesGoal_rcpp
bool seesGoal_rcpp(NumericVector p_n, NumericVector P_n, List objects);
RcppExport SEXP _m4ma_seesGoal_rcpp(SEXP p_nSEXP, SEXP P_nSEXP, SEXP objectsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p_n(p_nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P_n(P_nSEXP);
    Rcpp::traits::input_parameter< List >::type objects(objectsSEXP);
    rcpp_result_gen = Rcpp::wrap(seesGoal_rcpp(p_n, P_n, objects));
    return rcpp_result_gen;
END_RCPP
}
// seesCurrentGoal_rcpp
bool seesCurrentGoal_rcpp(int n, List state, List objects, int offset);
RcppExport SEXP _m4ma_seesCurrentGoal_rcpp(SEXP nSEXP, SEXP stateSEXP, SEXP objectsSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< List >::type objects(objectsSEXP);
    Rcpp::traits::input_parameter< int >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(seesCurrentGoal_rcpp(n, state, objects, offset));
    return rcpp_result_gen;
END_RCPP
}
// seesMany_rcpp
LogicalVector seesMany_rcpp(NumericVector p1, NumericMatrix ps, List objects);
RcppExport SEXP _m4ma_seesMany_rcpp(SEXP p1SEXP, SEXP psSEXP, SEXP objectsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ps(psSEXP);
    Rcpp::traits::input_parameter< List >::type objects(objectsSEXP);
    rcpp_result_gen = Rcpp::wrap(seesMany_rcpp(p1, ps, objects));
    return rcpp_result_gen;
END_RCPP
}
// seesGoalOK_rcpp
LogicalVector seesGoalOK_rcpp(int n, List objects, List state, NumericMatrix centres, LogicalVector ok);
RcppExport SEXP _m4ma_seesGoalOK_rcpp(SEXP nSEXP, SEXP objectsSEXP, SEXP stateSEXP, SEXP centresSEXP, SEXP okSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< List >::type objects(objectsSEXP);
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type centres(centresSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type ok(okSEXP);
    rcpp_result_gen = Rcpp::wrap(seesGoalOK_rcpp(n, objects, state, centres, ok));
    return rcpp_result_gen;
END_RCPP
}
// baUtility_rcpp
NumericVector baUtility_rcpp(double aBA, double bBA, NumericVector BA, IntegerVector idx_BA);
RcppExport SEXP _m4ma_baUtility_rcpp(SEXP aBASEXP, SEXP bBASEXP, SEXP BASEXP, SEXP idx_BASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type aBA(aBASEXP);
    Rcpp::traits::input_parameter< double >::type bBA(bBASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type BA(BASEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idx_BA(idx_BASEXP);
    rcpp_result_gen = Rcpp::wrap(baUtility_rcpp(aBA, bBA, BA, idx_BA));
    return rcpp_result_gen;
END_RCPP
}
// caUtility_rcpp
NumericVector caUtility_rcpp(double aCA, double bCA, double bCAlr);
RcppExport SEXP _m4ma_caUtility_rcpp(SEXP aCASEXP, SEXP bCASEXP, SEXP bCAlrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type aCA(aCASEXP);
    Rcpp::traits::input_parameter< double >::type bCA(bCASEXP);
    Rcpp::traits::input_parameter< double >::type bCAlr(bCAlrSEXP);
    rcpp_result_gen = Rcpp::wrap(caUtility_rcpp(aCA, bCA, bCAlr));
    return rcpp_result_gen;
END_RCPP
}
// flUtility_rcpp
NumericVector flUtility_rcpp(double aFL, double bFL, double dFL, NumericMatrix leaders, NumericMatrix dists);
RcppExport SEXP _m4ma_flUtility_rcpp(SEXP aFLSEXP, SEXP bFLSEXP, SEXP dFLSEXP, SEXP leadersSEXP, SEXP distsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type aFL(aFLSEXP);
    Rcpp::traits::input_parameter< double >::type bFL(bFLSEXP);
    Rcpp::traits::input_parameter< double >::type dFL(dFLSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type leaders(leadersSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dists(distsSEXP);
    rcpp_result_gen = Rcpp::wrap(flUtility_rcpp(aFL, bFL, dFL, leaders, dists));
    return rcpp_result_gen;
END_RCPP
}
// gaUtility_rcpp
NumericVector gaUtility_rcpp(double bGA, double aGA, NumericVector GA);
RcppExport SEXP _m4ma_gaUtility_rcpp(SEXP bGASEXP, SEXP aGASEXP, SEXP GASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type bGA(bGASEXP);
    Rcpp::traits::input_parameter< double >::type aGA(aGASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type GA(GASEXP);
    rcpp_result_gen = Rcpp::wrap(gaUtility_rcpp(bGA, aGA, GA));
    return rcpp_result_gen;
END_RCPP
}
// idUtility_rcpp
NumericVector idUtility_rcpp(double bID, double dID, double aID, double n, const LogicalMatrix ok, IntegerVector group, Nullable<NumericMatrix> ID_);
RcppExport SEXP _m4ma_idUtility_rcpp(SEXP bIDSEXP, SEXP dIDSEXP, SEXP aIDSEXP, SEXP nSEXP, SEXP okSEXP, SEXP groupSEXP, SEXP ID_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type bID(bIDSEXP);
    Rcpp::traits::input_parameter< double >::type dID(dIDSEXP);
    Rcpp::traits::input_parameter< double >::type aID(aIDSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type ok(okSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericMatrix> >::type ID_(ID_SEXP);
    rcpp_result_gen = Rcpp::wrap(idUtility_rcpp(bID, dID, aID, n, ok, group, ID_));
    return rcpp_result_gen;
END_RCPP
}
// psUtility_rcpp
NumericVector psUtility_rcpp(double aPS, double bPS, double sPref, double sSlow, double v, double d);
RcppExport SEXP _m4ma_psUtility_rcpp(SEXP aPSSEXP, SEXP bPSSEXP, SEXP sPrefSEXP, SEXP sSlowSEXP, SEXP vSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type aPS(aPSSEXP);
    Rcpp::traits::input_parameter< double >::type bPS(bPSSEXP);
    Rcpp::traits::input_parameter< double >::type sPref(sPrefSEXP);
    Rcpp::traits::input_parameter< double >::type sSlow(sSlowSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(psUtility_rcpp(aPS, bPS, sPref, sSlow, v, d));
    return rcpp_result_gen;
END_RCPP
}
// wbUtility_rcpp
NumericVector wbUtility_rcpp(double aWB, double bWB, NumericMatrix buddies, NumericMatrix dists);
RcppExport SEXP _m4ma_wbUtility_rcpp(SEXP aWBSEXP, SEXP bWBSEXP, SEXP buddiesSEXP, SEXP distsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type aWB(aWBSEXP);
    Rcpp::traits::input_parameter< double >::type bWB(bWBSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type buddies(buddiesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dists(distsSEXP);
    rcpp_result_gen = Rcpp::wrap(wbUtility_rcpp(aWB, bWB, buddies, dists));
    return rcpp_result_gen;
END_RCPP
}
// utility
NumericVector utility(NumericVector p, int n, double v, double d, Nullable<NumericVector> ba_, NumericVector ga, Nullable<NumericMatrix> id_, Nullable<List> fl_, Nullable<List> wb_, LogicalMatrix ok, IntegerVector group);
RcppExport SEXP _m4ma_utility(SEXP pSEXP, SEXP nSEXP, SEXP vSEXP, SEXP dSEXP, SEXP ba_SEXP, SEXP gaSEXP, SEXP id_SEXP, SEXP fl_SEXP, SEXP wb_SEXP, SEXP okSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type ba_(ba_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ga(gaSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericMatrix> >::type id_(id_SEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type fl_(fl_SEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type wb_(wb_SEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type ok(okSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(utility(p, n, v, d, ba_, ga, id_, fl_, wb_, ok, group));
    return rcpp_result_gen;
END_RCPP
}
// destinationAngle_rcpp
NumericVector destinationAngle_rcpp(double a, NumericMatrix p1, NumericMatrix P1);
RcppExport SEXP _m4ma_destinationAngle_rcpp(SEXP aSEXP, SEXP p1SEXP, SEXP P1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P1(P1SEXP);
    rcpp_result_gen = Rcpp::wrap(destinationAngle_rcpp(a, p1, P1));
    return rcpp_result_gen;
END_RCPP
}
// predClose_rcpp
Nullable<NumericMatrix> predClose_rcpp(int n, NumericMatrix p1, double a1, NumericMatrix p2, NumericVector r, NumericMatrix centres, NumericMatrix p_pred, List objects);
RcppExport SEXP _m4ma_predClose_rcpp(SEXP nSEXP, SEXP p1SEXP, SEXP a1SEXP, SEXP p2SEXP, SEXP rSEXP, SEXP centresSEXP, SEXP p_predSEXP, SEXP objectsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type centres(centresSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p_pred(p_predSEXP);
    Rcpp::traits::input_parameter< List >::type objects(objectsSEXP);
    rcpp_result_gen = Rcpp::wrap(predClose_rcpp(n, p1, a1, p2, r, centres, p_pred, objects));
    return rcpp_result_gen;
END_RCPP
}
// eObjects_rcpp
List eObjects_rcpp(NumericMatrix p1, NumericMatrix p2, NumericVector r);
RcppExport SEXP _m4ma_eObjects_rcpp(SEXP p1SEXP, SEXP p2SEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(eObjects_rcpp(p1, p2, r));
    return rcpp_result_gen;
END_RCPP
}
// iCones_rcpp
Nullable<NumericVector> iCones_rcpp(NumericMatrix p1, double a, NumericMatrix p2, NumericVector r, List objects);
RcppExport SEXP _m4ma_iCones_rcpp(SEXP p1SEXP, SEXP aSEXP, SEXP p2SEXP, SEXP rSEXP, SEXP objectsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< List >::type objects(objectsSEXP);
    rcpp_result_gen = Rcpp::wrap(iCones_rcpp(p1, a, p2, r, objects));
    return rcpp_result_gen;
END_RCPP
}
// iCones2Cells_rcpp
NumericVector iCones2Cells_rcpp(NumericVector iC, double v, double tStep);
RcppExport SEXP _m4ma_iCones2Cells_rcpp(SEXP iCSEXP, SEXP vSEXP, SEXP tStepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type iC(iCSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type tStep(tStepSEXP);
    rcpp_result_gen = Rcpp::wrap(iCones2Cells_rcpp(iC, v, tStep));
    return rcpp_result_gen;
END_RCPP
}
// blockedAngle_rcpp
NumericVector blockedAngle_rcpp(NumericMatrix p1, double a1, double v1, NumericMatrix p2, NumericVector r, List objects);
RcppExport SEXP _m4ma_blockedAngle_rcpp(SEXP p1SEXP, SEXP a1SEXP, SEXP v1SEXP, SEXP p2SEXP, SEXP rSEXP, SEXP objectsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< List >::type objects(objectsSEXP);
    rcpp_result_gen = Rcpp::wrap(blockedAngle_rcpp(p1, a1, v1, p2, r, objects));
    return rcpp_result_gen;
END_RCPP
}
// getLeaders_rcpp
Nullable<List> getLeaders_rcpp(int n, NumericMatrix p_mat, NumericVector a, NumericVector v, NumericMatrix P1, NumericVector group, NumericMatrix centres, List objects, bool onlyGroup, bool preferGroup, bool pickBest);
RcppExport SEXP _m4ma_getLeaders_rcpp(SEXP nSEXP, SEXP p_matSEXP, SEXP aSEXP, SEXP vSEXP, SEXP P1SEXP, SEXP groupSEXP, SEXP centresSEXP, SEXP objectsSEXP, SEXP onlyGroupSEXP, SEXP preferGroupSEXP, SEXP pickBestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p_mat(p_matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type centres(centresSEXP);
    Rcpp::traits::input_parameter< List >::type objects(objectsSEXP);
    Rcpp::traits::input_parameter< bool >::type onlyGroup(onlyGroupSEXP);
    Rcpp::traits::input_parameter< bool >::type preferGroup(preferGroupSEXP);
    Rcpp::traits::input_parameter< bool >::type pickBest(pickBestSEXP);
    rcpp_result_gen = Rcpp::wrap(getLeaders_rcpp(n, p_mat, a, v, P1, group, centres, objects, onlyGroup, preferGroup, pickBest));
    return rcpp_result_gen;
END_RCPP
}
// getBuddy_rcpp
Nullable<List> getBuddy_rcpp(int n, NumericMatrix p_mat, NumericVector v, NumericVector group, NumericVector a, NumericMatrix p_pred, NumericMatrix centres, List objects, bool pickBest);
RcppExport SEXP _m4ma_getBuddy_rcpp(SEXP nSEXP, SEXP p_matSEXP, SEXP vSEXP, SEXP groupSEXP, SEXP aSEXP, SEXP p_predSEXP, SEXP centresSEXP, SEXP objectsSEXP, SEXP pickBestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p_mat(p_matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p_pred(p_predSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type centres(centresSEXP);
    Rcpp::traits::input_parameter< List >::type objects(objectsSEXP);
    Rcpp::traits::input_parameter< bool >::type pickBest(pickBestSEXP);
    rcpp_result_gen = Rcpp::wrap(getBuddy_rcpp(n, p_mat, v, group, a, p_pred, centres, objects, pickBest));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_m4ma_object2lines_rcpp", (DL_FUNC) &_m4ma_object2lines_rcpp, 1},
    {"_m4ma_bodyObjectOverlap_rcpp", (DL_FUNC) &_m4ma_bodyObjectOverlap_rcpp, 3},
    {"_m4ma_bodyObjectOK_rcpp", (DL_FUNC) &_m4ma_bodyObjectOK_rcpp, 4},
    {"_m4ma_free_cells_rcpp", (DL_FUNC) &_m4ma_free_cells_rcpp, 3},
    {"_m4ma_dist_rcpp", (DL_FUNC) &_m4ma_dist_rcpp, 2},
    {"_m4ma_dist1_rcpp", (DL_FUNC) &_m4ma_dist1_rcpp, 2},
    {"_m4ma_angle2s_rcpp", (DL_FUNC) &_m4ma_angle2s_rcpp, 2},
    {"_m4ma_angle2_rcpp", (DL_FUNC) &_m4ma_angle2_rcpp, 2},
    {"_m4ma_aTOd_rcpp", (DL_FUNC) &_m4ma_aTOd_rcpp, 1},
    {"_m4ma_Iangle_rcpp", (DL_FUNC) &_m4ma_Iangle_rcpp, 3},
    {"_m4ma_Dn_rcpp", (DL_FUNC) &_m4ma_Dn_rcpp, 2},
    {"_m4ma_minAngle_rcpp", (DL_FUNC) &_m4ma_minAngle_rcpp, 2},
    {"_m4ma_headingAngle_rcpp", (DL_FUNC) &_m4ma_headingAngle_rcpp, 2},
    {"_m4ma_scaleVel_rcpp", (DL_FUNC) &_m4ma_scaleVel_rcpp, 2},
    {"_m4ma_c_vd_rcpp", (DL_FUNC) &_m4ma_c_vd_rcpp, 7},
    {"_m4ma_coneNum_rcpp", (DL_FUNC) &_m4ma_coneNum_rcpp, 1},
    {"_m4ma_ringNum_rcpp", (DL_FUNC) &_m4ma_ringNum_rcpp, 1},
    {"_m4ma_like_observation", (DL_FUNC) &_m4ma_like_observation, 7},
    {"_m4ma_like_state", (DL_FUNC) &_m4ma_like_state, 8},
    {"_m4ma_msumlogLike", (DL_FUNC) &_m4ma_msumlogLike, 7},
    {"_m4ma_pcnl_rcpp", (DL_FUNC) &_m4ma_pcnl_rcpp, 6},
    {"_m4ma_pmnl_rcpp", (DL_FUNC) &_m4ma_pmnl_rcpp, 3},
    {"_m4ma_line_line_intersection_rcpp", (DL_FUNC) &_m4ma_line_line_intersection_rcpp, 5},
    {"_m4ma_seesGoal_rcpp", (DL_FUNC) &_m4ma_seesGoal_rcpp, 3},
    {"_m4ma_seesCurrentGoal_rcpp", (DL_FUNC) &_m4ma_seesCurrentGoal_rcpp, 4},
    {"_m4ma_seesMany_rcpp", (DL_FUNC) &_m4ma_seesMany_rcpp, 3},
    {"_m4ma_seesGoalOK_rcpp", (DL_FUNC) &_m4ma_seesGoalOK_rcpp, 5},
    {"_m4ma_baUtility_rcpp", (DL_FUNC) &_m4ma_baUtility_rcpp, 4},
    {"_m4ma_caUtility_rcpp", (DL_FUNC) &_m4ma_caUtility_rcpp, 3},
    {"_m4ma_flUtility_rcpp", (DL_FUNC) &_m4ma_flUtility_rcpp, 5},
    {"_m4ma_gaUtility_rcpp", (DL_FUNC) &_m4ma_gaUtility_rcpp, 3},
    {"_m4ma_idUtility_rcpp", (DL_FUNC) &_m4ma_idUtility_rcpp, 7},
    {"_m4ma_psUtility_rcpp", (DL_FUNC) &_m4ma_psUtility_rcpp, 6},
    {"_m4ma_wbUtility_rcpp", (DL_FUNC) &_m4ma_wbUtility_rcpp, 4},
    {"_m4ma_utility", (DL_FUNC) &_m4ma_utility, 11},
    {"_m4ma_destinationAngle_rcpp", (DL_FUNC) &_m4ma_destinationAngle_rcpp, 3},
    {"_m4ma_predClose_rcpp", (DL_FUNC) &_m4ma_predClose_rcpp, 8},
    {"_m4ma_eObjects_rcpp", (DL_FUNC) &_m4ma_eObjects_rcpp, 3},
    {"_m4ma_iCones_rcpp", (DL_FUNC) &_m4ma_iCones_rcpp, 5},
    {"_m4ma_iCones2Cells_rcpp", (DL_FUNC) &_m4ma_iCones2Cells_rcpp, 3},
    {"_m4ma_blockedAngle_rcpp", (DL_FUNC) &_m4ma_blockedAngle_rcpp, 6},
    {"_m4ma_getLeaders_rcpp", (DL_FUNC) &_m4ma_getLeaders_rcpp, 11},
    {"_m4ma_getBuddy_rcpp", (DL_FUNC) &_m4ma_getBuddy_rcpp, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_m4ma(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
