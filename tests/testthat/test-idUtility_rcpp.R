test_filepath = file.path('data', 'trace_i.rda')
test_obj_name = load(test_filepath)

p = attr(get(test_obj_name), 'pMat')

list_obj = list(
ok = trace_i[[1]]$ok[[1]],
bID = p[,'bID'][1],
dID = p[,'dID'][1],
aID = p[,'aID'][1],
group = trace_i[[1]]$group,
ID = trace_i[[1]]$ID[2][[1]])

p <- c(
  'bID' = unname(p[,'bID'][1]),
  'dID' = unname(p[,'dID'][1]),
  'aID' = unname(p[,'aID'][1])
)


testthat::test_that("idUtility outputs a vector of 0s when ID = NULL", {
  ref = rep(0, 33)
  ll_state = m4ma::idUtility_rcpp(bID = list_obj$bID,
                                  dID = list_obj$dID,
                                  aID = list_obj$aID,
                                  n = 1,
                                  group = list_obj$group,
                                  ok = list_obj$ok,
                                  ID = NULL)
  testthat::expect_equal(ll_state, ref)
})


testthat::test_that("idUtility outputs a vector of negative doubles when ID != NULL", {
  ref = m4ma::idUtility(p, n = 1, ID = trace_i[[1]]$ID[2][[1]], trace_i[[1]]$ok[[1]], trace_i[[1]]$group)
  ll_state = m4ma::idUtility_rcpp(bID = list_obj$bID,
                                  dID = list_obj$dID,
                                  aID = list_obj$aID,
                                  n = 1,
                                  group = list_obj$group,
                                  ok = list_obj$ok,
                                  ID = list_obj$ID)
  testthat::expect_equal(ll_state, ref)
})


