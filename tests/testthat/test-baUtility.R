test_filepath = file.path('data', 'trace_i.rda')
test_obj_name = load(test_filepath)

p = attr(get(test_obj_name), 'pMat')

list_obj = list(
  aBA = p[,'aBA'][1],
  bBA = p[,'bBA'][1],
  BA = trace_i[[1]]$BA[[2]])

p <- c(
  'aBA' = unname(p[,'aBA'][1]),
  'bBA' = unname(p[,'bBA'][1])
)



testthat::test_that("baUtility works", {
  ref = m4ma::baUtility(p, trace_i[[1]]$BA[[2]])
  ll_state = m4ma::baUtility_rcpp(aBA = list_obj$aBA,
                                  bBA = list_obj$bBA,
                                  BA = list_obj$BA,
                                  idx_BA = as.numeric(names(trace_i[[1]]$BA[[2]])))
  testthat::expect_equal(ll_state, ref)
})
