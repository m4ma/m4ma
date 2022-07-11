test_filepath = file.path('data', 'trace_i.rda')
test_obj_name = load(test_filepath)

p = attr(get(test_obj_name), 'pMat')

list_obj = list(
  bGA = p[,'bGA'][1],
  aGA = p[,'aGA'][1],
  GA = trace_i[[1]]$GA[[1]])

p <- c(
  'bGA' = unname(p[,'bGA'][1]),
  'aGA' = unname(p[,'aGA'][1])
)


testthat::test_that("gaUtility works", {
  ref = m4ma::gaUtility(p, trace_i[[1]]$GA[[1]])
  ll_state = m4ma::gaUtility_rcpp(bGA = list_obj$bGA,
                                  aGA = list_obj$aGA,
                                  GA = list_obj$GA)
  testthat::expect_equal(ll_state, ref)
})
