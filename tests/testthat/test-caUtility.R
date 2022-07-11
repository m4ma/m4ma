test_filepath = file.path('data', 'trace_i.rda')
test_obj_name = load(test_filepath)

p = attr(get(test_obj_name), 'pMat')

list_obj = list(
  aCA = unname(p[,'aCA'][1]),
  bCA = unname(p[,'bCA'][1]),
  bCAlr = unname(p[,'bCAlr'][1]))

p <- c(
  'aCA' = unname(p[,'aCA'][1]),
  'bCA' = unname(p[,'bCA'][1]),
  'bCAlr' = unname(p[,'bCAlr'][1])
)



testthat::test_that("caUtility works", {
  ref = unname(m4ma::caUtility(p))
  ll_state = m4ma::caUtility_rcpp(aCA = list_obj$aCA,
                                  bCA = list_obj$bCA,
                                  bCAlr = list_obj$bCAlr)
  testthat::expect_equal(ll_state, ref)
})
