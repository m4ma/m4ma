
test_filepath = file.path('data', 'trace_i.rda')

test_obj_name = load(test_filepath)

p = attr(get(test_obj_name), 'pMat')

state = get(test_obj_name)[[1]]

testthat::test_that("Utility computation works", {
  idx = 1
  ref = -100632669
  utility = m4ma::utility_rcpp(
    p[idx, ], idx, state$v[idx], state$d[idx], state$GA[[idx]], state$BA[[idx]],
    state$ID[[idx]],state$FL[[idx]], state$WB[[idx]], state$ok[[idx]],
    state$group
  )
  testthat::expect_equal(sum(utility), ref, tolerance = 1e-5)
})