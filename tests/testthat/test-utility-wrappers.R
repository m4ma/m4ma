# Test rcpp utility functions against R counterparts on test trace

test_filepath = file.path('data', 'trace_i.rda')

test_obj_name = load(test_filepath)

p = attr(get(test_obj_name), 'pMat')


testthat::test_that("Blocked-angle utility wrapper works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      if (!is.null(state$FL[[i]])) {
        ref = m4ma::baUtility(p[i, ], state$BA[[i]], use = 'r')
        ba_utility = ref = m4ma::baUtility(p[i, ], state$BA[[i]], use = 'cpp')
        testthat::expect_equal(ba_utility, ref)
      }
    }
  })
})


testthat::test_that("Current-angle utility wrapper works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = unname(m4ma::caUtility(p[i, ], use = 'r'))
      ca_utility = unname(m4ma::caUtility(p[i, ], use = 'cpp'))
      testthat::expect_equal(ca_utility, ref)
    }
  })
})


testthat::test_that("Follow-leader utility wrapper works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      if (!is.null(state$FL[[i]])) {
        ref = unname(m4ma::flUtility(p[i, ], state$FL[[i]], use = 'r'))
        fl_utility = unname(m4ma::flUtility(p[i, ], state$FL[[i]], use = 'cpp'))
        testthat::expect_equal(fl_utility, ref)
      }
    }
  })
})


testthat::test_that("Goal-angle utility wrapper works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = m4ma::gaUtility(p[i, ], state$GA[[i]], use = 'r')
      ga_utility = m4ma::gaUtility(p[i, ], state$GA[[i]], use = 'cpp')
      testthat::expect_equal(ga_utility, ref)
    }
  })
})


testthat::test_that("Interpersonal distance utility wrapper works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = m4ma::idUtility(p[i, ], i, state$ID[[i]], state$ok[[i]],
                            state$group, use = 'r')
      id_utility = m4ma::idUtility(p[i, ], i, state$ID[[i]], state$ok[[i]],
                                   state$group, use = 'cpp')
      testthat::expect_equal(id_utility, ref)
    }
  })
})


testthat::test_that("Preferred speed utility wrapper works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = unname(m4ma::psUtility(p[i, ], state$v[i], state$d[i], use = 'r'))
      ps_utility = unname(m4ma::psUtility(p[i, ], state$v[i],
                                          state$d[i], use = 'cpp'))
      testthat::expect_equal(ps_utility, ref)
    }
  })
})


testthat::test_that("Walk-beside utility wrapper works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      if (!is.null(state$WB[[i]])) {
        ref = unname(m4ma::wbUtility(p[i, ], state$WB[[i]], use = 'r'))
        wb_utility = unname(m4ma::wbUtility(p[i, ], state$WB[[i]], use = 'cpp'))
        testthat::expect_equal(wb_utility, ref)
      } else {
        state$WB[[i]] = list(
          buddies = matrix(runif(12), 3, 4, dimnames = list(c("angleDisagree", "", ""), NULL)),
          dists = matrix(runif(66), 4, 33)
        )
        ref = unname(m4ma::wbUtility(p[i, ], state$WB[[i]], use = 'r'))
        wb_utility = unname(m4ma::wbUtility(p[i, ], state$WB[[i]], use = 'cpp'))
        testthat::expect_equal(wb_utility, ref)
      }
    }
  })
})
