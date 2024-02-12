# Test rcpp utility functions against R counterparts on test trace

test_filepath = testthat::test_path('data', 'trace_i.rda')

test_obj_name = load(test_filepath)

p = exp(attr(get(test_obj_name), 'pMat'))


testthat::test_that("Blocked-angle utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      if (!is.null(state$FL[[i]])) {
        ref = m4ma::baUtility_r(c(bBA = 4, aBA = 2), state$BA[[i]])
        ba_utility = m4ma::baUtility_rcpp(
          2,
          4,
          state$BA[[i]],
          as.numeric(names(state$BA[[i]])) - 1
        )
        testthat::expect_equal(ba_utility, ref)
      }
    }
  })
})

testthat::test_that("Blocked-angle utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      if (!is.null(state$FL[[i]])) {
        ref = m4ma::baUtility_r(exp(c(bBA = 4, aBA = 2)), state$BA[[i]])
        ba_utility = m4ma::baUtility_rcpp(
          exp(2),
          exp(4),
          state$BA[[i]],
          as.numeric(names(state$BA[[i]])) - 1
        )
        testthat::expect_equal(ba_utility, ref)
      }
    }
  })
})


testthat::test_that("Current-angle utility computation works (exponential scale)", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = unname(m4ma::caUtility_r(c(aCA = 1, bCA = 2, bCAlr = 10)))
      ca_utility = m4ma::caUtility_rcpp(
        1,
        2,
        10
      )
      testthat::expect_equal(ca_utility, ref)
    }
  })
})


testthat::test_that("Current-angle utility computation works (exponential scale)", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = unname(m4ma::caUtility_r(exp(c(aCA = 1, bCA = 2, bCAlr = 10))))
      ca_utility = m4ma::caUtility_rcpp(
        exp(1),
        exp(2),
        exp(10)
      )
      testthat::expect_equal(ca_utility, ref)
    }
  })
})


testthat::test_that("Follow-leader utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      if (!is.null(state$FL[[i]])) {
        ref = unname(m4ma::flUtility_r(c(aFL = 2, bFL = 1, dFL = 0), state$FL[[i]]))
        fl_utility = m4ma::flUtility_rcpp(
          2,
          1,
          0,
          state$FL[[i]]$leaders,
          state$FL[[i]]$dists
        )
        testthat::expect_equal(fl_utility, ref)
      }
    }
  })
})


testthat::test_that("Follow-leader utility computation works (exponential scale)", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      if (!is.null(state$FL[[i]])) {
        ref = unname(m4ma::flUtility_r(exp(c(aFL = 2, bFL = 1, dFL = 0)), state$FL[[i]]))
        fl_utility = m4ma::flUtility_rcpp(
          exp(2),
          exp(1),
          exp(0),
          state$FL[[i]]$leaders,
          state$FL[[i]]$dists
        )
        testthat::expect_equal(fl_utility, ref)
      }
    }
  })
})


testthat::test_that("Goal-angle utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = m4ma::gaUtility_r(c(bGA = 10, aGA = 2), state$GA[[i]])
      ga_utility = m4ma::gaUtility_rcpp(
        10,
        2,
        state$GA[[i]]
      )
      testthat::expect_equal(ga_utility, ref)
    }
  })
})


testthat::test_that("Goal-angle utility computation works (exponential scale)", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = m4ma::gaUtility_r(exp(c(bGA = 10, aGA = 2)), state$GA[[i]])
      ga_utility = m4ma::gaUtility_rcpp(
        exp(10),
        exp(2),
        state$GA[[i]]
      )
      testthat::expect_equal(ga_utility, ref)
    }
  })
})


testthat::test_that("Interpersonal distance utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = m4ma::idUtility_r(c(bID = 1, aID = 2, dID = 0), i, state$ID[[i]], state$ok[[i]], state$group)
      id_utility = m4ma::idUtility_rcpp(
        1,
        0,
        2,
        i - 1,
        state$ok[[i]],
        state$group,
        state$ID[[i]]
      )
      testthat::expect_equal(id_utility, ref)
    }
  })
})


testthat::test_that("Interpersonal distance utility computation works (exponential scale)", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = m4ma::idUtility_r(exp(c(bID = 1, aID = 2, dID = 0)), i, state$ID[[i]], state$ok[[i]], state$group)
      id_utility = m4ma::idUtility_rcpp(
        exp(1),
        exp(0),
        exp(2),
        i - 1,
        state$ok[[i]],
        state$group,
        state$ID[[i]]
      )
      testthat::expect_equal(id_utility, ref)
    }
  })
})


testthat::test_that("Preferred speed utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = unname(m4ma::psUtility_r(c(aPS = 2, bPS = 4, sPref = 1, sSlow = 1), state$v[i], state$d[i]))
      ps_utility = m4ma::psUtility_rcpp(
        2,
        4,
        1,
        1,
        state$v[i],
        state$d[i]
      )
      testthat::expect_equal(ps_utility, ref)
    }
  })
})


testthat::test_that("Preferred speed utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = unname(m4ma::psUtility_r(exp(c(aPS = 2, bPS = 4, sPref = 1, sSlow = 1)), state$v[i], state$d[i]))
      ps_utility = m4ma::psUtility_rcpp(
        exp(2),
        exp(4),
        exp(1),
        exp(1),
        state$v[i],
        state$d[i]
      )
      testthat::expect_equal(ps_utility, ref)
    }
  })
})


testthat::test_that("Walk-beside utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      if (!is.null(state$WB[[i]])) {
        ref = unname(m4ma::wbUtility_r(c(bWB = 0, aWB = 0), state$WB[[i]]))
        wb_utility = m4ma::wbUtility_rcpp(
          0,
          0,
          state$WB[[i]]$buddies,
          state$WB[[i]]$dists
        )
        testthat::expect_equal(wb_utility, ref)
      } else {
        state$WB[[i]] = list(
          buddies = matrix(runif(12), 3, 4, dimnames = list(c("angleDisagree", "", ""), NULL)),
          dists = matrix(runif(66), 4, 33)
        )
        ref = unname(m4ma::wbUtility(c(bWB = 0, aWB = 0), state$WB[[i]]))
        wb_utility = m4ma::wbUtility_rcpp(
          0,
          0,
          state$WB[[i]]$buddies,
          state$WB[[i]]$dists
        )
        testthat::expect_equal(wb_utility, ref)
      }
    }
  })
})


testthat::test_that("Walk-beside utility computation works (exponential scale)", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      if (!is.null(state$WB[[i]])) {
        ref = unname(m4ma::wbUtility_r(exp(c(bWB = 0, aWB = 0)), state$WB[[i]]))
        wb_utility = m4ma::wbUtility_rcpp(
          1,
          1,
          state$WB[[i]]$buddies,
          state$WB[[i]]$dists
        )
        testthat::expect_equal(wb_utility, ref)
      } else {
        state$WB[[i]] = list(
          buddies = matrix(runif(12), 3, 4, dimnames = list(c("angleDisagree", "", ""), NULL)),
          dists = matrix(runif(66), 4, 33)
        )
        ref = unname(m4ma::wbUtility(exp(c(bWB = 0, aWB = 0)), state$WB[[i]]))
        wb_utility = m4ma::wbUtility_rcpp(
          1,
          1,
          state$WB[[i]]$buddies,
          state$WB[[i]]$dists
        )
        testthat::expect_equal(wb_utility, ref)
      }
    }
  })
})


testthat::test_that("Utility computation works", {
  ref = c(-19915374, -29604211, -31746385)

  u = sapply(1:length(get(test_obj_name)), function(j) {
    state = get(test_obj_name)[[j]]
    out = sapply(1:length(state$v), function(i) {
      u = m4ma::utility(
        p[i, ], i - 1, state$v[i], state$d[i], state$BA[[i]], state$GA[[i]],
        state$ID[[i]], state$FL[[i]], state$WB[[i]], state$ok[[i]],
        state$group
      )

      return(sum(u[is.finite(u)]))
    })

    return(sum(out))
  })

  testthat::expect_equal(u, ref)
})
