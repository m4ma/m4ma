# Test rcpp utility functions against R counterparts on test trace

test_filepath = file.path('data', 'trace_i.rda')

test_obj_name = load(test_filepath)

p = attr(get(test_obj_name), 'pMat')


testthat::test_that("Blocked-angle utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      if (!is.null(state$FL[[i]])) {
        ref = m4ma::baUtility_r(p[i, ], state$BA[[i]])
        ba_utility = m4ma::baUtility_rcpp(
          p[i, "aBA"],
          p[i, "bBA"],
          state$BA[[i]],
          as.numeric(names(state$BA[[i]])) - 1
        )
        testthat::expect_equal(ba_utility, ref)
      }
    }
  })
})


testthat::test_that("Current-angle utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = unname(m4ma::caUtility_r(p[i, ]))
      ca_utility = m4ma::caUtility_rcpp(
        p[i, "aCA"],
        p[i, "bCA"],
        p[i, "bCAlr"]
      )
      testthat::expect_equal(ca_utility, ref)
    }
  })
})


testthat::test_that("Follow-leader utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      if (!is.null(state$FL[[i]])) {
        ref = unname(m4ma::flUtility_r(p[i, ], state$FL[[i]]))
        fl_utility = m4ma::flUtility_rcpp(
          p[i, "aFL"],
          p[i, "bFL"],
          p[i, "dFL"],
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
      ref = m4ma::gaUtility_r(p[i, ], state$GA[[i]])
      ga_utility = m4ma::gaUtility_rcpp(
        p[i, "bGA"],
        p[i, "aGA"],
        state$GA[[i]]
      )
      testthat::expect_equal(ga_utility, ref)
    }
  })
})


testthat::test_that("Interpersonal distance utility computation works", {
  lapply(get(test_obj_name), function(state) {
    for (i in 1:length(state$v)) {
      ref = m4ma::idUtility_r(p[i, ], i, state$ID[[i]], state$ok[[i]], state$group)
      id_utility = m4ma::idUtility_rcpp(
        p[i, "bID"],
        p[i, "dID"],
        p[i, "aID"],
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
      ref = unname(m4ma::psUtility_r(p[i, ], state$v[i], state$d[i]))
      ps_utility = m4ma::psUtility_rcpp(
        p[i, "aPS"],
        p[i, "bPS"],
        p[i, "sPref"],
        p[i, "sSlow"],
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
        ref = unname(m4ma::wbUtility_r(p[i, ], state$WB[[i]]))
        wb_utility = m4ma::wbUtility_rcpp(
          p[i, "aWB"],
          p[i, "bWB"],
          state$WB[[i]]$buddies,
          state$WB[[i]]$dists
        )
        testthat::expect_equal(wb_utility, ref)
      } else {
        state$WB[[i]] = list(
          buddies = matrix(runif(12), 3, 4, dimnames = list(c("angleDisagree", "", ""), NULL)),
          dists = matrix(runif(66), 4, 33)
        )
        ref = unname(m4ma::wbUtility(p[i, ], state$WB[[i]]))
        wb_utility = m4ma::wbUtility_rcpp(
          p[i, "aWB"],
          p[i, "bWB"],
          state$WB[[i]]$buddies,
          state$WB[[i]]$dists
        )
        testthat::expect_equal(wb_utility, ref)
      }
    }
  })
})


testthat::test_that("Utility computation works", {
  ref = c(-300907370, -301091391, -301180594)
  
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
