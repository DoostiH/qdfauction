# Tolerances are taken from the Project A runs (qdfauction 0.2.0, R 4.6).

Qc_U <- function(u) 1 + u                 # costs U[1, 2]
qc_U <- function(u) rep(1, length(u))

test_that("simulated procurement bids match the closed form b = c + (2 - c)/n", {
  set.seed(1)
  for (n in c(2, 3, 5)) {
    s <- proc_simulate(200, n, Qc_U, qc_U, M = 3)
    expect_lt(max(abs(s$bids - (s$costs + (2 - s$costs) / n))), 1e-5)  # observed 5.5e-7 at n = 3
  }
})

test_that("reflection is exact for the boundary-corrected kernel", {
  set.seed(2)
  s <- proc_simulate(100, 3, Qc_U, qc_U, M = 3)
  b <- sort(as.numeric(s$bids)); N <- length(b); u <- seq_len(N) / (N + 1)
  pv <- proc_values(s$bids, method = "kernel_corrected", bandwidth = 0.05, monotone = "none")
  direct <- b - (1 - u) * qdf(b, method = "kernel_corrected", bandwidth = 0.05, u = u)$q / 2
  expect_lt(max(abs(pv$costs_raw - direct)), 1e-10)                  # observed <= 1.2e-14
})

test_that("reflection is exact in the interior for the Jones kernel", {
  set.seed(3)
  s <- proc_simulate(100, 3, Qc_U, qc_U, M = 3)
  b <- sort(as.numeric(s$bids)); N <- length(b); u <- seq_len(N) / (N + 1)
  pv <- proc_values(s$bids, method = "kernel", bandwidth = 0.05, monotone = "none")
  direct <- b - (1 - u) * qdf(b, method = "kernel", bandwidth = 0.05, u = u)$q / 2
  k <- pv$u > 0.3 & pv$u < 0.7
  expect_lt(max(abs(pv$costs_raw - direct)[k]), 1e-5)                # observed 6.4e-8
})

test_that("proc_values returns ascending bids and costs on the original scale", {
  set.seed(4)
  s  <- proc_simulate(100, 3, Qc_U, qc_U, M = 3)
  pv <- proc_values(s$bids, method = "bernstein", bandwidth = "sqrt")
  expect_length(pv$bids, 300); expect_length(pv$costs, 300)
  expect_false(anyNA(pv$costs))
  expect_true(all(diff(pv$bids) >= 0))
  expect_true(all(diff(pv$costs) >= -1e-12))                         # isotonic default
  expect_equal(dim(pv$costs_mat), c(100, 3))
  expect_lt(mean((pv$costs_mat - s$costs)^2), 0.005)                 # MSEP ~ 4e-4 in the MC
})

test_that("with ceilings at or below v0, expected cost never exceeds v0", {
  set.seed(5)
  s  <- proc_simulate(100, 3, Qc_U, qc_U, M = 3)
  pv <- proc_values(s$bids, method = "bernstein", bandwidth = "sqrt")
  C  <- proc_cost(pv, r = seq(1.1, 1.5, by = 0.1), v0 = 1.5)   # every payment <= max(r, v0) = v0
  expect_true(all(C <= 1.5 + 1e-8))
})

test_that("no-ceiling expected cost is close to the mean winning bid (independence)", {
  set.seed(6)
  s  <- proc_simulate(200, 3, Qc_U, qc_U, M = 3)
  pv <- proc_values(s$bids, method = "bernstein", bandwidth = "sqrt")
  chk <- proc_check_winning_bid(pv, s$bids)
  expect_lt(abs(chk[["difference"]]), 0.03)                          # RMSE 0.007 at L = 200
})

test_that("parametric families without reflection symmetry are rejected", {
  s <- proc_simulate(50, 3, Qc_U, qc_U, M = 3)
  expect_error(proc_values(s$bids, family = "clayton"), "survival")
  expect_error(proc_values(s$bids, family = "gumbel"), "survival")
})

test_that("fpa_values gives the same result with and without explicit n = NULL", {
  set.seed(7)
  s <- proc_simulate(50, 3, Qc_U, qc_U, M = 3)
  B <- 3 - s$bids
  a <- fpa_values(B, method = "bernstein", bandwidth = "sqrt")
  b <- fpa_values(B, n = NULL, method = "bernstein", bandwidth = "sqrt")
  for (k in c("bids", "u", "q", "values")) expect_equal(a[[k]], b[[k]])
})

test_that("Frank ceiling estimate is not biased under independence at n = 5", {
  skip_on_cran()
  set.seed(8)
  ## true optimal ceiling for U[1,2] costs, v0 = 1.5: r + (r - 1) = 1.5 -> r* = 1.25
  bias <- replicate(50, {
    s  <- proc_simulate(100, 5, Qc_U, qc_U, M = 3)
    pv <- proc_values(s$bids, method = "bernstein", bandwidth = "sqrt", family = "frank")
    proc_optimal_ceiling(pv, v0 = 1.5)$r_star - 1.25
  })
  expect_lt(abs(mean(bias)), 0.05)        # 0.2.0: +0.23 (500 reps); fails until item 5 is fixed
})

test_that("bid_icc is near zero without and large with a common component", {
  set.seed(9)
  id <- rep(1:200, each = 3)
  expect_lt(abs(bid_icc(rnorm(600), id)), 0.1)
  expect_gt(bid_icc(rep(rnorm(200), each = 3) + rnorm(600, sd = 0.5), id), 0.6)
})
