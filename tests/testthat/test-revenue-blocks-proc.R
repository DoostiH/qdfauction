## Tests for R/revenue_patch.R and R/proc_cost_boot_blocks.R.
## NOT YET RUN: written in a session without R. Needs Project A's
## R/procurement.R (proc_values, proc_cost, proc_simulate) in the package.

test_that("fpa_revenue accepts a vector v_s of the same length as r", {
  set.seed(11); n <- 5
  Q <- function(u) stats::qgamma(u, 10, 2)
  q <- function(u) 1 / stats::dgamma(stats::qgamma(u, 10, 2), 10, 2)
  bf <- fpa_bid_function(Q, q, n)
  beta <- function(u) stats::approx(bf$U, bf$bid, xout = u, rule = 2)$y
  B <- matrix(beta(stats::runif(600 * n)), 600, n)
  fit <- fpa_values(B, method = "bernstein", bandwidth = "bcv", monotone = "spline")
  r <- c(4, 5, 6)
  expect_equal(fpa_revenue(fit, r, v_s = r),
               vapply(r, function(x) fpa_revenue(fit, x, v_s = x), numeric(1)))
  expect_equal(fpa_revenue(fit, r, v_s = 2),
               vapply(r, function(x) fpa_revenue(fit, x, v_s = 2), numeric(1)))
  expect_error(fpa_revenue(fit, r, v_s = c(1, 2)), "same length")
})

test_that("proc_cost_boot: 'auctions' is the default and 'blocks' with L = 1 matches it in law", {
  set.seed(21)
  sim <- proc_simulate(100, 3, Qc = function(u) 1 + u, qc = function(u) rep(1, length(u)), M = 3)
  a <- proc_cost_boot(sim$bids, r = c(Inf, 1.8), v0 = 2, B = 49)
  expect_equal(a$resample, "auctions"); expect_equal(a$L, 1L)
  b <- proc_cost_boot(sim$bids, r = c(Inf, 1.8), v0 = 2, B = 49, resample = "blocks", L = 1)
  expect_equal(b$L, 1L)
  expect_equal(a$estimate, b$estimate)
  expect_true(all(is.finite(b$se)) && all(b$se > 0))
  expect_equal(b$se, a$se, tolerance = 0.5)       # same scheme, different draws
})

test_that("proc_cost_boot: L = 'auto' runs and returns an integer block length", {
  skip_if_not_installed("np")
  set.seed(22)
  sim <- proc_simulate(120, 3, Qc = function(u) 1 + u, qc = function(u) rep(1, length(u)), M = 3)
  out <- suppressWarnings(proc_cost_boot(sim$bids, r = 1.8, v0 = 2, B = 19, resample = "blocks"))
  expect_true(is.integer(out$L) && out$L >= 1L)
})
