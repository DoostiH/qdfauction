## Tests for R/blocks.R. NOT YET RUN: written in a session without R.
## Run devtools::test() in the package session before release.

sim_ar1_bids <- function(TT, n, rho, seed) {
  set.seed(seed)
  Q <- function(u) stats::qgamma(u, 10, 2)
  q <- function(u) 1 / stats::dgamma(stats::qgamma(u, 10, 2), 10, 2)
  bf <- fpa_bid_function(Q, q, n)
  beta <- function(u) stats::approx(bf$U, bf$bid, xout = u, rule = 2)$y
  Z <- matrix(0, TT, n); Z[1, ] <- stats::rnorm(n)
  for (t in 2:TT) Z[t, ] <- rho * Z[t - 1, ] + sqrt(1 - rho^2) * stats::rnorm(n)
  matrix(beta(as.numeric(stats::pnorm(Z))), TT, n)
}

test_that("resample_blocks returns valid row indices", {
  set.seed(1)
  idx <- resample_blocks(50, 1)
  expect_length(idx, 50); expect_true(all(idx >= 1 & idx <= 50))
  idx <- resample_blocks(50, 7)
  expect_length(idx, 50); expect_true(all(idx >= 1 & idx <= 50))
  ## within a block, consecutive rows (circularly)
  d <- diff(idx[1:7]) %% 50
  expect_true(all(d == 1))
})

test_that("block_length_auto is small under independence", {
  skip_if_not_installed("np")
  B <- sim_ar1_bids(200, 5, 0, seed = 2)
  L <- block_length_auto(B)
  expect_true(L >= 1 && L <= 8)
  expect_false(is.null(attr(L, "L_pw")))
  expect_false(is.null(attr(L, "T_over_Lpw")))
})

test_that("block_length_auto grows with dependence and is capped at T/2", {
  skip_if_not_installed("np")
  B <- sim_ar1_bids(200, 5, 0.7, seed = 3)
  L <- suppressWarnings(block_length_auto(B))
  expect_true(L >= 4 && L <= 100)
})

test_that("block_length_auto warns when the sample is short for its dependence", {
  skip_if_not_installed("np")
  B <- sim_ar1_bids(40, 5, 0.95, seed = 4)
  expect_warning(block_length_auto(B), "T / L_pw")
})

test_that("winner = 'min' uses the lowest bid", {
  skip_if_not_installed("np")
  B <- sim_ar1_bids(150, 4, 0.5, seed = 5)
  K <- 1.05 * max(B)
  ## reflecting the matrix and using max must equal using min on the original
  expect_equal(attr(block_length_auto(K - B, winner = "max"), "L_pw"),
               attr(block_length_auto(B, winner = "min"), "L_pw"))
})

test_that("fpa_block_boot runs, holds the bandwidth fixed, and reports L", {
  skip_if_not_installed("np")
  B <- sim_ar1_bids(100, 5, 0.5, seed = 6)
  out <- fpa_revenue_boot_blocks(B, r = c(4, 5), v_s = 2, B_boot = 19)
  expect_length(out$point, 2)
  expect_true(all(is.finite(out$boot_sd)) && all(out$boot_sd > 0))
  expect_equal(out$failures, 0)
  expect_true(out$L >= 1)
  ## L = 1 is the i.i.d. auction bootstrap
  out1 <- fpa_revenue_boot_blocks(B, r = 4, v_s = 2, L = 1, B_boot = 19)
  expect_equal(out1$L, 1L)
})
