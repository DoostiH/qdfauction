test_that("vectorised variance equals the brute-force triple sum of eq. (4.2)", {
  set.seed(3)
  b <- sort(runif(25)) * 0.8; n <- 5; h_g <- 0.15; h_f <- 0.2; grid <- c(0.35, 0.45)
  fs <- qdfauction:::gpv_first_stage(b, n, h_g)
  vec <- qdfauction:::gpv_variance(fs, grid, n, h_g, h_f)
  tw4 <- qdfauction:::.tw4; twd <- qdfauction:::.tw_d; N <- length(b)
  bf <- sapply(grid, function(v) {
    eta <- function(i, j) if (!fs$trim[j]) 0 else
      twd((fs$V[j] - v) / h_f) * fs$G[j] / fs$g[j]^2 * tw4((b[i] - b[j]) / h_g)
    s <- 0
    for (i in 1:N) for (j in 1:N) if (j != i) for (jp in 1:N) if (jp != i && jp != j)
      s <- s + eta(i, j) * eta(i, jp)
    max(s, 0) / (N * (N - 1) * (N - 2)) / (n * (n - 1)^2 * h_f^2 * h_g) })
  expect_equal(vec, bf, tolerance = 1e-10)
})

test_that("weighted bootstrap first stage equals recomputation on the resample", {
  set.seed(4)
  b <- runif(60) * 0.8; n <- 3; h_g <- 0.12
  fs <- qdfauction:::gpv_first_stage(b, n, h_g)
  idx <- sample.int(60, 60, replace = TRUE); w <- tabulate(idx, 60)
  fw <- qdfauction:::gpv_first_stage_w(w, b, n, h_g, fs$Kg, order(b), rep(TRUE, 60), FALSE)
  fd <- qdfauction:::gpv_first_stage(b[idx], n, h_g)
  expect_equal(fw$V[idx], fd$V)
})

test_that("fpa_density_ci returns a valid estimate, intervals and band", {
  set.seed(5)
  sim <- fpa_simulate(100, 5, function(u) u)
  ci <- fpa_density_ci(sim$bids, grid = seq(0.25, 0.75, by = 0.05), B = 49)
  expect_s3_class(ci, "fpa_density_ci")
  expect_equal(ci$pointwise, "percentile")
  expect_true(all(ci$se > 0))
  expect_true(all(ci$band_lower <= ci$f & ci$band_upper >= ci$f))
  expect_true(ci$crit > stats::qnorm(0.975))
  expect_lt(max(abs(ci$f - 1)), 0.5)
  expect_output(print(ci), "uniform band")
  ci2 <- fpa_density_ci(sim$bids, grid = c(0.4, 0.6), B = 19, pointwise = "normal",
                        constants = c(3.15, 3.72), sd_v = "all")
  expect_equal(ci2$upper - ci2$f, stats::qnorm(0.975) * ci2$se)
  expect_error(fpa_density_ci(sim$bids[1:10], n = 5), "20 bids")
})
