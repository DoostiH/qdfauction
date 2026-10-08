test_that("periodised wavelet basis is orthonormal", {
  z <- (seq_len(4000) - 0.5) / 4000
  B <- cbind(phi_per(z, 3, 0:7), psi_per(z, 3, 0:7))
  expect_lt(max(abs(crossprod(B) / length(z) - diag(16))), 1e-4)
})

test_that("nonlinear wavelet QDF estimators (CSDA 2016 Thms 3.2/3.4; Shirazi-Doosti 2022)", {
  set.seed(1)
  x <- rbeta(400, 0.5, 0.5); u <- seq(0.1, 0.9, by = 0.05)
  qt <- 1 / dbeta(qbeta(u, 0.5, 0.5), 0.5, 0.5)
  for (r in c("hard_lipschitz", "hard_general", "block")) {
    q <- qdf_wavelet_nl(u, x, r, smooth = TRUE)
    expect_true(all(is.finite(q)), label = r)
    expect_lt(mean((q - qt)^2), 0.1, label = r)
    expect_false(is.null(attr(qdf_wavelet_nl(u, x, r), "kept")))
  }
  expect_error(qdf_wavelet_nl(u, x, tau = 1), "tau")
})

test_that("Theorem 2 estimator of the wavelet density paper runs", {
  set.seed(2)
  sim <- fpa_simulate(40, 5, function(u) qbeta(u, 2, 2))
  fd <- fpa_density(sim$bids, method = "wavelet_theory", x = seq(0.1, 0.9, by = 0.1))
  expect_true(all(is.finite(fd$f)))
  expect_lt(mean((fd$f - dbeta(fd$x, 2, 2))^2), 0.3)
})

test_that("nonparametric copula multiplier is consistent under independence", {
  set.seed(3)
  U <- rcopula_archimedean(1500, 3, "independence")
  u <- c(0.3, 0.5, 0.7)
  expect_lt(max(abs(copula_psi_np(U, u) - u / 2)), 0.05)
})

test_that("copula specification test rejects a wrong family and keeps a right one", {
  set.seed(4)
  U <- rcopula_archimedean(200, 3, "gumbel", theta = 2.5)
  expect_lt(copula_spec_test(U, "independence", B = 49)$p.value, 0.1)
  t_ok <- copula_spec_test(U, "gumbel", B = 49)
  expect_s3_class(t_ok, "htest"); expect_gt(t_ok$p.value, 0.05)
})

test_that("flexible (nonparametric) copula inversion and its revenue guard", {
  set.seed(5)
  sim <- fpa_simulate(60, 3, function(u) 1 + 2 * u, family = "clayton", theta = 2)
  fv <- fpa_values(sim$bids, method = "bernstein", bandwidth = 0.05, family = "nonparametric")
  expect_equal(fv$copula$family, "nonparametric")
  tv <- as.numeric(sim$values)[fv$order]
  expect_lt(median(abs(fv$values - tv)), 0.15)
  expect_error(fpa_revenue_qd(fv, p = 0.5), "diagonal section")
  expect_error(fpa_values(as.numeric(sim$bids), n = 3, family = "nonparametric"), "matrix")
})

test_that("bootstrap bands for any density estimator", {
  set.seed(6)
  sim <- fpa_simulate(40, 4, function(u) qbeta(u, 2, 2))
  bb <- fpa_density_boot(sim$bids, method = "marmer_shneyerov", x = seq(0.2, 0.8, by = 0.1), B = 19)
  expect_s3_class(bb, "fpa_density_boot")
  expect_true(all(bb$band_lower <= bb$f & bb$band_upper >= bb$f, na.rm = TRUE))
  expect_output(print(bb), "critical value")
  bb2 <- fpa_density_boot(sim$bids, method = "kde_hh", x = c(0.4, 0.6), B = 9, resample = "auctions")
  expect_equal(bb2$resample, "auctions")
})
