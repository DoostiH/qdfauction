test_that("filters are orthonormal scaling filters", {
  for (f in c("haar", "db2", "db3", "db4", "sym4", "coif1", "coif2")) {
    h <- wavelet_filter(f)
    expect_equal(sum(h), sqrt(2), tolerance = 1e-8, label = f)
    expect_equal(sum(h^2), 1, tolerance = 1e-8, label = f)
    # double-shift orthogonality
    for (s in seq_len(length(h) / 2 - 1))
      expect_equal(sum(h[-(1:(2 * s))] * h[1:(length(h) - 2 * s)]), 0, tolerance = 1e-7,
                   label = paste(f, "shift", s))
  }
})

test_that("Daubechies-Lagarias scaling function integrates to one and partitions unity", {
  for (f in c("db3", "sym4", "coif1")) {
    tab <- scaling_function(f)
    expect_equal(sum(diff(tab$x) * (tab$phi[-1] + tab$phi[-length(tab$phi)]) / 2), 1,
                 tolerance = 1e-6)
    z <- seq(0.1, 0.9, by = 0.1)
    pu <- rowSums(phi_jk(z + 3, 0, -(tab$N):(tab$N + 3), f))  # sum_k phi(x - k) = 1
    expect_equal(pu, rep(1, length(z)), tolerance = 1e-3, label = f)
  }
})

test_that("periodic DWT is invertible and orthonormal", {
  set.seed(1)
  for (f in c("haar", "coif1", "coif2", "db3")) {
    x <- rnorm(64)
    w <- dwt_periodic(x, 3, f)
    expect_equal(idwt_periodic(w, 3, f), x, tolerance = 1e-10, label = f)
    expect_equal(sum(w^2), sum(x^2), tolerance = 1e-10, label = f)
  }
})

test_that("thresholding rules behave as expected", {
  set.seed(2)
  w <- dwt_periodic(rnorm(32), 3, "coif1")
  wh <- threshold_coefficients(w, 3, n = 200, type = "hard")
  wb <- threshold_coefficients(w, 3, n = 200, type = "block")
  expect_equal(wh[1:4], w[1:4]); expect_equal(wb[1:4], w[1:4])   # smooth part untouched
  expect_true(all(wh == w | wh == 0)); expect_true(all(wb == w | wb == 0))
  wj <- threshold_coefficients(w, 3, n = 200, type = "block", rule = "james-stein")
  expect_true(all(abs(wj) <= abs(w) + 1e-12))
})

test_that("wavelet QDF estimators recover Beta(0.5,0.5) and GLD quantile densities", {
  set.seed(3)
  # Beta(0.5, 0.5): design of Chesneau, Dewan & Doosti (2016), Table 1
  x <- rbeta(500, 0.5, 0.5)
  u <- seq(0.05, 0.95, by = 0.05)
  qt <- 1 / dbeta(qbeta(u, 0.5, 0.5), 0.5, 0.5)
  for (thr in c("none", "hard", "block")) {
    q <- qdf_wavelet(u, x, j0 = 5, h = 0.15, threshold = thr, smooth = TRUE)
    mise <- mean((q - qt)^2)
    expect_lt(mise, 0.05, label = paste("beta", thr))
  }
  # GLD(0.5, 1, 2, 6): quantile density has a closed form
  set.seed(4)
  v <- runif(500); l <- c(0.5, 1, 2, 6)
  y <- l[1] + (v^l[3] - (1 - v)^l[4]) / l[2]
  qg <- (l[3] * u^(l[3] - 1) + l[4] * (1 - u)^(l[4] - 1)) / l[2]
  q <- qdf_wavelet(u, y, j0 = 5, h = 0.15, threshold = "block", smooth = TRUE)
  expect_lt(median(abs(q - qg) / qg), 0.25)
  # exact and Simpson coefficient estimates agree closely
  q1 <- qdf_wavelet(u, y, integration = "exact", smooth = FALSE)
  q2 <- qdf_wavelet(u, y, integration = "simpson", smooth = FALSE)
  expect_lt(max(abs(q1 - q2)) / max(qg), 0.1)
})

test_that("qdf() wavelet branch and CV selection work", {
  set.seed(5)
  x <- rbeta(200, 0.5, 0.5)
  fit <- qdf(x, method = "wavelet_block", bandwidth = 0.15, u = c(0.25, 0.5, 0.75))
  expect_s3_class(fit, "qdf"); expect_equal(fit$j0, 5L); expect_true(all(is.finite(fit$q)))
  expect_equal(predict(fit, 0.5), fit$q[2])
  expect_output(print(fit), "wavelet_block")
  fit2 <- qdf(x, method = "wavelet", bandwidth = "cv", u = c(0.25, 0.5, 0.75),
              j0_grid = 4:5, h_grid = c(0.1, 0.2), folds = 3)
  expect_true(fit2$j0 %in% 4:5); expect_true(fit2$h %in% c(0.1, 0.2))
  expect_error(qdf(x, "wavelet", bandwidth = "bcv"), "cv")
})
