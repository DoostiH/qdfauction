test_that("kde_hh reproduces the Hickman-Hubbard bandwidth and integrates to about one", {
  set.seed(1)
  x <- rbeta(300, 2, 1)
  n <- length(x)
  g <- seq(0.001, 0.999, length.out = 500)
  f <- kde_hh(g, x)
  expect_equal(attr(f, "h"), 2.978 * (4 / 3)^(1 / 5) * sd(x) * n^(-1 / 5))
  expect_equal(sum(diff(g) * (f[-1] + f[-length(f)]) / 2), 1, tolerance = 0.1)
  expect_lt(mean((f - dbeta(g, 2, 1))^2), 0.08)
  fb <- kde_hh(g, x, side = "both")
  expect_true(all(is.finite(fb)))
  # the lower-boundary correction removes the bias of an uncorrected kernel at 0
  x0 <- rbeta(2000, 1, 1)
  expect_lt(abs(kde_hh(0.01, x0)[1] - 1), 0.3)
})

test_that("Marmer-Shneyerov estimator is sensible on uniform values", {
  set.seed(2)
  sim <- fpa_simulate(100, 5, Q = function(u) u)
  g <- seq(0.1, 0.9, by = 0.05)
  f <- density_marmer_shneyerov(g, sim$bids, 5)
  expect_true(all(is.finite(f)))
  expect_lt(median(abs(f - 1)), 0.3)
})

test_that("fpa_density recovers Beta densities with every second-stage method", {
  set.seed(3)
  sim <- fpa_simulate(40, 5, Q = function(u) qbeta(u, 2, 2))
  g <- seq(0.1, 0.9, by = 0.02); ft <- dbeta(g, 2, 2)
  for (m in c("wavelet_ll", "wavelet_linear", "wavelet_hard", "wavelet_block",
              "kde_hh", "kde", "empirical_quantile", "marmer_shneyerov")) {
    fd <- fpa_density(sim$bids, method = m, x = g, j0 = 5, h = 0.1)
    expect_s3_class(fd, "fpa_density")
    expect_true(all(is.finite(fd$f)), label = m)
    expect_lt(mean((fd$f - ft)^2), switch(m, empirical_quantile = 3, marmer_shneyerov = 0.4, 0.15), label = m)
    expect_equal(fd$fun(g), fd$f)
  }
  expect_equal(fpa_density(sim$bids, x = g)$map$mode, "none")   # values in [0, 1] -> identity map
  expect_output(print(fpa_density(sim$bids, x = g)), "wavelet_ll")
})

test_that("log-scale map integrates correctly on dollar-scale data", {
  data(timber)
  fd <- fpa_density(timber$n3, method = "wavelet_ll", j0 = 4, h = 0.15)
  expect_equal(fd$map$mode, "log")
  xx <- seq(min(fd$values), max(fd$values), length.out = 2000)
  ff <- fd$fun(xx)
  expect_equal(sum(diff(xx) * (ff[-1] + ff[-length(ff)]) / 2), 1, tolerance = 0.15)
  expect_gt(fd$x[which.max(fd$f)], 20); expect_lt(fd$x[which.max(fd$f)], 120)
  fd2 <- fpa_density(timber$n3, method = "kde_hh", scale = "minmax")
  expect_equal(fd2$map$mode, "minmax")
  expect_true(all(fd2$f >= 0))
})

test_that("first stage can come from fpa_values or a qdf method", {
  set.seed(4)
  sim <- fpa_simulate(40, 5, Q = function(u) qbeta(u, 2, 2))
  fv <- fpa_values(sim$bids, method = "bernstein", bandwidth = "bcv")
  fd1 <- fpa_density(sim$bids, method = "kde", first_stage = fv, x = c(0.3, 0.5, 0.7))
  fd2 <- fpa_density(sim$bids, method = "kde", first_stage = "bernstein", bandwidth = "bcv",
                     x = c(0.3, 0.5, 0.7))
  expect_equal(fd1$f, fd2$f)
  expect_equal(fd1$first_stage, "fpa_values"); expect_equal(fd2$first_stage, "bernstein")
})

test_that("cross-validated selection and method comparison run", {
  set.seed(5)
  sim <- fpa_simulate(40, 5, Q = function(u) qbeta(u, 3, 1))
  sel <- select_fpa_density(sim$bids, 5, "wavelet_ll", j0_grid = 3:5, h_grid = c(0.1, 0.3), folds = 3)
  expect_true(sel$j0 %in% 3:5); expect_true(sel$h %in% c(0.1, 0.3))
  fd <- fpa_density(sim$bids, method = "wavelet_ll", j0 = "cv", j0_grid = 3:4, h_grid = 0.2, folds = 2)
  expect_true(fd$j0 %in% 3:4); expect_equal(fd$h, 0.2)
  cv <- fpa_density_cv(sim$bids, methods = c("wavelet_ll", "kde_hh"), folds = 2, j0 = 4)
  expect_named(cv, c("wavelet_ll", "kde_hh")); expect_true(all(is.finite(cv)))
})

test_that("LSCV criterion is available, trimmed, and does not collapse to the smallest h", {
  set.seed(11)
  bf <- fpa_bid_function(function(u) qbeta(u, 3, 1), n = 5)
  H <- c(0.01, 0.03, 0.08, 0.12, 0.3)
  sc <- 0
  for (p in 1:3) {
    b <- as.numeric(fpa_simulate(40, 5, function(u) qbeta(u, 3, 1), bid_function = bf)$bids)
    s <- select_fpa_density(b, 5, "wavelet_ll", "gpv_hh", 4:6, H, folds = 3, criterion = "lscv", scale = "none")
    expect_equal(s$criterion, "lscv"); expect_equal(dim(s$ll), c(3, 5))
    sc <- sc + s$ll
  }
  best <- which(sc == max(sc), arr.ind = TRUE)[1, ]
  expect_gt(H[best[2]], 0.01)
  fd <- fpa_density(b, 5, "wavelet_ll", j0 = "cv", criterion = "lscv", j0_grid = 4:5, h_grid = c(0.08, 0.2),
                    folds = 2, scale = "none")
  expect_true(fd$j0 %in% 4:5)
  expect_equal(fd$selection$criterion, "lscv")
  s2 <- select_fpa_density(b, 5, "wavelet_ll", "gpv_hh", 4:5, c(0.08, 0.2), folds = 2, criterion = "loglik", scale = "none")
  expect_equal(s2$criterion, "loglik")
})
