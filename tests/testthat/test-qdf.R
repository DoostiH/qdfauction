set.seed(1)
x <- rgamma(200, 5, 1)
u <- seq(0.05, 0.95, by = 0.05)
q_true <- 1 / dgamma(qgamma(u, 5, 1), 5, 1)

methods <- c("kernel", "kernel_corrected", "poisson", "bernstein", "jones",
             "soni", "indirect_poisson")

test_that("all estimators are close to the truth away from the boundaries", {
  for (m in methods) {
    h <- if (m %in% c("jones", "soni")) 0.5 else 0.05
    est <- switch(m,
                  soni = qdf_soni(u, x, h, H = 0.05),
                  do.call(paste0("qdf_", m), list(u, x, h)))
    rel <- abs(est - q_true) / q_true
    expect_lt(median(rel), 0.25, label = paste(m, "median relative error"))
  }
})

test_that("qdf() runs for every method and selector and returns a valid object", {
  set.seed(2)
  xs <- rgamma(40, 5, 1)
  for (m in methods) for (s in c("bcv", "rlcv", "wbcv")) {
    fit <- qdf(xs, method = m, bandwidth = s, u = c(0.25, 0.5, 0.75))
    expect_s3_class(fit, "qdf")
    expect_length(fit$q, 3)
    expect_true(all(is.finite(fit$q)))
    expect_true(fit$h > 0)
    expect_equal(fit$selection$selector, s)
  }
})

test_that("fixed bandwidth, predict and print work", {
  fit <- qdf(x, "bernstein", bandwidth = 0.04)
  expect_null(fit$selection)
  expect_equal(fit$h, 0.04)
  expect_equal(predict(fit), fit$q)
  expect_equal(predict(fit, u = 0.5), qdf_bernstein(0.5, x, 0.04))
  expect_output(print(fit), "bernstein")
})

test_that("input validation", {
  expect_error(qdf(c(x, -1), "poisson"), "positive")
  expect_error(qdf(x, "kernel", bandwidth = -1), "positive number")
  expect_error(qdf_kernel(c(0.5, 1.5), x, 0.05), "\\[0, 1\\]")
  expect_error(qdf(x[1:3]), "at least 5")
})

test_that("loo = FALSE approximation is close to exact for a linear smoother", {
  set.seed(3)
  xs <- rgamma(100, 5, 1)
  a <- bcv_criterion(20, xs, qdf_bernstein, loo = TRUE)
  b <- bcv_criterion(20, xs, qdf_bernstein, loo = FALSE)
  expect_lt(abs(a - b) / abs(a), 0.15)
})

test_that("grid search never does worse than plain optimize", {
  set.seed(4)
  xs <- rgamma(60, 5, 1)
  s0 <- select_bandwidth(xs, qdf_kernel, "bcv", c(10, 50), grid = 0)
  s1 <- select_bandwidth(xs, qdf_kernel, "bcv", c(10, 50), grid = 15)
  expect_lte(s1$criterion, s0$criterion + 1e-8)
})
