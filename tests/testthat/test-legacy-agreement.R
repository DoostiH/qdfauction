set.seed(20250925)
n  <- 50
Y  <- rgamma(n + 1, 5, 1)
X  <- sort((Y[1:n] + Y[2:(n + 1)]) / 2)   # the DGP of the original script
u  <- c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99)
ni <- (1:n) / n
Xsd <- diff(c(0, X))

test_that("direct kernel estimators reproduce the original functions", {
  expect_equal(qdf_kernel(u, X, 0.05), legacy$dqk(u, X, 0.05))
  expect_equal(qdf_kernel_corrected(u, X, 0.05), legacy$dqkc(u, X, 0.05))
  expect_equal(qdf_kernel(ni, X, 1 / 23.7), legacy$dqk(ni, X, 1 / 23.7))
})

test_that("series estimators reproduce the original functions", {
  for (h in c(1 / 12, 1 / 25.5, 1 / 40)) {
    expect_equal(qdf_poisson(u, X, h),   legacy$dqp(u, X, h))
    expect_equal(qdf_bernstein(u, X, h), legacy$dqb(u, X, h))
  }
})

test_that("indirect estimators reproduce the original functions", {
  expect_equal(qdf_jones(u, X, 0.4), legacy$idqfj(u, X, 0.4))
  # Soni: closed-form Gaussian integral vs numerical integrate()
  expect_equal(qdf_soni(u, X, 0.4, H = 0.07), legacy$idqfST(u, X, 0.4, 0.07),
               tolerance = 1e-8)
  expect_equal(qdf_indirect_poisson(u, X, 1 / 20, invert = "uniroot"), legacy$idqfp(u, X, 1 / 20),
               tolerance = 1e-6)
  # Newton inversion: same estimator, more accurate root than uniroot's default tol
  expect_equal(qdf_indirect_poisson(u, X, 1 / 20, invert = "newton"), legacy$idqfp(u, X, 1 / 20), tolerance = 1e-3)
  expect_equal(qdf_indirect_poisson(u, X, 1 / 20, invert = "newton"),
               qdf_indirect_poisson(u, X, 1 / 20, tol = 1e-12, invert = "uniroot"), tolerance = 1e-8)
})

test_that("BCV criteria reproduce the original functions", {
  m <- 22.3
  expect_equal(bcv_criterion(m, X, qdf_kernel),           legacy$BCVfK(m, X, Xsd, ni, n))
  expect_equal(bcv_criterion(m, X, qdf_kernel_corrected), legacy$BCVfKC(m, X, Xsd, ni, n))
  expect_equal(bcv_criterion(m, X, qdf_poisson),          legacy$BCVfP(m, X, Xsd, ni, n))
  expect_equal(bcv_criterion(m, X, qdf_bernstein),        legacy$BCVfB(m, X, Xsd, ni, n))
  expect_equal(bcv_criterion(4, X, qdf_jones),            legacy$BCVfJ(4, X, Xsd, ni, n))
  soni <- function(u, x, h) qdf_soni(u, x, h, H = 0.07)
  expect_equal(bcv_criterion(4, X, soni), legacy$BCVfST(4, 0.07, X, Xsd, ni, n),
               tolerance = 1e-7)
  ipu <- function(u, x, h) qdf_indirect_poisson(u, x, h, invert = "uniroot")
  expect_equal(bcv_criterion(m, X, ipu), legacy$BCVfIP(m, X, Xsd, ni, n), tolerance = 1e-5)
  expect_equal(bcv_criterion(m, X, qdf_indirect_poisson), legacy$BCVfIP(m, X, Xsd, ni, n), tolerance = 1e-3)
})

test_that("RLCV criteria reproduce the original (_N) functions", {
  m <- 22.3
  expect_equal(rlcv_criterion(m, X, qdf_kernel),           legacy$RLCVK_N(m, X))
  expect_equal(rlcv_criterion(m, X, qdf_kernel_corrected), legacy$RLCVKC_N(m, X))
  expect_equal(rlcv_criterion(m, X, qdf_poisson),          legacy$RLCVP_N(m, X))
  expect_equal(rlcv_criterion(m, X, qdf_bernstein),        legacy$RLCVB_N(m, X))
  expect_equal(rlcv_criterion(4, X, qdf_jones),            legacy$RLCVJ_N(4, X))
  soni <- function(u, x, h) qdf_soni(u, x, h, H = 0.07)
  expect_equal(rlcv_criterion(4, X, soni), legacy$RLCVST_N(4, 0.07, X), tolerance = 1e-7)
  ipu <- function(u, x, h) qdf_indirect_poisson(u, x, h, invert = "uniroot")
  expect_equal(rlcv_criterion(m, X, ipu), legacy$RLCVIP_N(m, X), tolerance = 1e-5)
})

test_that("qdf() with selectors reproduces the original optimize() workflow", {
  bh <- 1 / optimize(legacy$BCVfB, c(10, 50), X = X, Xsd = Xsd, ni = ni, n = n)$minimum
  fit <- qdf(X, method = "bernstein", bandwidth = "bcv", u = u)
  expect_equal(fit$h, bh)
  expect_equal(fit$q, legacy$dqb(u, X, bh))

  lh <- 1 / optimize(legacy$RLCVK_N, c(10, 50), y = X)$minimum
  fit2 <- qdf(X, method = "kernel", bandwidth = "rlcv", u = u)
  expect_equal(fit2$h, lh)
  expect_equal(fit2$q, legacy$dqk(u, X, lh))
})
