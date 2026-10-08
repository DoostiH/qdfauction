set.seed(2012)
n <- 50; z <- rexp(n); h <- 0.15
u <- seq(1 / (n + 1), n / (n + 1), 1 / (n + 1))

test_that("qdf_soni(form = 'riemann', triangular) reproduces Soni's dentquant()", {
  expect_equal(qdf_soni(u, z, h = h, H = h, kernel = "triangular", form = "riemann"),
               soni$dentquant(z, n, h))
})

test_that("qdf_jones(triangular) reproduces Soni's jon1() (at u = (i-1)/n)", {
  expect_equal(qdf_jones((seq_len(n) - 1) / n + 1e-9, z, h = h, kernel = "triangular"),
               soni$jon1(z, n, h))
})

test_that("qdf_kernel(triangular) equals Soni's jon2() up to the omitted i = 0 term", {
  full <- qdf_kernel(u, z, h = h, kernel = "triangular")
  K <- function(a) pmax(1 - abs(a) / h, 0) / h
  term0 <- sort(z)[1] * (K(-u) - K(1 / n - u))
  expect_equal(full - term0, soni$jon2(z, n, h))
})

test_that("Epanechnikov kernel and censoring work", {
  q <- qdf_soni(u, z, h = h, H = h, kernel = "epanechnikov", form = "riemann")
  expect_true(all(is.finite(q)))
  st <- rbinom(n, 1, 0.8)
  qc <- qdf_soni(u, z, h = h, H = h, kernel = "triangular", status = st)
  expect_true(all(is.finite(qc)))
  expect_equal(qdf_soni(u, z, h = h, H = h, kernel = "triangular", status = rep(1, n)),
               qdf_soni(u, z, h = h, H = h, kernel = "triangular"))
  fit <- qdf(z, "soni", bandwidth = "bcv", H = 0.1, kernel = "epanechnikov", u = c(.3, .6))
  expect_equal(fit$kernel, "epanechnikov")
  expect_equal(predict(fit, 0.3), fit$q[1])
})

test_that("kernel CDFs integrate their densities", {
  for (k in c("triangular", "epanechnikov")) {
    K <- qdfauction:::kernel_fun(k)
    for (a in c(-0.7, 0, 0.4)) expect_equal(K$p(a), integrate(K$d, -1, a)$value, tolerance = 1e-5)
  }
})
