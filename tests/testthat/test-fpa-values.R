# --- reproduce one replication of the Auction_g.Rmd pipeline (Gamma(10, 2), N = 5, A = 20)
set.seed(2025)
N <- 5; A <- 20; a <- 10; b0 <- 2
G <- function(x) pgamma(x, a, b0)
beta_fun <- function(x) x - sapply(x, function(xi)
  integrate(function(y) (G(y)^(N - 1)) / (G(xi)^(N - 1)), 0, xi)$value)
bids <- c(replicate(A, beta_fun(rgamma(N, a, b0))))
bids <- sort(bids); n <- length(bids); u <- (1:n) / (n + 1)
Xsd <- diff(c(0, bids)); ni <- (1:n) / n

test_that("fpa_values() reproduces the 2025 pseudo-values (Bernstein, BCV)", {
  skip_if(!exists("legacy"), "reference code not shipped in this build")
  bh <- 1 / optimize(legacy$BCVfB, c(10, 50), X = bids, Xsd = Xsd, ni = ni, n = n)$minimum
  v_legacy <- bids + u / (N - 1) * legacy$dqb(u, bids, bh)
  fit <- fpa_values(bids, n = N, method = "bernstein", bandwidth = "bcv")
  expect_equal(fit$values_raw, v_legacy)
  expect_equal(fit$values, v_legacy)          # monotone = "none"
  # spline monotonicity step, as in nonparametric_regression() (replace = "all")
  nr <- function(x, y) {
    o <- order(x); xs <- x[o]; ys <- y[o]; nn <- length(xs)
    i0 <- floor(0.10 * nn) + 1; i1 <- ceiling(0.90 * nn)
    xt <- xs[i0:i1]; yt <- ys[i0:i1]
    sf <- lm(yt ~ splines::ns(xt, df = 6))
    fv <- predict(sf, newdata = data.frame(xt = xs)); out <- numeric(nn); out[o] <- fv; out
  }
  fit_m <- fpa_values(bids, n = N, method = "bernstein", bandwidth = "bcv", monotone = "spline")
  expect_equal(fit_m$values, nr(bids, v_legacy), tolerance = 1e-8)
  expect_true(all(diff(fit_m$values) > -1e-8))
  fit_b <- fpa_values(bids, n = N, method = "bernstein", bandwidth = "bcv",
                      monotone = "spline", replace = "boundary")
  expect_equal(fit_b$values[11:90], v_legacy[11:90])
})

test_that("GPV benchmark and isotonic option work; CRRA scales the markup", {
  fit <- fpa_values(bids, n = N, method = "gpv")
  expect_true(all(fit$values > bids))
  fit_i <- fpa_values(bids, n = N, method = "gpv_bc", monotone = "isotonic")
  expect_true(all(diff(fit_i$values) >= 0))
  f0 <- fpa_values(bids, n = N, method = "bernstein", bandwidth = 0.05)
  f1 <- fpa_values(bids, n = N, method = "bernstein", bandwidth = 0.05, eta = 0.25)
  expect_equal(f1$markup, 0.75 * f0$markup)
  expect_output(print(f1), "CRRA")
})

test_that("recovered values are close to the truth on simulated IPV data", {
  set.seed(7)
  sim <- fpa_simulate(40, 5, Q = function(u) qgamma(u, 2, 10))
  fit <- fpa_values(sim$bids, method = "indirect_poisson", bandwidth = "bcv", monotone = "spline")
  tv <- as.numeric(sim$values)[fit$order]
  mid <- 20:180
  expect_lt(median(abs(fit$values[mid] - tv[mid]) / tv[mid]), 0.05)
  expect_equal(dim(fit$values_mat), c(40, 5))
  expect_equal(sort(as.numeric(fit$values_mat)), sort(fit$values))
})

# --- copula module
test_that("psi nests GPV and delta' is the derivative of delta", {
  uu <- seq(0.05, 0.95, by = 0.05); n <- 3
  for (fm in c("clayton", "frank", "gumbel")) {
    t0 <- switch(fm, clayton = 1e-4, gumbel = 1 + 1e-4, frank = 1e-4)
    expect_lt(max(abs(copula_psi(uu, fm, t0, n) - uu / (n - 1))), 1e-3, label = fm)
    th <- copula_tau_to_theta(0.5, fm)
    eps <- 1e-5
    num <- (copula_delta(uu + eps, fm, th, n) - copula_delta(uu - eps, fm, th, n)) / (2 * eps)
    expect_equal(copula_delta_prime(uu, fm, th, n), num, tolerance = 1e-5, label = fm)
    expect_equal(copula_theta_to_tau(th, fm), 0.5, tolerance = 1e-6)
  }
  expect_error(copula_psi(0.5, "clayton", n = 3), "theta")
})

test_that("copula samplers and PML fit recover theta", {
  set.seed(11)
  for (fm in c("clayton", "frank", "gumbel")) {
    th <- copula_tau_to_theta(0.5, fm)
    U <- rcopula_archimedean(2000, 3, fm, th)
    expect_true(all(U > 0 & U < 1))
    tau_hat <- cor(U[, 1], U[, 2], method = "kendall")
    expect_lt(abs(tau_hat - 0.5), 0.06, label = fm)
    ft <- copula_fit(U, fm)
    expect_lt(abs(ft$tau - 0.5), 0.06, label = fm)
  }
})

test_that("the quantile-space bid function equals the GPV integral under IPV", {
  bf <- fpa_bid_function(Q = function(u) qgamma(u, 2, 10),
                         q = function(u) 1 / dgamma(qgamma(u, 2, 10), 2, 10), n = 5)
  v <- c(0.05, 0.1, 0.2, 0.4)
  direct <- v - sapply(v, function(vi)
    integrate(function(t) pgamma(t, 2, 10)^4, 0, vi)$value / pgamma(vi, 2, 10)^4)
  expect_equal(bf$beta(v), direct, tolerance = 1e-4)
  # numerical q by default
  bf2 <- fpa_bid_function(Q = function(u) qgamma(u, 2, 10), n = 5)
  expect_equal(bf2$beta(v), direct, tolerance = 1e-3)
})

test_that("APV inversion with the true copula recovers values", {
  set.seed(3)
  Q <- function(u) 1 + 2 * u^2; q <- function(u) 4 * u
  sim <- fpa_simulate(60, 3, Q = Q, q = q, family = "gumbel", theta = 2)
  fit <- fpa_values(sim$bids, method = "bernstein", bandwidth = "bcv",
                    family = "gumbel", theta = 2)
  tv <- as.numeric(sim$values)[fit$order]
  mid <- 20:160
  expect_lt(median(abs(fit$values[mid] - tv[mid])), 0.06)
  fit2 <- fpa_values(sim$bids, method = "bernstein", bandwidth = 0.05, family = "gumbel")
  expect_true(is.finite(fit2$copula$theta))
  expect_error(fpa_values(as.numeric(sim$bids), n = 3, family = "gumbel"), "matrix")
})
