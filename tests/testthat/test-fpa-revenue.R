# --- legacy revenue estimator from Revenue26April.R (verbatim)
rev_est_legacy <- function(rpr, vals_mat, eta = 0, V_S = 2) {
  nb <- ncol(vals_mat); vals_all <- c(vals_mat); N <- length(vals_all)
  vals_sorted <- sort(vals_all); expo <- (nb - 1) / (1 - eta); exp2 <- (nb + eta - 2) / (1 - eta)
  cdf <- function(v) findInterval(v, vals_sorted) / N
  Fvals <- cdf(vals_sorted); Fr <- cdf(rpr); weighted_vals <- vals_sorted * Fvals^exp2
  V_check <- apply(vals_mat, 1, max); Fv_check <- cdf(V_check)
  L <- nrow(vals_mat); out <- numeric(length(rpr))
  for (j in seq_along(rpr)) {
    r_j <- rpr[j]; Fr_j <- Fr[j]; contrib <- 0
    for (l in 1:L) {
      if (V_check[l] <= r_j) next
      idx <- which(vals_sorted >= r_j & vals_sorted <= V_check[l])
      b_hat <- if (length(idx) == 0) r_j else {
        s0 <- sum(weighted_vals[idx]) / N; (r_j * Fr_j^expo + expo * s0) / Fv_check[l]^expo }
      contrib <- contrib + b_hat
    }
    out[j] <- V_S * Fr_j^nb + contrib / L
  }
  out
}

test_that("fpa_revenue reproduces the legacy plug-in estimator", {
  set.seed(1)
  vm <- matrix(runif(150, 0, 10), 50, 3)
  r <- c(3, 4, 5, 6, 7, 8)
  expect_equal(fpa_revenue(vm, r, eta = 0, v_s = 2), rev_est_legacy(r, vm, 0, 2))
  expect_equal(fpa_revenue(vm, r, eta = 0.25, v_s = 2), rev_est_legacy(r, vm, 0.25, 2))
})

test_that("model revenue matches the U[0,10] closed form and the plug-in converges to it", {
  # Zincenko (2024) design: U[0, 10], n = 3, v_s = 2
  Q <- function(u) 10 * u; n <- 3; v_s <- 2
  r <- c(3, 5, 7)
  # closed form: R(r) = v_s (r/10)^n + n \int_r^10 beta(v; r) (v/10)^(n-1) dv / 10,
  # beta(v; r) = (r^n + (n-1)/n (v^n - r^n)) / v^(n-1)
  Rtrue <- vapply(r, function(rr) {
    beta <- function(v) (rr^n + (n - 1) / n * (v^n - rr^n)) / v^(n - 1)
    v_s * (rr / 10)^n + integrate(function(v) beta(v) * n * (v / 10)^(n - 1) / 10, rr, 10)$value
  }, numeric(1))
  expect_equal(fpa_revenue_model(Q, n, r = r, v_s = v_s), Rtrue, tolerance = 1e-3)
  set.seed(2)
  sim <- fpa_simulate(2000, n, Q)
  expect_equal(fpa_revenue(sim$values, r, v_s = v_s), Rtrue, tolerance = 0.02)
  # optimal reserve for U[0,10], v_s = 2 is r* = (10 + 2)/2 = 6
  o <- fpa_optimal_reserve_model(Q, n, v_s = v_s)
  expect_true(o$interior); expect_equal(o$r_star, 6, tolerance = 0.02)
  o2 <- fpa_optimal_reserve(sim$values, v_s = v_s)
  expect_lt(abs(o2$r_star - 6), 0.5)
})

test_that("revenue from recovered values is close to the truth", {
  set.seed(3)
  Q <- function(u) 10 * u; r <- c(3, 5, 7)
  sim <- fpa_simulate(100, 3, Q)
  Rt <- fpa_revenue_model(Q, 3, r = r, v_s = 2)
  fv <- fpa_values(sim$bids, method = "bernstein", bandwidth = "wbcv")
  expect_lt(max(abs(fpa_revenue(fv, r, v_s = 2) - Rt) / Rt), 0.1)
  fz <- fpa_values(sim$bids, method = "zincenko")
  expect_lt(max(abs(fpa_revenue(fz, r, v_s = 2) - Rt) / Rt), 0.15)
  expect_error(fpa_values(sim$bids, method = "zincenko", family = "clayton", theta = 1), "independent")
})

test_that("bootstrap inference returns valid intervals and bands", {
  set.seed(4)
  sim <- fpa_simulate(40, 3, Q = function(u) 10 * u)
  rb <- fpa_revenue_boot(sim$bids, r = c(4, 6), B = 25, v_s = 2, method = "bernstein", bandwidth = 0.05)
  expect_s3_class(rb, "fpa_revenue_boot")
  expect_true(all(rb$se > 0))
  expect_true(all(rb$ci[, 1, ] <= rb$estimate & rb$ci[, 2, ] >= rb$estimate))
  expect_true(all(rb$band[, 1, ] <= rb$ci[, 1, ] + 1e-8))
  expect_equal(dim(rb$boot), c(25, 2))
  expect_output(print(rb), "bootstrap")
  rb2 <- fpa_revenue_boot(sim$bids, r = 5, B = 10, v_s = 2, resample = "auctions",
                          method = "gpv")
  expect_true(is.finite(rb2$se))
})

test_that("evaluation metrics and Monte Carlo runner work", {
  set.seed(5)
  pdf <- function(x) dgamma(x, 10, 2); cdf <- function(x) pgamma(x, 10, 2)
  sim <- fpa_simulate(20, 5, Q = function(u) qgamma(u, 10, 2))
  fv <- fpa_values(sim$bids, method = "bernstein", bandwidth = "bcv")
  ev <- fpa_evaluate_values(fv, pdf, cdf, true = sim$values)
  expect_named(ev, c("mise", "loglik", "ks", "msep", "msep_trim"))
  expect_true(all(is.finite(ev)))
  expect_equal(unname(ev["ks"]), unname(ks.test(fv$values, "pgamma", 10, 2)$statistic), tolerance = 1e-10)
  fd <- fpa_density(sim$bids, method = "kde_hh", scale = "minmax")
  ed <- fpa_evaluate_density(fd, pdf)
  expect_named(ed, c("mise", "mise_tail", "eloglik"))
  mc <- fpa_montecarlo(R = 2, L = 20, n = 5, Q = function(u) qgamma(u, 10, 2),
    estimators = list(B = function(b, v) fpa_values(b, method = "bernstein", bandwidth = "bcv"),
                      G = function(b, v) fpa_values(b, method = "gpv")),
    score = function(fit, b, v) fpa_evaluate_values(fit, pdf, cdf, v), seed = 1)
  expect_equal(nrow(mc), 4); expect_true(all(c("rep", "estimator", "mise", "msep") %in% names(mc)))
})

test_that("CRRA bids and revenue match the integral form of Zincenko (2024)", {
  eta <- 0.25; n <- 3
  e1 <- (n - 1) / (1 - eta); e2 <- (n + eta - 2) / (1 - eta)
  Fc <- function(x) pmin(pmax(x / 10, 0), 1); fc <- function(x) ifelse(x >= 0 & x <= 10, 0.1, 0)
  bid_int <- function(v) e1 * integrate(function(x) x * Fc(x)^e2 * fc(x), 0, v)$value / Fc(v)^e1
  bf <- fpa_bid_function(function(u) 10 * u, n = n, eta = eta)
  v <- c(1, 3, 5, 8)
  expect_equal(bf$beta(v), sapply(v, bid_int), tolerance = 1e-3)
  # revenue with reserve r under CRRA: beta(v; r) = (r F(r)^e1 + e1 int_r^v x F^e2 f) / F(v)^e1
  r <- 4
  b_r <- function(v) (r * Fc(r)^e1 + e1 * integrate(function(x) x * Fc(x)^e2 * fc(x), r, v)$value) / Fc(v)^e1
  Rt <- 2 * Fc(r)^n + integrate(Vectorize(function(v) b_r(v) * n * Fc(v)^(n - 1) * fc(v)), r, 10)$value
  expect_equal(fpa_revenue_model(function(u) 10 * u, n, r = r, v_s = 2, eta = eta), Rt, tolerance = 1e-3)
  set.seed(9)
  sim <- fpa_simulate(60, n, function(u) 10 * u, eta = eta)
  expect_equal(sim$eta, eta)
  fv <- fpa_values(sim$bids, method = "bernstein", bandwidth = "wbcv", eta = eta)
  expect_lt(abs(fpa_revenue(fv, r, v_s = 2) - Rt) / Rt, 0.1)
})
