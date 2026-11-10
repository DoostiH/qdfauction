# Regression tests for the changes merged in 0.2.0 (EL/JNS revision fixes)

test_that("Frank copula objects are positive and stable at large theta", {
  th <- copula_tau_to_theta(0.9, "frank")
  expect_gt(copula_psi(0.999, "frank", th, 3), 0)
  expect_true(all(is.finite(copula_delta(c(0.5, 0.999), "frank", th, 3))))
  expect_true(all(copula_delta_prime(c(0.5, 0.999), "frank", th, 3) > 0))
  # agrees with the direct formula where that formula is accurate
  old_psi <- function(u, theta, n) {
    a <- expm1(-theta * u); D <- expm1(-theta); y <- 1 + a^n / D^(n - 1)
    -y * a / ((n - 1) * theta * exp(-theta * u))
  }
  t3 <- copula_tau_to_theta(0.3, "frank"); u <- c(0.1, 0.5, 0.9)
  expect_equal(copula_psi(u, "frank", t3, 4), old_psi(u, t3, 4), tolerance = 1e-12)
})

test_that("revenue engine tail integral equals a reverse-cumsum reference (Clayton tau = 0.75)", {
  set.seed(31)
  Q <- function(u) (1 - 8 * u / 9)^(-1/2)
  th <- copula_tau_to_theta(0.75, "clayton")
  sim <- fpa_simulate(60, 3, Q, family = "clayton", theta = th)
  fv <- fpa_values(sim$bids, method = "poisson", bandwidth = 0.05,
                   family = "clayton", theta = th)
  eng <- qdfauction:::qd_revenue_engine(fv, "poisson", 0.05)
  U <- eng$U
  psiv <- copula_psi(U, "clayton", th, 3); lam <- qdfauction:::make_lambda(U, psiv)
  y <- exp(-(lam - max(lam))) * copula_delta_prime(U, "clayton", th, 3)
  seg <- 0.5 * diff(U) * (y[-1] + y[-length(y)])
  A_ref <- c(rev(cumsum(rev(seg))), 0)
  expect_equal(eng$A, A_ref, tolerance = 1e-12)
  expect_true(all(eng$A >= 0))
  expect_true(all(diff(eng$A) <= 0))
})

test_that("qdf_bernstein_loo() equals brute-force leave-one-out refits", {
  set.seed(32)
  x <- sort(rgamma(50, 4, 2)); u <- seq_along(x) / length(x); h <- 1 / 18
  for (z in c(TRUE, FALSE)) {
    bf <- sapply(seq_along(x), function(i) qdf_bernstein(u[i], x[-i], h, zero_at_origin = z))
    expect_equal(qdf_bernstein_loo(u, x, h, zero_at_origin = z), bf, tolerance = 1e-12)
  }
})

test_that("cross-validation uses exact leave-one-out, with the Bernstein shortcut", {
  set.seed(33)
  x <- rgamma(40, 4, 2)
  wrapped <- function(u, x, h) qdf_bernstein(u, x, h)     # forces brute-force refits
  for (crit in list(bcv_criterion, rlcv_criterion, wbcv_criterion))
    expect_equal(suppressWarnings(crit(20, x, qdf_bernstein, loo = TRUE)),
                 suppressWarnings(crit(20, x, wrapped, loo = TRUE)), tolerance = 1e-10)
  # the criteria can be called directly with their default u_scheme
  expect_true(is.finite(bcv_criterion(20, x, qdf_bernstein)))
  expect_true(is.finite(suppressWarnings(wbcv_criterion(20, x, qdf_bernstein))))
})

test_that("select_bandwidth(): optimize() by default, exhaustive on request; a user optimizer is honoured", {
  set.seed(34)
  x <- rgamma(50, 4, 2)
  for (sel in c("bcv", "rlcv")) {
    sb <- suppressWarnings(select_bandwidth(x, qdf_bernstein, sel, m_range = c(10, 30)))
    s0 <- suppressWarnings(select_bandwidth(x, qdf_bernstein, sel, m_range = c(10, 30), grid = 0L))
    expect_null(sb$grid_m)                         # optimize(), as in 0.1.0
    expect_equal(sb$m, s0$m)
  }
  se <- suppressWarnings(select_bandwidth(x, qdf_bernstein, "bcv", m_range = c(10, 30),
                                          grid = "exhaustive"))
  expect_equal(se$grid_m, 10:30)
  expect_equal(se$m, se$grid_m[which.min(se$grid_cv)])
  expect_equal(se$h, 1 / se$m)
  calls <- 0
  opt <- function(f, interval) { calls <<- calls + 1; list(minimum = 15, objective = f(15)) }
  s2 <- select_bandwidth(x, qdf_bernstein, "bcv", m_range = c(10, 30), optimizer = opt)
  expect_equal(calls, 1)
  expect_equal(s2$m, 15)
  expect_false(s2$at_boundary)
})

test_that("WBCV warns, once per selection, and uses optimize() by default", {
  set.seed(44)
  x <- rgamma(50, 4, 2)
  expect_warning(wbcv_criterion(20, x, qdf_bernstein), "units of the data")
  n_warn <- 0
  s <- withCallingHandlers(select_bandwidth(x, qdf_bernstein, "wbcv", m_range = c(10, 30)),
    warning = function(w) { if (grepl("units of the data", conditionMessage(w))) n_warn <<- n_warn + 1
                            invokeRestart("muffleWarning") })
  expect_equal(n_warn, 1)
  expect_null(s$grid_m)
  expect_warning(qdf(x, "bernstein", bandwidth = "wbcv"), "units of the data")
})

test_that("boundary minimisers are flagged and warned about, whatever the search", {
  set.seed(40)
  x <- rgamma(60, 4, 2)
  opt_hi <- function(f, interval) list(minimum = interval[2], objective = f(interval[2]))
  expect_warning(s1 <- select_bandwidth(x, qdf_bernstein, "bcv", m_range = c(10, 30),
                                        optimizer = opt_hi), "upper end")
  expect_true(s1$at_boundary); expect_identical(s1$boundary, "upper")
  ex <- withCallingHandlers(
    select_bandwidth(x, qdf_bernstein, "bcv", m_range = c(10, 30), grid = "exhaustive"),
    warning = function(w) invokeRestart("muffleWarning"))
  k <- which.min(ex$grid_cv)
  expect_identical(ex$at_boundary, k %in% c(1L, length(ex$grid_m)))
  expect_warning(q <- qdf(x, "bernstein", bandwidth = "bcv", m_range = c(10, 30),
                          optimizer = opt_hi), "upper end")
  expect_true(q$selection$at_boundary)
  expect_output(print(q), "no bandwidth selected")
  # EL setting on the timber bids: RLCV, Poisson, optimize() over [10, 50] is interior
  data("timber", package = "qdfauction")
  q0 <- qdf(as.numeric(timber$n3), "poisson", bandwidth = "rlcv")
  expect_false(q0$selection$at_boundary)
  expect_equal(q0$h, 0.02272727, tolerance = 1e-5)
})

test_that("fpa_optimal_reserve_model() finds interior optima below the 2nd percentile", {
  Q <- function(u) (1 - 8 * u / 9)^(-1/2)
  q <- function(u) (4 / 9) * (1 - 8 * u / 9)^(-3/2)
  o <- fpa_optimal_reserve_model(Q, 3, q, "gumbel", copula_tau_to_theta(0.9, "gumbel"), v_s = 1)
  expect_true(o$interior)
  expect_lt(o$p_star, 0.02)
})

test_that("affiliated virtual value reduces to the Myerson virtual value under independence", {
  vv <- fpa_virtual_value_model(function(u) u, function(u) rep(1, length(u)), n = 3)
  k <- vv$U > 0.01 & vv$U < 0.99
  expect_equal(vv$omega[k], 2 * vv$U[k] - 1, tolerance = 1e-4)
  set.seed(35)
  sim <- fpa_simulate(40, 3, function(u) qbeta(u, 2, 2), family = "clayton", theta = 1)
  fv <- fpa_values(sim$bids, method = "poisson", bandwidth = 0.05, family = "clayton")
  ev <- fpa_virtual_value(fv, bandwidth = 0.05)
  expect_true(all(is.finite(ev$omega)))
  expect_named(ev, c("U", "omega", "H", "Q", "q", "h"))
})

test_that("fpa_corner_test() returns a classed decision, with and without bootstrap", {
  set.seed(36)
  Q <- function(u) (1 - 8 * u / 9)^(-1/2)
  sim <- fpa_simulate(40, 3, Q, family = "clayton", theta = copula_tau_to_theta(0.5, "clayton"))
  ct <- fpa_corner_test(sim$bids, v_s = 1, family = "clayton", bandwidth = 0.05)
  expect_s3_class(ct, "fpa_corner_test")
  expect_type(ct$corner, "logical")
  expect_output(print(ct), "Decision")
  ctb <- fpa_corner_test(sim$bids, v_s = 1, family = "clayton", bandwidth = 0.05, B = 3)
  expect_length(ctb$boot_stat, 3)
  expect_true(ctb$boot_corner_freq >= 0 && ctb$boot_corner_freq <= 1)
})

test_that("fpa_montecarlo() passes p0 and eta to the simulator", {
  set.seed(37)
  mc <- fpa_montecarlo(R = 2, L = 10, n = 3, Q = function(u) qbeta(u, 2, 2), eta = 0.3,
                       estimators = list(gpv = function(b, v) fpa_values(b, method = "gpv")),
                       score = function(fit, b, v) c(max_bid_share = max(b) / max(v)))
  expect_equal(nrow(mc), 2)
  expect_true(all(mc$max_bid_share < 1))
})

# --- wavelet functions without tests so far -----------------------------------

test_that("psi_jk and the periodised basis are orthonormal on [0, 1]", {
  z <- seq(0, 1, length.out = 2^12 + 1); w <- c(0.5, rep(1, 2^12 - 1), 0.5) / 2^12
  P <- cbind(phi_per(z, 4, 0:15), psi_per(z, 4, 0:15))
  G <- crossprod(P * w, P)
  expect_equal(G, diag(32), tolerance = 2e-2, ignore_attr = TRUE)
  v <- psi_jk(seq(-2, 7, by = 0.01), 0, 0)
  expect_equal(sum(v) * 0.01, 0, tolerance = 1e-3)            # zero mean
  expect_equal(sum(v^2) * 0.01, 1, tolerance = 1e-2)          # unit norm
})

test_that("loclin_smooth() reproduces linear functions", {
  g <- seq(0, 1, length.out = 101)
  expect_equal(loclin_smooth(c(0.2, 0.5, 0.9), g, 3 * g + 1, h = 0.1),
               3 * c(0.2, 0.5, 0.9) + 1, tolerance = 1e-8)
})

test_that("select_wavelet() returns tuning parameters from its grids", {
  set.seed(38)
  x <- rgamma(150, 4, 2)
  s <- select_wavelet(x, j0_grid = 3:4, h_grid = c(0.1, 0.2), folds = 3)
  expect_true(s$j0 %in% 3:4)
  expect_true(s$h %in% c(0.1, 0.2))
  expect_identical(s$criterion, "lscv")
})

test_that("exact LOO and exhaustive search do not affect the wavelet estimators", {
  set.seed(39)
  x <- rgamma(120, 4, 2)
  a <- qdf(x, "wavelet", bandwidth = 0.15, j0 = 4)
  b <- qdf(x, "wavelet", bandwidth = 0.15, j0 = 4, grid = 0L, loo = FALSE)
  expect_equal(a$q, b$q)
})

# --- 0.2.0 follow-up: requested changes ------------------------------------------

test_that("wbcv_criterion() uses exact leave-one-out by default", {
  set.seed(41)
  x <- rgamma(40, 4, 2)
  w <- function(...) suppressWarnings(wbcv_criterion(...))
  expect_equal(w(20, x, qdf_bernstein), w(20, x, qdf_bernstein, loo = TRUE))
  expect_false(isTRUE(all.equal(w(20, x, qdf_bernstein), w(20, x, qdf_bernstein, loo = FALSE))))
})

test_that("select_bandwidth() rejects invalid grid values", {
  x <- rgamma(40, 4, 2)
  expect_error(select_bandwidth(x, qdf_bernstein, "bcv", grid = "exhaustve"), "exhaustive")
  expect_error(select_bandwidth(x, qdf_bernstein, "bcv", grid = -1), "non-negative")
  expect_error(select_bandwidth(x, qdf_bernstein, "bcv", grid = 2.5), "integer")
  expect_error(qdf(x, "bernstein", bandwidth = "bcv", grid = TRUE), "grid")
  expect_silent(suppressWarnings(select_bandwidth(x, qdf_bernstein, "bcv", grid = 5)))
})

test_that("qdf_bernstein_loo() checks its inputs", {
  x <- rgamma(30, 4, 2); u <- seq_along(x) / (length(x) + 1)
  expect_error(qdf_bernstein_loo(u[-1], x, 0.1), "one value per observation")
  expect_error(qdf_bernstein_loo(rev(u), x, 0.1), "increasing")
  expect_error(qdf_bernstein_loo(u, x, 0.7), "floor")
  expect_error(qdf_bernstein_loo(u, c(-1, x[-1]), 0.1), "positive")
})

test_that("fpa_corner_test() requires the copula family", {
  set.seed(42)
  sim <- fpa_simulate(30, 3, function(u) 1 + u, family = "clayton", theta = 1)
  expect_error(fpa_corner_test(sim$bids, v_s = 1), "family")
})

test_that("fpa_revenue_boot() selects the bandwidth once and reuses it", {
  set.seed(43)
  sim <- fpa_simulate(40, 3, Q = function(u) 10 * u)
  n_warn <- 0
  rb <- withCallingHandlers(
    fpa_revenue_boot(sim$bids, r = c(3, 5), B = 5, v_s = 2, method = "bernstein",
                     bandwidth = "bcv"),
    warning = function(w) { n_warn <<- n_warn + 1; invokeRestart("muffleWarning") })
  expect_false(is.null(rb$bandwidth))
  expect_equal(rb$bandwidth$m, 1 / rb$bandwidth$h)
  expect_lte(n_warn, 1)                                    # at most the full-sample fit
  expect_identical(n_warn == 1, rb$bandwidth$at_boundary)
  expect_output(print(rb), "fixed in all replications")
  # with a numeric bandwidth nothing is selected
  rb2 <- fpa_revenue_boot(sim$bids, r = c(3, 5), B = 3, v_s = 2, method = "bernstein",
                          bandwidth = 0.05)
  expect_null(rb2$bandwidth)
})

test_that("Stage 3 falls back to [10, 50] with optimize(), or stops, at a boundary", {
  set.seed(1)                                             # small sample: N = 105 bids
  sim <- fpa_simulate(35, 3, function(u) 1 + 10 * qbeta(u, 2, 2), family = "gumbel", theta = 1.5)
  fv <- fpa_values(sim$bids, method = "poisson", bandwidth = 0.05, family = "gumbel")
  # default (optimize over [10, 50]) is interior: no warning
  expect_silent(r0 <- fpa_revenue_qd(fv, p = 0.2, v_s = 1))
  # a search that ends at the edge of a non-published range: refit and warn
  edge <- function(f, interval) list(minimum = interval[2], objective = f(interval[2]))
  expect_warning(r1 <- fpa_revenue_qd(fv, p = 0.2, v_s = 1, optimizer = edge,
                                      m_range = c(10, 100)), "refitted with the published range")
  expect_equal(r1, r0)
  # an edge solution although the published range and optimize() were used: stop
  expect_error(withCallingHandlers(
    qdfauction:::.stage3_qdf(sort(fv$values), "poisson", "rlcv", fpa_ugrid(),
                             m_range = c(10, 50), grid = 0L, optimizer = NULL),
    error = function(e) stop(e)), NA)                      # interior here: no error

})

test_that("Stage 3 stops when the published search is also at a boundary", {
  fake <- function(x, method, bandwidth, u, m_range = c(10, 50), grid = 0L, ...)
    list(q = rep(1, length(u)), h = 1 / m_range[2],
         selection = list(selector = "rlcv", m = m_range[2], m_range = m_range,
                          at_boundary = TRUE, boundary = "upper"))
  ## `.package` is required when the tests run against the installed package
  ## (as under R CMD check) rather than a pkgload-loaded source tree.
  local_mocked_bindings(qdf = fake, .package = "qdfauction")
  expect_error(.stage3_qdf(1:50, "poisson", "rlcv", fpa_ugrid()),
               "does not select a bandwidth")
  expect_error(.stage3_qdf(1:50, "poisson", "rlcv", fpa_ugrid(), m_range = c(10, 100)),
               "also with the published range")
})

test_that("bandwidth = \"sqrt\" gives m = ceiling(sqrt(N)) for the Bernstein estimator", {
  data("timber", package = "qdfauction")
  q2 <- qdf(as.numeric(timber$n2), "bernstein", bandwidth = "sqrt")
  q3 <- qdf(as.numeric(timber$n3), "bernstein", bandwidth = "sqrt")
  expect_equal(c(q2$m, q3$m), c(15, 18))
  expect_null(q3$selection)
  expect_output(print(q3), "ceiling\\(sqrt\\(n\\)\\)")
  expect_error(qdf(as.numeric(timber$n3), "poisson", bandwidth = "sqrt"), "bernstein")
  fv <- fpa_values(timber$n3, method = "bernstein", bandwidth = "sqrt")
  expect_equal(fv$qdf$m, 18)
  set.seed(45)
  rb <- fpa_revenue_boot(timber$n3, r = c(40, 60), B = 3, v_s = 20,
                         method = "bernstein", bandwidth = "sqrt")
  expect_null(rb$bandwidth)                                # a rule, nothing selected
})

test_that("fpa_corner_test(): the point rule is the default decision", {
  set.seed(46)
  Q <- function(u) (1 - 8 * u / 9)^(-1/2)
  sim <- fpa_simulate(40, 3, Q, family = "clayton", theta = copula_tau_to_theta(0.5, "clayton"))
  ct <- fpa_corner_test(sim$bids, v_s = 1, family = "clayton", bandwidth = 0.05, B = 4)
  expect_identical(ct$decision, "point")
  expect_identical(ct$corner, ct$statistic >= 0)
  expect_output(print(ct), "point rule")
  ctb <- fpa_corner_test(sim$bids, v_s = 1, family = "clayton", bandwidth = 0.05, B = 4,
                         decision = "bootstrap")
  expect_identical(ctb$corner, ctb$boot_corner_freq >= 0.5)
  expect_error(fpa_corner_test(sim$bids, v_s = 1, family = "clayton", decision = "bootstrap"), "B > 0")
})

test_that("fpa_revenue_boot() resamples whole auctions by default", {
  expect_identical(eval(formals(fpa_revenue_boot)$resample)[1], "auctions")
  set.seed(47)
  sim <- fpa_simulate(30, 3, Q = function(u) 10 * u)
  rb <- fpa_revenue_boot(sim$bids, r = c(3, 5), B = 3, v_s = 2, method = "bernstein",
                         bandwidth = "sqrt")
  expect_identical(rb$resample, "auctions")
  expect_output(print(rb), "resampling auctions")
})
