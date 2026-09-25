# Expected revenue and optimal reserve price -----------------------------------

#' Plug-in estimator of the seller's expected revenue
#'
#' @description
#' Estimates the seller's expected revenue \eqn{R(r)} in a symmetric
#' first-price auction with reserve price \eqn{r}, from recovered private
#' values (Zincenko 2024, eq. 6; Doosti 2026b, eq. 11). With \eqn{L} auctions
#' of \eqn{n} bidders, pooled pseudo-values \eqn{\hat v_1 \le \ldots \le \hat v_N},
#' their empirical distribution \eqn{\hat F}, and the winner's value
#' \eqn{\check V_\ell = \max_i \hat v_{i\ell}} in auction \eqn{\ell},
#' \deqn{\hat R(r) = v_s \hat F(r)^n + \frac1L\sum_{\ell=1}^L 1\{\check V_\ell > r\}\,\hat\beta(\check V_\ell; r),}
#' \deqn{\hat\beta(v; r) = \frac{r\,\hat F(r)^{e_1} + e_1 N^{-1}\sum_{r \le \hat v_j \le v} \hat v_j \hat F(\hat v_j)^{e_2}}{\hat F(v)^{e_1}},}
#' with \eqn{e_1 = (n-1)/(1-\eta)}, \eqn{e_2 = (n+\eta-2)/(1-\eta)} for CRRA
#' bidders with parameter \eqn{\eta} (\eqn{\eta = 0}: risk neutral), and
#' \eqn{v_s} the seller's own valuation.
#'
#' @param values an [fpa_values()] object fitted to a bid *matrix*, or an
#'   `L` by `n` matrix of private values.
#' @param r reserve prices at which to evaluate the revenue.
#' @param eta CRRA parameter; taken from `values` when it is an
#'   [fpa_values()] object.
#' @param v_s seller's valuation of the object (revenue when unsold).
#' @return numeric vector \eqn{\hat R(r)}.
#' @references
#' Zincenko, F. (2024). Estimation and inference of seller's expected revenue
#' in first-price auctions. *Journal of Econometrics*, 241(1), 105734.
#' \doi{10.1016/j.jeconom.2024.105734}
#'
#' Doosti, H. (2026b). Estimating seller's expected revenue in first-price
#' sealed-bid auctions via quantile density functions. Working paper, SSRN
#' 7071380.
#' @examples
#' set.seed(1)
#' sim <- fpa_simulate(100, 3, Q = function(u) 10 * u)      # values U[0, 10]
#' fv  <- fpa_values(sim$bids, method = "bernstein", bandwidth = "wbcv")
#' r   <- 3:8
#' cbind(r, estimate = fpa_revenue(fv, r, v_s = 2),
#'       truth = fpa_revenue_model(function(u) 10 * u, n = 3, r = r, v_s = 2))
#' @export
fpa_revenue <- function(values, r, eta = NULL, v_s = 0) {
  if (inherits(values, "fpa_values")) {
    if (is.null(values$values_mat)) stop("`values` must come from fpa_values() on a bid matrix.", call. = FALSE)
    if (is.null(eta)) eta <- values$eta
    vm <- values$values_mat
  } else {
    if (!is.matrix(values)) stop("`values` must be a matrix or an fpa_values object.", call. = FALSE)
    if (is.null(eta)) eta <- 0
    vm <- values
  }
  n <- ncol(vm); L <- nrow(vm); N <- length(vm)
  vs <- sort(as.numeric(vm))
  e1 <- (n - 1) / (1 - eta); e2 <- (n + eta - 2) / (1 - eta)
  Fs <- seq_len(N) / N                    # F-hat at the sorted values
  w  <- vs * Fs^e2
  cw <- c(0, cumsum(w))                   # cumulative weights
  Vc <- apply(vm, 1, max)
  iV <- findInterval(Vc, vs)              # index of last value <= V_check
  FV <- iV / N
  vapply(r, function(rj) {
    Fr <- findInterval(rj, vs) / N
    ir <- findInterval(rj, vs, left.open = TRUE) + 1L   # first value >= r
    keep <- Vc > rj
    if (!any(keep)) return(v_s * Fr^n)
    s0 <- (cw[iV[keep] + 1L] - cw[ir]) / N              # sum over r <= v_j <= V_check
    s0[iV[keep] < ir] <- 0
    bhat <- (rj * Fr^e1 + e1 * s0) / FV[keep]^e1
    v_s * Fr^n + sum(bhat) / L
  }, numeric(1))
}

#' Expected revenue and optimal reserve price of a known model
#'
#' @description
#' Computes the seller's expected revenue \eqn{R(r) = v_s\,\delta(p) +
#' \int_p^1 \beta(u; p)\,d\delta(u)} for a private-value model with marginal
#' quantile function `Q`, copula multiplier \eqn{\psi} and diagonal
#' \eqn{\delta} ([copula_psi()], [copula_delta()]), where
#' \eqn{\beta(\cdot; p)} is the equilibrium bid with reserve price
#' \eqn{r = Q(p)} ([fpa_bid_function()]). Under independence and risk
#' neutrality this is the revenue used as the truth in Zincenko (2024).
#' `fpa_optimal_reserve_model()` finds the interior maximiser \eqn{p^*} of
#' \eqn{R} from the first-order condition
#' \eqn{\delta'(p)(v_s - Q(p)) + q(p)\int_p^1 e^{\Lambda(p)-\Lambda(u)}\delta'(u)\,du = 0}.
#'
#' @inheritParams fpa_bid_function
#' @param r reserve prices (on the value scale); alternatively give `p`.
#' @param p reserve quantiles.
#' @param v_s seller's valuation.
#' @param eta CRRA risk-aversion parameter of the bidders.
#' @return `fpa_revenue_model()`: numeric vector of expected revenues.
#'   `fpa_optimal_reserve_model()`: a list with `p_star`, `r_star`,
#'   `revenue` and `interior`.
#' @examples
#' Q <- function(u) (1 - 8 * u / 9)^(-1/2)             # truncated Pareto on [1, 3]
#' fpa_revenue_model(Q, n = 3, p = c(0.1, 0.3, 0.5), v_s = 1)
#' fpa_optimal_reserve_model(Q, n = 3, v_s = 1)
#' fpa_optimal_reserve_model(Q, n = 3, v_s = 1, family = "clayton", theta = 2)
#' @export
fpa_revenue_model <- function(Q, n, r = NULL, p = NULL, q = NULL, family = "independence",
                              theta = NULL, v_s = 0, eta = 0, U = fpa_ugrid()) {
  if (is.null(p) && is.null(r)) stop("Give `r` or `p`.", call. = FALSE)
  Qv <- Q(U)
  if (is.null(p)) p <- vapply(r, function(rr) {
    if (rr <= Qv[1]) return(U[1])
    if (rr >= Qv[length(Qv)]) return(U[length(U)])
    stats::approx(Qv, U, xout = rr, ties = "ordered")$y
  }, numeric(1))
  if (is.null(q)) {
    eps <- 1e-6
    q <- function(u) (Q(pmin(u + eps, 1 - 1e-9)) - Q(pmax(u - eps, 1e-9))) /
      (pmin(u + eps, 1 - 1e-9) - pmax(u - eps, 1e-9))
  }
  psiv <- (1 - eta) * copula_psi(U, family, theta, n); lam <- make_lambda(U, psiv)
  qv <- q(U); dv <- copula_delta(U, family, theta, n); dpv <- copula_delta_prime(U, family, theta, n)
  NU <- length(U)
  vapply(p, function(pp) {
    j0 <- which.min(abs(U - pp))
    m  <- markup_grid(U, qv, lam, j0)
    idx <- j0:NU
    v_s * dv[j0] + trapz(U[idx], (Qv[idx] - m[idx]) * dpv[idx])
  }, numeric(1))
}

#' @rdname fpa_revenue_model
#' @param p_lo,p_hi search interval for the optimal reserve quantile.
#' @export
fpa_optimal_reserve_model <- function(Q, n, q = NULL, family = "independence", theta = NULL,
                                      v_s = 0, eta = 0, p_lo = 0.02, p_hi = 0.98, U = fpa_ugrid()) {
  if (is.null(q)) {
    eps <- 1e-6
    q <- function(u) (Q(pmin(u + eps, 1 - 1e-9)) - Q(pmax(u - eps, 1e-9))) /
      (pmin(u + eps, 1 - 1e-9) - pmax(u - eps, 1e-9))
  }
  Qv <- Q(U); qv <- q(U)
  psiv <- (1 - eta) * copula_psi(U, family, theta, n); lam <- make_lambda(U, psiv)
  dpv <- copula_delta_prime(U, family, theta, n); NU <- length(U)
  jset <- which(U >= p_lo & U <= p_hi)
  Phi <- vapply(jset, function(j0) {
    idx <- j0:NU; Kv <- exp(lam[j0] - lam[idx])
    dpv[j0] * (v_s - Qv[j0]) + qv[j0] * trapz(U[idx], Kv * dpv[idx])
  }, numeric(1))
  sc <- which(Phi[-length(Phi)] > 0 & Phi[-1] <= 0)
  if (!length(sc)) return(list(p_star = NA_real_, r_star = NA_real_, revenue = NA_real_, interior = FALSE))
  k <- max(sc); w <- Phi[k] / (Phi[k] - Phi[k + 1])
  ps <- U[jset[k]] + w * (U[jset[k + 1]] - U[jset[k]])
  rs <- stats::approx(U, Qv, xout = ps)$y
  list(p_star = ps, r_star = rs,
       revenue = fpa_revenue_model(Q, n, p = ps, q = q, family = family, theta = theta, v_s = v_s,
                                   eta = eta, U = U),
       interior = TRUE)
}

#' Optimal reserve price from recovered private values
#'
#' Maximises the plug-in revenue [fpa_revenue()] over a grid of reserve prices.
#' @inheritParams fpa_revenue
#' @param grid candidate reserve prices; defaults to 200 points between the
#'   2nd and 98th percentiles of the values.
#' @return a list with `r_star`, `revenue` and the evaluated `grid` and `R`.
#' @export
fpa_optimal_reserve <- function(values, eta = NULL, v_s = 0, grid = NULL) {
  vm <- if (inherits(values, "fpa_values")) values$values_mat else values
  if (is.null(grid)) grid <- seq(stats::quantile(vm, 0.02), stats::quantile(vm, 0.98), length.out = 200)
  R <- fpa_revenue(values, grid, eta, v_s)
  k <- which.max(R)
  list(r_star = grid[k], revenue = R[k], grid = grid, R = R)
}

#' Bootstrap inference for the expected-revenue estimator
#'
#' @description
#' Recomputes the pseudo-values and the plug-in revenue on `B` bootstrap
#' samples of the pooled bids (Zincenko 2024, Section 4.3; Doosti 2026b) and
#' returns the bootstrap standard error
#' \eqn{\hat\sigma^*(r) = \{B^{-1}\sum_b (\hat R^*_b(r) - \hat R(r))^2\}^{1/2}},
#' symmetric percentile-\eqn{t} pointwise confidence intervals (critical value:
#' the \eqn{1-\alpha} quantile of \eqn{|Z^*_b(r)|}), and a uniform confidence
#' band over `r` based on the bootstrap distribution of
#' \eqn{\sup_r |Z^*_b(r)|}, \eqn{Z^*_b(r) = \{\hat R^*_b(r) - \hat R(r)\}/\hat\sigma^*(r)}.
#' (The original simulation code used the \eqn{1-\alpha/2} quantile of
#' \eqn{Z^*} for the pointwise intervals; the symmetric version guarantees
#' that the band contains the intervals.)
#'
#' @param bids `L` by `n` bid matrix.
#' @param r reserve prices.
#' @param B number of bootstrap replications.
#' @param eta,v_s as in [fpa_revenue()].
#' @param resample `"bids"` (pooled bids with replacement, as in the papers)
#'   or `"auctions"` (whole auctions with replacement).
#' @param alpha significance levels for the intervals and band.
#' @param trim discard bootstrap deviations larger than this constant
#'   (\eqn{c_b} in Zincenko's eq. 24; default no trimming).
#' @param ... arguments to [fpa_values()] (e.g. `method`, `bandwidth`,
#'   `monotone`).
#' @return an object of class `"fpa_revenue_boot"`: a list with `r`,
#'   `estimate`, `se`, `ci` (array `r` by `c(lower, upper)` by `alpha`),
#'   `band` (same layout), the bootstrap critical values and the `B` by
#'   `length(r)` matrix of bootstrap estimates.
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- fpa_simulate(60, 3, Q = function(u) 10 * u)
#' rb <- fpa_revenue_boot(sim$bids, r = c(3, 5, 7), B = 50, v_s = 2,
#'                        method = "bernstein", bandwidth = "wbcv")
#' rb
#' }
#' @export
fpa_revenue_boot <- function(bids, r, B = 200, eta = 0, v_s = 0,
                             resample = c("bids", "auctions"), alpha = c(0.10, 0.05, 0.01),
                             trim = Inf, ...) {
  resample <- match.arg(resample)
  if (!is.matrix(bids)) stop("`bids` must be an L x n matrix.", call. = FALSE)
  L <- nrow(bids); n <- ncol(bids); N <- length(bids)
  fit <- function(bm) {
    fv <- fpa_values(bm, eta = eta, ...)
    fpa_revenue(fv, r, eta, v_s)
  }
  est <- fit(bids)
  boot <- matrix(NA_real_, B, length(r))
  for (b in seq_len(B)) {
    bm <- if (resample == "bids") matrix(sample(bids, N, replace = TRUE), L, n)
          else bids[sample.int(L, L, replace = TRUE), , drop = FALSE]
    boot[b, ] <- tryCatch(fit(bm), error = function(e) rep(NA_real_, length(r)))
  }
  d <- sweep(boot, 2, est)
  d[abs(d) > trim] <- NA
  se <- sqrt(colMeans(d^2, na.rm = TRUE))
  Z <- sweep(d, 2, se, "/")
  ci <- array(NA_real_, c(length(r), 2, length(alpha)),
              dimnames = list(r = r, c("lower", "upper"), alpha = alpha))
  band <- ci
  sup <- apply(abs(Z), 1, function(z) if (all(is.na(z))) NA else max(z, na.rm = TRUE))
  z_pt <- matrix(NA_real_, length(r), length(alpha)); z_band <- numeric(length(alpha))
  for (a in seq_along(alpha)) {
    z_pt[, a] <- apply(abs(Z), 2, stats::quantile, probs = 1 - alpha[a], na.rm = TRUE)
    z_band[a] <- stats::quantile(sup, 1 - alpha[a], na.rm = TRUE)
    ci[, 1, a] <- est - z_pt[, a] * se; ci[, 2, a] <- est + z_pt[, a] * se
    band[, 1, a] <- est - z_band[a] * se; band[, 2, a] <- est + z_band[a] * se
  }
  structure(list(r = r, estimate = est, se = se, ci = ci, band = band, z_pointwise = z_pt,
                 z_band = z_band, alpha = alpha, boot = boot, B = B, resample = resample,
                 call = match.call()), class = "fpa_revenue_boot")
}

#' @export
print.fpa_revenue_boot <- function(x, digits = 4, ...) {
  cat("Expected revenue with bootstrap inference (B =", x$B, ", resampling", x$resample, ")\n")
  a <- x$alpha[which.min(abs(x$alpha - 0.05))]; k <- which(x$alpha == a)
  tab <- cbind(r = x$r, estimate = x$estimate, se = x$se,
               ci_lower = x$ci[, 1, k], ci_upper = x$ci[, 2, k],
               band_lower = x$band[, 1, k], band_upper = x$band[, 2, k])
  cat(sprintf("  pointwise CI and uniform band at level %g\n", 1 - a))
  print(round(tab, digits))
  invisible(x)
}

#' Plot a revenue curve with bootstrap bands
#' @param x an `"fpa_revenue_boot"` object.
#' @param alpha which level to draw.
#' @param ... passed to [graphics::plot()].
#' @export
plot.fpa_revenue_boot <- function(x, alpha = 0.05, ...) {
  k <- which.min(abs(x$alpha - alpha))
  args <- list(x = x$r, y = x$estimate, type = "l", lwd = 2, xlab = "reserve price",
               ylab = "expected revenue", ylim = range(x$band[, , k], na.rm = TRUE))
  args <- utils::modifyList(args, list(...))
  do.call(graphics::plot, args)
  graphics::polygon(c(x$r, rev(x$r)), c(x$band[, 1, k], rev(x$band[, 2, k])),
                    col = grDevices::adjustcolor("grey", 0.4), border = NA)
  graphics::lines(x$r, x$ci[, 1, k], lty = 2); graphics::lines(x$r, x$ci[, 2, k], lty = 2)
  graphics::lines(x$r, x$estimate, lwd = 2)
  invisible(x)
}

# Zincenko (2024) kernel-BC pseudo-values with linear extrapolation near the boundaries.
# Adapted from code provided by Federico Zincenko, (C) Federico Zincenko;
# see inst/COPYRIGHTS.
values_zincenko <- function(b, n, eta = 0, hg_co = 1.06, hgder_co = 2.83, tau = 0.99) {
  b <- as.numeric(b); N <- length(b)
  hg <- hg_co * stats::sd(b) * N^(-1 / 5); hgder <- hgder_co * stats::sd(b) * N^(-1 / 7)
  sxi <- hg^tau; minB <- min(b); maxB <- max(b)
  tw  <- function(z) 35 / 32 * (1 - z^2)^3 * (abs(z) < 1)
  twI <- function(z) ifelse(z <= -1, 0, ifelse(z >= 1, 1, (-5 * z^7 + 21 * z^5 - 35 * z^3 + 35 * z + 16) / 32))
  twD <- function(z) -105 / 16 * z * (1 - z^2)^2 * (abs(z) < 1)
  g   <- function(x) rowMeans(tw(outer(x, b, "-") / hg)) / hg
  G   <- function(x) rowMeans(twI(outer(x, b, "-") / hg))
  gd  <- function(x) rowMeans(twD(outer(x, b, "-") / hgder)) / hgder^2
  xi  <- function(x) x + (1 - eta) * G(x) / ((n - 1) * g(x))
  xid <- function(x) (n - eta) / (n - 1) - (1 - eta) * G(x) * gd(x) / ((n - 1) * g(x)^2)
  lo <- minB + sxi; hi <- maxB - sxi
  out <- numeric(N)
  i1 <- b < lo; i3 <- b > hi; i2 <- !(i1 | i3)
  if (any(i2)) out[i2] <- xi(b[i2])
  if (any(i1)) out[i1] <- xi(lo) + xid(lo) * (b[i1] - lo)
  if (any(i3)) out[i3] <- xi(hi) + xid(hi) * (b[i3] - hi)
  pmax(out, minB)
}
