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
#' sealed-bid auctions via quantile density functions. SSRN preprint.
#' \doi{10.2139/ssrn.7071380}
#' @examples
#' set.seed(1)
#' sim <- fpa_simulate(100, 3, Q = function(u) 10 * u)      # values U[0, 10]
#' fv  <- fpa_values(sim$bids, method = "bernstein", bandwidth = "bcv")
#' r   <- 3:8
#' cbind(r, estimate = fpa_revenue(fv, r, v_s = 2),
#'       truth = fpa_revenue_model(function(u) 10 * u, n = 3, r = r, v_s = 2))
#' @export
fpa_revenue <- function(values, r, eta = NULL, v_s = 0) {
  if (inherits(values, "fpa_values")) {
    if (is.null(values$values_mat))
      stop("Expected revenue needs the auction structure: call fpa_values() with the bids as ",
           "a matrix (one row per auction, one column per bidder), not a pooled vector. ",
           "For procurement data use proc_values() and proc_cost().", call. = FALSE)
    if (is.null(eta)) eta <- values$eta
    vm <- values$values_mat
  } else {
    values <- as_bid_matrix(values, allow_vector = FALSE, arg = "values")$bids
    if (is.null(eta)) eta <- 0
    vm <- values
  }
  if (length(v_s) != 1L && length(v_s) != length(r))
    stop("`v_s` must be a scalar or have the same length as `r`.", call. = FALSE)
  v_s <- rep_len(v_s, length(r))
  n <- ncol(vm); L <- nrow(vm); N <- length(vm)
  vs <- sort(as.numeric(vm))
  e1 <- (n - 1) / (1 - eta); e2 <- (n + eta - 2) / (1 - eta)
  Fs <- seq_len(N) / N                    # F-hat at the sorted values
  w  <- vs * Fs^e2
  cw <- c(0, cumsum(w))                   # cumulative weights
  Vc <- apply(vm, 1, max)
  iV <- findInterval(Vc, vs)              # index of last value <= V_check
  FV <- iV / N
  vapply(seq_along(r), function(j) {
    rj <- r[j]; vsj <- v_s[j]
    Fr <- findInterval(rj, vs) / N
    ir <- findInterval(rj, vs, left.open = TRUE) + 1L   # first value >= r
    keep <- Vc > rj
    if (!any(keep)) return(vsj * Fr^n)
    s0 <- (cw[iV[keep] + 1L] - cw[ir]) / N              # sum over r <= v_j <= V_check
    s0[iV[keep] < ir] <- 0
    bhat <- (rj * Fr^e1 + e1 * s0) / FV[keep]^e1
    vsj * Fr^n + sum(bhat) / L
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
                                      v_s = 0, eta = 0, p_lo = 1e-4, p_hi = 0.98, U = fpa_ugrid()) {
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
#' samples of whole auctions (by default; see `resample`) (Zincenko 2024,
#' Section 4.3; Doosti 2026b) and
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
#' @section Reproducibility:
#' Fold assignment (or bootstrap resampling) uses R's random number generator;
#' call [set.seed()] beforehand to make results exactly reproducible.
#' @param B number of bootstrap replications.
#' @param eta,v_s as in [fpa_revenue()].
#' @details If the auctions have a time order and the same bidders recur,
#'   resampling auctions or bids independently understates the sampling
#'   uncertainty, in proportion to the dependence; use
#'   `resample = "blocks"` in that case.
#' @param L block length when `resample = "blocks"`: `"auto"` for the rule of
#'   [block_length_auto()], or an integer. `L = 1` is the i.i.d. auction
#'   bootstrap. Ignored by the other schemes.
#' @param resample `"auctions"` (whole auctions with replacement; the
#'   default) or `"bids"` (pooled bids with replacement). The revenue
#'   estimator uses the highest pseudo-value in each auction, and pooled
#'   resampling breaks that structure: it shifts the bootstrap distribution
#'   by an amount of the same order as the estimation error, so the
#'   intervals under-cover (about 0.65-0.80 at nominal 0.90 for the
#'   Bernstein estimator with \eqn{m = \lceil\sqrt N\rceil} in the design of
#'   Zincenko, 2024). `"bids"` is kept for comparison only.
#' @param alpha significance levels for the intervals and band.
#' @param trim discard bootstrap deviations larger than this constant
#'   (\eqn{c_b} in Zincenko's eq. 24; default no trimming).
#' @param ... arguments to [fpa_values()] (e.g. `method`, `bandwidth`,
#'   `monotone`). A data-driven `bandwidth` is selected once, on the full
#'   sample, and held fixed in every bootstrap replication, as assumed by the
#'   theory (Doosti 2026b, Theorems 6--7).
#' @return an object of class `"fpa_revenue_boot"`: a list with `r`,
#'   `estimate`, `se`, `ci` (array `r` by `c(lower, upper)` by `alpha`),
#'   `band` (same layout), the bootstrap critical values, the `B` by
#'   `length(r)` matrix of bootstrap estimates and, when the bandwidth was
#'   selected from the data, `bandwidth`: a list with the selected `h` and `m`,
#'   the `selector`, and `at_boundary`/`boundary` (see [select_bandwidth()]).
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- fpa_simulate(60, 3, Q = function(u) 10 * u)
#' rb <- fpa_revenue_boot(sim$bids, r = c(3, 5, 7), B = 50, v_s = 2,
#'                        method = "bernstein", bandwidth = "sqrt")
#' rb
#' }
#' @export
fpa_revenue_boot <- function(bids, r, B = 200, eta = 0, v_s = 0,
                             resample = c("auctions", "bids", "blocks"), L = "auto",
                             alpha = c(0.10, 0.05, 0.01),
                             trim = Inf, ...) {
  resample <- match.arg(resample)
  L_block <- L            # public argument name; `L` below is the auction count
  bids <- as_bid_matrix(bids, allow_vector = FALSE)$bids
  L <- nrow(bids); n <- ncol(bids); N <- length(bids)
  if (resample == "blocks") {
    if (identical(L_block, "auto")) L_block <- block_length_auto(bids, winner = "max")
    L_block <- as.integer(L_block)
    if (is.na(L_block) || L_block < 1L || L_block > L)
      stop("`L` must be \"auto\" or an integer between 1 and the number of auctions.",
           call. = FALSE)
  } else L_block <- 1L
  # The bandwidth is selected once, on the full sample, and held fixed in every
  # bootstrap replication (Doosti 2026b, Theorems 6-7): re-selecting it in each
  # replicate would add selector noise that the theory does not cover.
  dots <- list(...)
  fv0 <- fpa_values(bids, eta = eta, ...)
  est <- fpa_revenue(fv0, r, eta, v_s)
  bw <- NULL
  if (!is.null(fv0$qdf) && !is.null(fv0$qdf$selection)) {
    s <- fv0$qdf$selection
    bw <- list(h = fv0$qdf$h, m = 1 / fv0$qdf$h, selector = s$selector,
               at_boundary = isTRUE(s$at_boundary), boundary = s$boundary)
    dots$bandwidth <- fv0$qdf$h
  }
  fit <- function(bm) {
    fv <- do.call(fpa_values, c(list(bm, eta = eta), dots))
    fpa_revenue(fv, r, eta, v_s)
  }
  boot <- matrix(NA_real_, B, length(r))
  for (b in seq_len(B)) {
    bm <- if (resample == "bids") matrix(sample(bids, N, replace = TRUE), L, n)
          else if (resample == "blocks") bids[resample_blocks(L, L_block), , drop = FALSE]
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
                 L = L_block,
                 bandwidth = bw, call = match.call()), class = "fpa_revenue_boot")
}

#' @export
print.fpa_revenue_boot <- function(x, digits = 4, ...) {
  cat("Expected revenue with bootstrap inference (B =", x$B, ", resampling", x$resample,
      if (identical(x$resample, "blocks")) paste0(", L = ", x$L) else "", ")\n")
  if (!is.null(x$bandwidth))
    cat(sprintf("  bandwidth    : m = %s selected by %s on the full sample, fixed in all replications%s\n",
                format(signif(x$bandwidth$m, 4)), x$bandwidth$selector,
                if (isTRUE(x$bandwidth$at_boundary))
                  sprintf("\n                 (minimum at the %s end of m_range: no bandwidth selected)",
                          x$bandwidth$boundary) else ""))
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
#' @return No return value, called for its side effect of drawing the
#'   estimated expected-revenue curve with its pointwise confidence interval
#'   (dashed) and uniform confidence band (shaded). The object `x` is
#'   returned invisibly.
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

# ---- quantile-density (Stage 3) revenue estimator ----------------------------

# Stage-3 quantile density of the pseudo-values. If a data-driven bandwidth is
# minimised at an end of m_range, refit with the published range [10, 50] and
# stats::optimize() (Doosti 2026c), and warn; stop if that fit is also at an end.
# q_V enters the revenue functional only below the reserve quantile, where the
# exact m matters little, but the edge of an arbitrary range should not decide it.
.stage3_qdf <- function(v, method, bandwidth, U, ...) {
  quiet <- function(expr) withCallingHandlers(expr, warning = function(w)
    if (grepl("did not select a bandwidth", conditionMessage(w))) invokeRestart("muffleWarning"))
  fit <- quiet(qdf(v, method = method, bandwidth = bandwidth, u = U, ...))
  if (is.null(fit$selection) || !isTRUE(fit$selection$at_boundary)) return(fit)
  sel <- fit$selection
  first <- sprintf("the %s criterion for the pseudo-values is minimised at the %s end of m_range = [%s, %s] (m = %s)",
                   toupper(sel$selector), sel$boundary, format(sel$m_range[1]),
                   format(sel$m_range[2]), format(signif(sel$m, 4)))
  dots <- list(...)
  used_grid <- if (is.null(dots$grid)) 0L else dots$grid
  published <- isTRUE(all.equal(as.numeric(sel$m_range), c(10, 50))) &&
    identical(as.numeric(used_grid), 0) && is.null(dots$optimizer)
  if (!published) {
    dots$m_range <- c(10, 50); dots$grid <- 0L; dots$optimizer <- NULL
    fit2 <- quiet(do.call(qdf, c(list(v, method = method, bandwidth = bandwidth, u = U), dots)))
    if (!isTRUE(fit2$selection$at_boundary)) {
      warning("Stage 3: ", first, "; refitted with the published range m_range = [10, 50] and ",
              "stats::optimize(), which selects m = ", format(signif(fit2$selection$m, 4)), ".",
              call. = FALSE)
      return(fit2)
    }
  }
  stop("Stage 3: ", first, if (!published) ", and also with the published range [10, 50] and stats::optimize()",
       ". The criterion does not select a bandwidth for these pseudo-values; give `bandwidth` numerically.",
       call. = FALSE)
}

# Q-hat and q-hat of the pseudo-values on the quantile grid U
qd_value_functions <- function(values, method, bandwidth, U, ...) {
  v <- sort(as.numeric(values$values))
  if (method %in% c("poisson", "bernstein", "indirect_poisson") && any(v <= 0))
    stop("The '", method, "' estimator needs positive pseudo-values; choose `qdf_method = \"kernel\"`.",
         call. = FALSE)
  N   <- length(v)
  fit <- .stage3_qdf(v, method, bandwidth, U, ...)
  qv  <- pmax(fit$q, 0)
  # Q-hat: the (continuous, type 7) empirical quantile function of the
  # pseudo-values; q-hat enters only through the markup integral of beta and the
  # Myerson condition. Building Q-hat by integrating q-hat instead is fragile when
  # the pseudo-values have a long upper tail.
  Qv  <- stats::quantile(v, U, type = 7, names = FALSE)
  list(Q = function(u) stats::approx(U, Qv, xout = u, rule = 2)$y,
       q = function(u) stats::approx(U, qv, xout = u, rule = 2)$y, fit = fit)
}

#' Expected revenue and optimal reserve by the quantile-density functional
#'
#' @description
#' Estimates the seller's expected revenue \eqn{R(p)} at reserve quantile
#' \eqn{p} (reserve price \eqn{r = Q_V(p)}) and the optimal reserve from the
#' revenue functional of Doosti (2026c, Proposition 1),
#' \deqn{R(p) = v_s\,\delta(p) + \int_p^1 \beta(u; p)\,\delta'(u)\,du,}
#' where \eqn{\delta} is the diagonal section of the copula of the private
#' values (\eqn{\delta(u) = u^n} under independence) and \eqn{\beta(u;p)} is
#' the equilibrium bid in quantile space with reserve \eqn{Q_V(p)}. Writing
#' \eqn{K(w,u) = e^{\Lambda(w) - \Lambda(u)}}, \eqn{\Lambda' = 1/\psi}, the
#' proposition's \eqn{\beta(u;p) = Q_V(u) - \int_p^u K(w,u) q_V(w)\,dw} is
#' equivalently
#' \deqn{\beta(u; p) = Q_B(u) + \int_0^p K(w, u)\, q_V(w)\,dw,}
#' because without a reserve the equilibrium bid quantile function is the
#' observed bid quantile function \eqn{Q_B}. The estimator uses this second
#' form: the empirical bid quantile function carries the upper tail, and the
#' quantile density \eqn{\hat q_V} of the pseudo-values is needed only below
#' the reserve quantile, where it is estimated reliably. The optimal reserve
#' maximises \eqn{\hat R} over the quantile grid; at an interior optimum it
#' satisfies the generalised Myerson condition
#' \eqn{r^* = v_s + q_V(p^*) H(p^*; \theta)}, which reduces to
#' \eqn{r^* = v_s + (1 - F(r^*))/f(r^*)} under independence.
#' This is Stage 3 of the copula--quantile density estimator; Stages 1--2 are
#' [fpa_values()] with the same `family`.
#'
#' @param values an [fpa_values()] object; its bids, pseudo-values, copula
#'   specification (fitted \eqn{\hat\theta}, or independence) and CRRA
#'   parameter are used.
#' @param r reserve prices at which to evaluate \eqn{\hat R}; alternatively `p`.
#' @param p reserve quantiles.
#' @param v_s seller's valuation.
#' @param qdf_method,bandwidth quantile density estimator for the
#'   pseudo-values, passed to [qdf()] (default: Poisson weights, robust
#'   likelihood cross-validation, as in Doosti 2026c).
#' @param p_max largest reserve quantile searched by `fpa_optimal_reserve_qd()`.
#' @param ... further arguments to [qdf()].
#' @return `fpa_revenue_qd()`: numeric vector \eqn{\hat R}.
#'   `fpa_optimal_reserve_qd()`: a list with `r_star`, `p_star`, `revenue`,
#'   `interior` (`FALSE` when \eqn{\hat R} is maximised at the smallest reserve,
#'   i.e. no binding reserve is optimal) and the evaluated curve `p`, `R`.
#' @references
#' Doosti, H. (2026c). Optimal reserve prices in affiliated private value
#' auctions: a copula--quantile density approach. SSRN preprint.
#' \doi{10.2139/ssrn.7273869}
#' @examples
#' \donttest{
#' set.seed(1)
#' Q <- function(u) (1 - 8 * u / 9)^(-1/2)          # truncated Pareto on [1, 3]
#' th <- copula_tau_to_theta(0.25, "clayton")
#' sim <- fpa_simulate(50, 3, Q, family = "clayton", theta = th)
#' fv  <- fpa_values(sim$bids, method = "poisson", bandwidth = "rlcv", family = "clayton")
#' fpa_optimal_reserve_qd(fv, v_s = 1)[c("r_star", "p_star", "interior")]
#' fpa_optimal_reserve_model(Q, 3, family = "clayton", theta = th, v_s = 1)$r_star
#' }
#' @export
fpa_revenue_qd <- function(values, r = NULL, p = NULL, v_s = 0, qdf_method = "poisson",
                           bandwidth = "rlcv", ...) {
  if (!inherits(values, "fpa_values")) stop("`values` must be an fpa_values() object.", call. = FALSE)
  if (is.null(p) && is.null(r)) stop("Give `r` or `p`.", call. = FALSE)
  eng <- qd_revenue_engine(values, qdf_method, bandwidth, ...)
  if (is.null(p)) p <- vapply(r, function(rr) {
    k <- findInterval(rr, eng$Qv); eng$U[max(1L, min(k, length(eng$U)))] }, numeric(1))
  eng$R_at(p, v_s)
}

#' @rdname fpa_revenue_qd
#' @export
fpa_optimal_reserve_qd <- function(values, v_s = 0, qdf_method = "poisson", bandwidth = "rlcv",
                                   p_max = 0.98, ...) {
  if (!inherits(values, "fpa_values")) stop("`values` must be an fpa_values() object.", call. = FALSE)
  eng <- qd_revenue_engine(values, qdf_method, bandwidth, ...)
  pg <- eng$U[eng$U <= p_max]
  R  <- eng$R_at(pg, v_s)
  k  <- which.max(R)
  list(r_star = eng$Qv[match(pg[k], eng$U)], p_star = pg[k], revenue = R[k],
       interior = k > 1 && pg[k] > 0.005, p = pg, R = R)
}

# grid engine for the revenue functional with beta(u;p) = Q_B(u) + int_0^p K(w,u) q_V(w) dw
qd_revenue_engine <- function(values, qdf_method, bandwidth, ...) {
  if (identical(values$copula$family, "nonparametric"))
    stop("Expected revenue needs the diagonal section of a copula; refit fpa_values() with a ",
         "parametric `family` (the specification test copula_spec_test() can guide the choice).",
         call. = FALSE)
  U  <- fpa_ugrid()
  vf <- qd_value_functions(values, qdf_method, bandwidth, U, ...)
  cop <- values$copula; n <- values$n; eta <- values$eta
  QB <- stats::quantile(values$bids, U, type = 7, names = FALSE)
  psiv <- (1 - eta) * copula_psi(U, cop$family, cop$theta, n)
  lam  <- make_lambda(U, psiv)
  dv   <- copula_delta(U, cop$family, cop$theta, n)
  dpv  <- copula_delta_prime(U, cop$family, cop$theta, n)
  qv   <- if (is.null(values$.q_true)) vf$q(U) else values$.q_true(U)   # .q_true: testing hook
  cum_trap <- function(x, y) c(0, cumsum(0.5 * diff(x) * (y[-1] + y[-length(y)])))
  # int_p^1 y, accumulated from the top: forming cy[length(cy)] - cy subtracts
  # terms of size exp(-Lambda(1e-6)) (1e28 for Clayton at tau = 0.75) and
  # loses all precision in A below.
  tail_int <- function(y) { seg <- 0.5 * diff(U) * (y[-1] + y[-length(y)]); c(rev(cumsum(rev(seg))), 0) }
  c0 <- max(lam)                                     # constant cancels in K(w,u)
  S  <- cum_trap(U, exp(lam - c0) * qv)              # int_0^p e^{Lambda(w)} q(w) dw
  A  <- tail_int(exp(-(lam - c0)) * dpv)             # int_p^1 e^{-Lambda(u)} delta'(u) du
  B  <- tail_int(QB * dpv)                           # int_p^1 Q_B(u) delta'(u) du
  Qv <- vf$Q(U)
  R_at <- function(p, v_s) {
    k <- vapply(p, function(pp) which.min(abs(U - pp)), integer(1))
    v_s * dv[k] + B[k] + S[k] * A[k]
  }
  list(U = U, Qv = Qv, R_at = R_at, A = A, S = S, B = B)   # A, S, B exposed for tests
}

# Affiliated virtual value and corner test ---------------------------------------
# (EL revision: Doosti 2026c). H(p) = int_p^1 K(p,u) delta'(u) du / delta'(p),
# computed with the tail integral accumulated from the top (see qd_revenue_engine()).
.H_grid <- function(U, psiv, dpv) {
  lam <- make_lambda(U, psiv); NU <- length(U)
  c0 <- max(lam)
  # A(p) = int_p^1 e^{-Lambda(u)} delta'(u) du  (one pass, cumulative trapezoid)
  # NB: integrate from the top (reverse cumsum). Forming cy[NU] - cy, as in
  # qd_revenue_engine(), subtracts numbers of size e^{-Lambda(1e-6)} ~ 1e12
  # and loses ~1e-4 relative precision at Clayton theta >= 2.
  y <- exp(-(lam - c0)) * dpv
  seg <- 0.5 * diff(U) * (y[-1] + y[-NU])
  A <- c(rev(cumsum(rev(seg))), 0)
  exp(lam - c0) * A / dpv
}

#' Affiliated virtual value
#'
#' @description
#' The virtual value of affiliated private values in quantile space,
#' \deqn{\omega(u; \theta) = Q_V(u) - q_V(u)\, H(u; \theta), \qquad
#'   H(p) = \frac{1}{\delta'(p)} \int_p^1 K(p, u)\, \delta'(u)\, du,}
#' with \eqn{K(w, u) = e^{\Lambda(w) - \Lambda(u)}}, \eqn{\Lambda' = 1/\psi}
#' and \eqn{\delta} the diagonal section of the copula (Doosti 2026c). The
#' derivative of expected revenue is \eqn{R'(p) = \delta'(p)\{v_s - \omega(p)\}},
#' so no binding reserve is optimal (a corner) when \eqn{\omega(u) \ge v_s}
#' for all \eqn{u}; under independence \eqn{\omega} is the Myerson virtual
#' value \eqn{v - (1 - F(v))/f(v)}.
#'
#' `fpa_virtual_value_model()` evaluates \eqn{\omega} for a known model;
#' `fpa_virtual_value()` estimates it from an [fpa_values()] object, with the
#' quantile function and quantile density of the pseudo-values estimated as in
#' [fpa_revenue_qd()].
#'
#' @param Q,q quantile function and quantile density of the private values.
#' @param n number of bidders.
#' @param family,theta copula family and parameter, as in [copula_psi()].
#' @param U quantile grid (default [fpa_ugrid()]).
#' @param values an [fpa_values()] object fitted with a parametric copula.
#' @param qdf_method,bandwidth quantile density estimator for the
#'   pseudo-values, passed to [qdf()].
#' @param ... further arguments to [qdf()].
#' @return A list with the grid `U`, the virtual value `omega` and the factor
#'   `H`; `fpa_virtual_value()` also returns the estimated `Q` and `q` on `U`
#'   and the bandwidth `h` used for the pseudo-values.
#' @references
#' Doosti, H. (2026c). Optimal reserve prices in affiliated private value
#' auctions: a copula--quantile density approach. SSRN preprint.
#' \doi{10.2139/ssrn.7273869}
#' @examples
#' Q <- function(u) (1 - 8 * u / 9)^(-1/2)              # truncated Pareto on [1, 3]
#' q <- function(u) (4 / 9) * (1 - 8 * u / 9)^(-3/2)
#' vv <- fpa_virtual_value_model(Q, q, n = 3, family = "clayton", theta = 2)
#' range(vv$omega)
#' @export
fpa_virtual_value_model <- function(Q, q, n, family = "independence", theta = NULL,
                                    U = fpa_ugrid()) {
  psiv <- copula_psi(U, family, theta, n); dpv <- copula_delta_prime(U, family, theta, n)
  list(U = U, omega = Q(U) - q(U) * .H_grid(U, psiv, dpv), H = .H_grid(U, psiv, dpv))
}

#' @rdname fpa_virtual_value_model
#' @export
fpa_virtual_value <- function(values, qdf_method = "poisson", bandwidth = "rlcv", ...) {
  U <- fpa_ugrid(); n <- values$n; cop <- values$copula
  vf <- qd_value_functions(values, qdf_method, bandwidth, U, ...)
  psiv <- (1 - values$eta) * copula_psi(U, cop$family, cop$theta, n)
  dpv  <- copula_delta_prime(U, cop$family, cop$theta, n)
  H    <- .H_grid(U, psiv, dpv)
  list(U = U, omega = vf$Q(U) - vf$q(U) * H, H = H, Q = vf$Q(U), q = vf$q(U), h = vf$fit$h)
}

#' Test for a corner solution of the optimal reserve under affiliation
#'
#' @description
#' Decides whether any binding reserve price is optimal. The point rule
#' classifies the sample as a corner (no binding reserve optimal) when
#' \deqn{T = \min_{u \in [p_{lo},\, p_{hi}]} \hat\omega(u) - v_s \ge 0,}
#' with \eqn{\hat\omega} the estimated affiliated virtual value of
#' [fpa_virtual_value()]. The window starts at `p_lo = 0.10` rather than at
#' the boundary: with \eqn{v_s} equal to the lower end of the support,
#' \eqn{\omega(u) \to v_s} as \eqn{u \to 0} for every copula, so the region
#' near zero carries no information about the corner, while the quantile
#' density of the pseudo-values is least accurately estimated there.
#'
#' The decision is the point rule on the full sample (`decision = "point"`,
#' the default and the rule used in Doosti 2026c): a corner is reported when
#' \eqn{T \ge 0}. The null tested is "no reserve above the `p_lo` quantile
#' beats no reserve". In applications the copula family should be chosen with
#' [copula_spec_test()].
#'
#' With `B > 0`, auctions are resampled with replacement, with the bandwidths
#' of both stages selected once on the full sample and held fixed in every
#' replication, and the frequency with which the point rule reports a corner is
#' returned. Two alternative uses of the bootstrap are available but are not
#' the default: `decision = "bootstrap"` reports a corner when that frequency
#' is at least `level`; and `corner_not_rejected` records whether the upper
#' \eqn{1 - \alpha} bootstrap percentile of \eqn{T} is non-negative (a
#' conservative rule, since \eqn{T}, a minimum of a noisy curve, is biased
#' downward).
#'
#' @param bids an `L` by `n` matrix of bids, one row per auction.
#' @param v_s seller's valuation.
#' @param family copula family used to recover the values (see [fpa_values()]);
#'   required, because the corner result depends on it.
#' @param method,bandwidth first-stage quantile density estimator and
#'   bandwidth rule for [fpa_values()]; `bandwidth` is also used for the
#'   quantile density of the pseudo-values.
#' @param qdf_method quantile density estimator for the pseudo-values.
#' @param p_lo,p_hi window of quantiles over which the minimum is taken.
#' @param B number of bootstrap replications (0: point rule only).
#' @param decision `"point"` (default): the point rule on the full sample;
#'   `"bootstrap"`: a corner when the bootstrap frequency is at least `level`
#'   (requires `B > 0`).
#' @param level bootstrap corner frequency above which a corner is reported
#'   when `decision = "bootstrap"`.
#' @param alpha level for the upper bootstrap percentile of the statistic.
#' @param ... further arguments to [fpa_values()].
#' @param x an `"fpa_corner_test"` object.
#' @return An object of class `"fpa_corner_test"`: a list with the
#'   `statistic` \eqn{T}, the decision `corner`, the estimated virtual value
#'   `omega` on the grid `U`, `v_s`, `p_range`, the `bandwidth`s used in
#'   the two stages and the `decision` rule; with `B > 0` also
#'   `boot_stat`, `boot_corner_freq`, `boot_upper` and `corner_not_rejected`.
#'   The `print()` method is called for its side effect of printing the
#'   decision and returns `x` invisibly.
#' @references
#' Doosti, H. (2026c). Optimal reserve prices in affiliated private value
#' auctions: a copula--quantile density approach. SSRN preprint.
#' \doi{10.2139/ssrn.7273869}
#' @section Reproducibility:
#' Call `set.seed()` beforehand to make the bootstrap reproducible.
#' @examples
#' \donttest{
#' set.seed(1)
#' Q <- function(u) (1 - 8 * u / 9)^(-1/2)
#' sim <- fpa_simulate(60, 3, Q, family = "clayton",
#'                     theta = copula_tau_to_theta(0.5, "clayton"))
#' fpa_corner_test(sim$bids, v_s = 1, family = "clayton")
#' }
#' @export
fpa_corner_test <- function(bids, v_s = 0, family, method = "poisson",
                            bandwidth = "rlcv", qdf_method = "poisson",
                            p_lo = 0.10, p_hi = 0.98, B = 0, decision = c("point", "bootstrap"),
                            level = 0.5, alpha = 0.1, ...) {
  decision <- match.arg(decision)
  if (decision == "bootstrap" && B <= 0)
    stop("`decision = \"bootstrap\"` needs `B > 0` bootstrap replications.", call. = FALSE)
  if (missing(family))
    stop("`family` is required: the corner result depends on the copula family ",
         "(e.g. \"clayton\", \"gumbel\", \"frank\", \"independence\").", call. = FALSE)
  one <- function(bm, bw1, bw3) {
    fv <- fpa_values(bm, method = method, bandwidth = bw1, family = family, ...)
    vv <- fpa_virtual_value(fv, qdf_method = qdf_method, bandwidth = bw3)
    idx <- vv$U >= p_lo & vv$U <= p_hi
    list(stat = min(vv$omega[idx]) - v_s, omega = vv$omega, U = vv$U,
         h1 = if (is.null(fv$qdf)) NULL else fv$qdf$h, h3 = vv$h)
  }
  r0 <- one(bids, bandwidth, bandwidth)
  # bandwidths selected once on the full sample and held fixed in the bootstrap
  bw1 <- if (is.null(r0$h1)) bandwidth else r0$h1
  bw3 <- if (is.null(r0$h3)) bandwidth else r0$h3
  out <- list(statistic = r0$stat, corner = r0$stat >= 0, omega = r0$omega, U = r0$U,
              v_s = v_s, p_range = c(p_lo, p_hi), bandwidth = list(values = bw1, pseudo_values = bw3),
              decision = decision)
  if (B > 0) {
    L <- nrow(bids)
    bs <- vapply(seq_len(B), function(b) {
      s <- tryCatch(one(bids[sample.int(L, L, replace = TRUE), , drop = FALSE], bw1, bw3)$stat,
                    error = function(e) NA_real_)
      s }, numeric(1))
    out$boot_stat <- bs
    out$boot_corner_freq <- mean(bs >= 0, na.rm = TRUE)
    # "corner not rejected": the upper (1 - alpha) bootstrap percentile of the
    # statistic is non-negative.  The point statistic is biased downward (a
    # minimum over u of a noisy curve), so this is the conservative rule.
    out$boot_upper <- unname(stats::quantile(bs, 1 - alpha, na.rm = TRUE))
    out$corner_not_rejected <- out$boot_upper >= 0
    if (decision == "bootstrap") out$corner <- out$boot_corner_freq >= level
  }
  class(out) <- "fpa_corner_test"
  out
}

#' @rdname fpa_corner_test
#' @export
print.fpa_corner_test <- function(x, ...) {
  cat(sprintf("Corner test: min omega_hat - v_s on [%.2f, %.2f] = %.4f\n",
              x$p_range[1], x$p_range[2], x$statistic))
  if (!is.null(x$boot_corner_freq))
    cat(sprintf("Bootstrap frequency of the corner decision: %.3f (B = %d)\n",
                x$boot_corner_freq, length(x$boot_stat)))
  rule <- if (identical(x$decision, "bootstrap")) "bootstrap frequency rule" else "point rule"
  cat("Decision (", rule, "): ",
      if (x$corner) "CORNER (no binding reserve optimal)" else "INTERIOR optimum", "\n", sep = "")
  invisible(x)
}
