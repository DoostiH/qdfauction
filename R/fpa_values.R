# Private-value recovery ------------------------------------------------------

#' Recover private values from first-price auction bids
#'
#' @description
#' Nonparametric inversion of the equilibrium bid function of a symmetric
#' first-price sealed-bid auction through the quantile density of bids:
#' \deqn{\hat v = b + \lambda^{-1}\{\psi(\hat u)\,\hat q(\hat u)\},\qquad \hat u = \hat G(b),}
#' where \eqn{\hat q} is a quantile density estimate ([qdf()]), \eqn{\psi} is
#' the copula multiplier (equal to \eqn{u/(n-1)} under independent private
#' values, so that \eqn{\psi q = G/\{(n-1)g\}} and the formula is that of
#' Guerre, Perrigne and Vuong 2000), and \eqn{\lambda^{-1}} accounts for risk
#' aversion (Guerre, Perrigne and Vuong 2009): the identity for risk-neutral
#' bidders and \eqn{y \mapsto (1-\eta) y} for CRRA utility \eqn{U(t) = t^{1-\eta}}.
#'
#' `method = "zincenko"` gives the kernel estimator of Zincenko (2024): a
#' triweight kernel for \eqn{\hat G} and \eqn{\hat g} with normal-reference
#' bandwidths, and linear extrapolation of the inverse bid function within
#' \eqn{h_g^{0.99}} of the boundaries instead of trimming.
#' `method = "gpv"` gives the conventional benchmark
#' \eqn{\hat v = b + \hat G(b)/\{(n-1)\hat g(b)\}} with a Gaussian kernel and
#' Sheather--Jones bandwidth for \eqn{\hat g}; `"gpv_bc"` renormalises a
#' triweight kernel by its mass inside the bid support (a simple boundary
#' correction; the Hickman--Hubbard correction is provided by
#' `fpa_density()`).
#'
#' @section Monotonicity:
#' The equilibrium bid function is increasing, so \eqn{\hat v} should be
#' increasing in \eqn{b}. Departures occur mostly near the boundaries.
#' `monotone = "spline"` follows Doosti, Dewan and Talebian (2025): a natural
#' cubic spline of \eqn{\hat v} on \eqn{b} is fitted to the central
#' \eqn{1 - 2\,\mathrm{trim}} share of the bids and its fitted values replace
#' \eqn{\hat v} either everywhere (`replace = "all"`, the published
#' implementation) or only in the trimmed tails (`replace = "boundary"`).
#' `monotone = "isotonic"` applies the pool-adjacent-violators rearrangement
#' ([stats::isoreg()]). `"none"` returns the raw pseudo-values.
#'
#' @param bids numeric vector of pooled bids, or an `L` by `n` matrix (one row
#'   per auction). A matrix is required for `family != "independence"` with
#'   `theta = NULL`, so that the copula parameter can be estimated from
#'   within-auction pairs.
#' @param n number of bidders per auction (taken from `ncol(bids)` for a matrix).
#' @param method quantile-density estimator passed to [qdf()], or `"gpv"`,
#'   `"gpv_bc"`, `"zincenko"` for the kernel benchmarks.
#' @param bandwidth passed to [qdf()]; a selector name or a numeric bandwidth.
#' @param eta CRRA risk-aversion parameter in \eqn{[0, 1)}; `0` is risk
#'   neutrality.
#' @param family copula family of the private values; see [copula_psi()].
#' @param theta copula parameter; `NULL` estimates it by [copula_fit()].
#' @param u_scheme scores \eqn{\hat u} attached to the sorted bids:
#'   `"i/(n+1)"` (Hazen ranks, default) or `"i/n"`.
#' @param monotone `"none"`, `"spline"` or `"isotonic"`.
#' @param trim tail share excluded from the spline fit (default 0.1, i.e. the
#'   central 80 percent is used).
#' @param replace for `monotone = "spline"`: `"all"` or `"boundary"`.
#' @param df degrees of freedom of the natural spline (default 6).
#' @param ... further arguments to [qdf()] (e.g. `kernel`, `H`, `j0`, `m_range`).
#' @return An object of class `"fpa_values"`: a list with
#' \describe{
#'   \item{bids, u}{sorted pooled bids and their scores}
#'   \item{q, psi, markup}{the quantile density, the copula multiplier and
#'     the estimated markup \eqn{\lambda^{-1}(\psi q)} at each bid}
#'   \item{values_raw, values}{pseudo-values before and after the
#'     monotonicity step (identical when `monotone = "none"`)}
#'   \item{values_mat}{when `bids` is a matrix, `values` rearranged to the
#'     original `L` by `n` layout}
#'   \item{qdf}{the fitted [qdf()] object (or `NULL` for the kernel benchmarks)}
#'   \item{copula}{the copula specification / fit}
#' }
#' @references
#' Guerre, E., Perrigne, I. and Vuong, Q. (2000). Optimal nonparametric
#' estimation of first-price auctions. *Econometrica*, 68, 525--574.
#'
#' Guerre, E., Perrigne, I. and Vuong, Q. (2009). Nonparametric identification
#' of risk aversion in first-price auctions under exclusion restrictions.
#' *Econometrica*, 77, 1193--1227.
#'
#' Doosti, H., Dewan, I. and Talebian, M. (2025). *Economics Letters*, 257, 112670.
#'
#' Doosti, H. (2026a). *Economics Letters*, 261, 112837.
#' @examples
#' set.seed(1)
#' sim <- fpa_simulate(20, 5, Q = function(u) qgamma(u, 2, 10))
#' fit <- fpa_values(sim$bids, method = "bernstein", bandwidth = "bcv")
#' fit
#' plot(fit, true = sim$values)
#' # risk-averse bidders with CRRA eta = 0.25
#' fit_ra <- fpa_values(sim$bids, method = "bernstein", bandwidth = 0.05, eta = 0.25)
#' # affiliated values, Clayton copula estimated from the data
#' \donttest{
#' sim2 <- fpa_simulate(50, 3, Q = function(u) 1 + 2 * u^2, family = "clayton", theta = 2)
#' fit2 <- fpa_values(sim2$bids, method = "indirect_poisson", family = "clayton")
#' fit2$copula$theta
#' }
#' @export
fpa_values <- function(bids, n = NULL, method = "indirect_poisson", bandwidth = "bcv",
                       eta = 0, family = "independence", theta = NULL,
                       u_scheme = c("i/(n+1)", "i/n"),
                       monotone = c("none", "spline", "isotonic"),
                       trim = 0.1, replace = c("all", "boundary"), df = 6, ...) {
  u_scheme <- match.arg(u_scheme); monotone <- match.arg(monotone); replace <- match.arg(replace)
  family <- match.arg(family, .copula_families)
  if (eta < 0 || eta >= 1) stop("`eta` must be in [0, 1).", call. = FALSE)
  is_mat <- is.matrix(bids)
  if (is_mat) { n <- ncol(bids); L <- nrow(bids) }
  if (is.null(n) || n < 2) stop("`n` (bidders per auction, >= 2) is required.", call. = FALSE)
  b_all <- as.numeric(bids)
  if (any(!is.finite(b_all))) stop("`bids` must be finite.", call. = FALSE)
  o <- order(b_all); b <- b_all[o]; N <- length(b)
  u <- if (u_scheme == "i/(n+1)") seq_len(N) / (N + 1) else seq_len(N) / N

  # copula multiplier
  cop <- list(family = family, theta = theta)
  if (family != "independence" && is.null(theta)) {
    if (!is_mat) stop("Estimating `theta` requires `bids` as an L x n matrix.", call. = FALSE)
    cop <- copula_fit(bids, family)
  }
  psi <- copula_psi(u, family, cop$theta, n)

  # quantile density of bids at u
  qfit <- NULL
  if (method == "zincenko") {
    if (family != "independence") stop("`method = \"zincenko\"` is defined for independent private values.", call. = FALSE)
    v_raw <- values_zincenko(b, n, eta)
    markup <- v_raw - b; q <- rep(NA_real_, N); psi_q <- markup / (1 - eta)
  } else if (method %in% c("gpv", "gpv_bc")) {
    G <- u
    if (method == "gpv") {
      g <- kde_at(b, b, stats::bw.SJ(b), "gaussian")
    } else {
      g <- kde_bc_at(b, b)
    }
    q <- 1 / g
    psi_q <- G / ((n - 1) * g)       # exact GPV form; equals psi * q under independence
    if (family != "independence") psi_q <- psi * q
  } else {
    qfit <- qdf(b, method = method, bandwidth = bandwidth, u = u, ...)
    q <- qfit$q
    psi_q <- psi * q
  }
  if (method != "zincenko") {
    markup <- (1 - eta) * psi_q
    v_raw <- b + markup
  }

  v <- switch(monotone,
              none = v_raw,
              spline = monotone_spline(b, v_raw, trim, replace, df),
              isotonic = stats::isoreg(b, v_raw)$yf)

  out <- list(bids = b, u = u, n = n, q = q, psi = psi, markup = markup, eta = eta,
              values_raw = v_raw, values = v, method = method, monotone = monotone,
              qdf = qfit, copula = cop, order = o, call = match.call())
  if (is_mat) {
    vm <- numeric(N); vm[o] <- v
    out$values_mat <- matrix(vm, nrow = L, ncol = n)
  }
  class(out) <- "fpa_values"
  out
}

# Gaussian/other KDE at points
kde_at <- function(x, data, h, kernel = "gaussian") {
  K <- kernel_fun(kernel)
  rowMeans(K$d(outer(x, data, "-") / h)) / h
}

# triweight KDE with Silverman bandwidth, renormalised by the kernel mass inside the support
kde_bc_at <- function(x, data) {
  n <- length(data)
  h <- 1.06 * min(stats::sd(data), stats::IQR(data) / 1.34) * n^(-1 / 5)
  tw <- function(z) (35 / 32) * (1 - z^2)^3 * (abs(z) < 1)
  raw <- rowMeans(tw(outer(x, data, "-") / h)) / h
  lo <- min(data); hi <- max(data)
  Ptw <- function(z) { z <- pmin(pmax(z, -1), 1); 0.5 + (35 / 32) * (z - z^3 + 3 * z^5 / 5 - z^7 / 7) }
  nrm <- Ptw((x - lo) / h) - Ptw((x - hi) / h)  # mass of K((x - .)/h)/h on [lo, hi]
  pmax(raw / pmax(nrm, 1e-3), 1e-8)
}

# natural-spline monotonicity fix of Doosti, Dewan & Talebian (2025)
monotone_spline <- function(b, v, trim = 0.1, replace = "all", df = 6) {
  N <- length(b)
  i0 <- floor(trim * N) + 1L; i1 <- ceiling((1 - trim) * N)
  idx <- i0:i1
  fit <- stats::lm(v[idx] ~ splines::ns(b[idx], df = df))
  basis <- splines::ns(b[idx], df = df)
  fitted_all <- as.numeric(cbind(1, stats::predict(basis, b)) %*% stats::coef(fit))
  if (replace == "all") return(fitted_all)
  out <- v; keep <- setdiff(seq_len(N), idx); out[keep] <- fitted_all[keep]
  out
}

#' @export
print.fpa_values <- function(x, digits = 4, ...) {
  cat("Private values recovered from first-price auction bids\n")
  cat("  bids      :", length(x$bids), "pooled,", x$n, "bidders per auction\n")
  cat("  method    :", x$method, if (!is.null(x$qdf)) sprintf("(h = %s)", format(x$qdf$h, digits = digits)), "\n")
  cat("  utility   :", if (x$eta == 0) "risk neutral" else sprintf("CRRA, eta = %s", x$eta), "\n")
  cat("  copula    :", x$copula$family,
      if (x$copula$family != "independence") sprintf("(theta = %s)", format(x$copula$theta, digits = digits)), "\n")
  cat("  monotone  :", x$monotone, "\n")
  mk <- x$values - x$bids
  cat("  markup    : median", format(stats::median(mk), digits = digits),
      "(", format(100 * stats::median(mk / x$bids), digits = 3), "% of bid )\n")
  invisible(x)
}

#' Plot recovered private values against bids
#' @param x an `"fpa_values"` object.
#' @param true optional true values (vector or matrix in the layout of `bids`)
#'   for simulated data.
#' @param ... passed to [graphics::plot()].
#' @export
plot.fpa_values <- function(x, true = NULL, ...) {
  args <- list(x = x$bids, y = x$values_raw, pch = 20, col = "grey50",
               xlab = "bid", ylab = "private value", main = "Inverse bid function")
  args <- utils::modifyList(args, list(...))
  do.call(graphics::plot, args)
  if (x$monotone != "none") graphics::lines(x$bids, x$values, col = 4, lwd = 2)
  if (!is.null(true)) {
    tv <- as.numeric(true)[x$order]
    graphics::lines(x$bids, tv, col = 2, lwd = 2, lty = 2)
    graphics::legend("topleft", c("pseudo-values", if (x$monotone != "none") "monotone fit", "true"),
                     col = c("grey50", if (x$monotone != "none") 4, 2), pch = c(20, NA, NA)[c(TRUE, x$monotone != "none", TRUE)],
                     lty = c(NA, 1, 2)[c(TRUE, x$monotone != "none", TRUE)], bty = "n")
  }
  invisible(x)
}
