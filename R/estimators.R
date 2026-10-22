# Quantile density estimators -------------------------------------------------
#
# Revised and vectorised from R code by Y. P. Chaubey, I. Dewan and J. Li
# (Chaubey, Dewan and Li 2021, 2024), (C) the original authors; see
# inst/COPYRIGHTS.
#
# All functions take an evaluation grid `u` in [0, 1], a *sorted* sample `x`
# and a smoothing parameter `h`. They return the estimated quantile density
# q(u) = Q'(u) = 1 / g(Q(u)) at each u. They are the workhorses behind
# [qdf()] and are exported so that users can call them directly with a fixed
# bandwidth.
#
# Notation follows Doosti, Dewan and Talebian (2025, Section 2) and
# Chaubey, Dewan and Li (2024). Where an estimator is a "series" estimator
# (Poisson, Bernstein, indirect Poisson) the smoothing parameter is the
# integer `m = floor(1/h)`; passing `h` keeps a single interface.

#' Direct kernel estimators of the quantile density
#'
#' @description
#' Kernel-type estimators that smooth the empirical quantile function directly.
#'
#' * `qdf_kernel()` is the Jones (1992) estimator
#'   \deqn{\hat q(u) = \sum_{i=1}^n X_{(i)} \{ k_h((i-1)/n - u) - k_h(i/n - u) \}}
#'   with a Gaussian kernel \eqn{k_h(a) = k(a/h)/h}.
#' * `qdf_kernel_corrected()` is the boundary-corrected version (estimator 2 of
#'   Doosti, Dewan and Talebian 2025), which renormalises the kernel weights by
#'   the mass falling inside \eqn{[0,1]} and adds a correction proportional to
#'   the kernel quantile estimator.
#'
#' @param u numeric vector of evaluation points in \eqn{[0, 1]}.
#' @param x numeric sample (sorted internally).
#' @param h bandwidth on the probability scale (\eqn{0 < h < 1}).
#' @param kernel `"gaussian"` (default), `"triangular"` or `"epanechnikov"`, or
#'   a user-defined kernel from [make_kernel()] (or a kernel density function).
#'   Soni, Dewan and Jain (2012) used the triangular and Epanechnikov kernels.
#' @return numeric vector of quantile density estimates at `u`.
#' @references
#' Jones, M. C. (1992). Estimating densities, quantiles, quantile densities and
#' density quantiles. *Annals of the Institute of Statistical Mathematics*, 44,
#' 721--727.
#'
#' Chaubey, Y. P., Dewan, I. and Li, J. (2024). On some non parametric
#' estimators of the quantile density function for a stationary associated
#' process. *Communications in Statistics -- Theory and Methods*, 53, 5553--5573.
#' @seealso [qdf()] for automatic bandwidth selection.
#' @examples
#' x <- rgamma(100, 5, 1)
#' u <- seq(0.05, 0.95, by = 0.05)
#' qdf_kernel(u, x, h = 0.05)
#' qdf_kernel_corrected(u, x, h = 0.05)
#' @export
qdf_kernel <- function(u, x, h, kernel = "gaussian") {
  x <- prep_sample(x); u <- check_u(u); n <- length(x); K <- kernel_fun(kernel)
  lo <- (seq_len(n) - 1) / n
  hi <- seq_len(n) / n
  # rows: u, cols: i
  W <- (K$d(outer(u, lo, "-") / h) - K$d(outer(u, hi, "-") / h)) / h
  as.numeric(W %*% x)
}

#' @rdname qdf_kernel
#' @export
qdf_kernel_corrected <- function(u, x, h, kernel = "gaussian") {
  x <- prep_sample(x); u <- check_u(u); n <- length(x); K <- kernel_fun(kernel)
  s <- (0:n) / n
  # mass of the kernel inside [0, 1] for each interval [s_j, s_{j+1}]
  Wi <- K$p(outer(s[-1], u, "-") / h) - K$p(outer(s[-(n + 1)], u, "-") / h) # n x |u|
  sumW <- colSums(Wi)
  Qk   <- as.numeric(x %*% Wi) / sumW                   # kernel quantile estimator
  kh   <- function(a) K$d(a / h) / h
  L    <- seq_len(n - 1) / n
  term1 <- as.numeric(crossprod(kh(outer(L, u, "-")), diff(x))) # sum (X_(i+1)-X_(i)) k_h(i/n - u)
  k1 <- kh(1 - u); k0 <- kh(-u)
  (term1 - x[n] * k1 + x[1] * k0 + (k1 - k0) * Qk) / sumW
}

#' Series estimators of the quantile density
#'
#' @description
#' Estimators that smooth the empirical quantile function \eqn{Q_n} with
#' probability weights on a lattice \eqn{k/m}, \eqn{m = \lfloor 1/h \rfloor}.
#'
#' * `qdf_poisson()` uses Poisson weights (Chaubey, Dewan and Li 2021; estimator
#'   3 of Doosti, Dewan and Talebian 2025):
#'   \deqn{\hat q(u) = m \sum_{k=0}^{m} Q_n(k/m)\, w_k'(m u)}
#'   where \eqn{w_k(\mu) = p_k(\mu)/P_m(\mu)} are truncated Poisson weights.
#' * `qdf_bernstein()` uses Bernstein polynomials (estimator 4):
#'   \deqn{\hat q(u) = m \sum_{k=0}^{m-1} \{Q_n((k+1)/m) - Q_n(k/m)\} b_k(u; m-1).}
#'
#' Both require non-negative data because \eqn{Q_n(0)} is taken to be 0;
#' `qdf_bernstein()` can instead use \eqn{Q_n(0) = X_{(1)}}
#' (`zero_at_origin = FALSE`), which is preferable when the data do not start
#' near zero (for example bids in dollars).
#'
#' `qdf_bernstein_loo()` returns the exact leave-one-out Bernstein estimates
#' \eqn{\hat q_{-i}(u_i)}, \eqn{i = 1, \dots, n}, used by the cross-validation
#' criteria: removing \eqn{X_{(i)}} shifts every higher order statistic down
#' by one position, so the leave-one-out lattice values are obtained by an
#' index shift in \eqn{O(nK)} operations, without refitting.
#'
#' @inheritParams qdf_kernel
#' @param zero_at_origin logical; `TRUE` (default) takes \eqn{Q_n(0) = 0},
#'   `FALSE` takes \eqn{Q_n(0) = X_{(1)}}.
#' @return numeric vector of quantile density estimates at `u`; for
#'   `qdf_bernstein_loo()`, the vector
#'   \eqn{(\hat q_{-1}(u_1), \dots, \hat q_{-n}(u_n))}, where `u` has the
#'   length of `x` and is paired with the sorted sample.
#' @references
#' Chaubey, Y. P., Dewan, I. and Li, J. (2021). On some smooth estimators of
#' the quantile function for a stationary associated process. *Sankhya B*,
#' 83 (Suppl 1), 114--139. \doi{10.1007/s13571-020-00242-x}
#' @examples
#' x <- rgamma(100, 5, 1)
#' qdf_poisson(c(0.25, 0.5, 0.75), x, h = 0.05)
#' qdf_bernstein(c(0.25, 0.5, 0.75), x, h = 0.05)
#' # exact leave-one-out estimates at u_i = i/n agree with refitting
#' xs <- sort(x); u <- seq_along(xs) / length(xs)
#' loo <- qdf_bernstein_loo(u, xs, h = 0.05)
#' all.equal(loo[10], qdf_bernstein(u[10], xs[-10], h = 0.05))
#' @export
qdf_poisson <- function(u, x, h) {
  x <- prep_sample(x, positive = TRUE); u <- check_u(u)
  m  <- 1 / h
  l  <- floor(m)
  eq <- eq_lattice((0:l) / m, x)
  vapply(u, function(ui) {
    W   <- pois_w(ui * m, l + 1L)                       # p_0..p_l at mu = m u
    wpp <- c(-W[1L], W[1:l] - W[2:(l + 1L)])            # p_k'(mu) = p_{k-1} - p_k
    Pn  <- sum(W)
    wp  <- (wpp + W * W[l + 1L] / Pn) / Pn              # w_k'(mu)
    m * sum(wp * eq)
  }, numeric(1))
}

#' @rdname qdf_poisson
#' @export
qdf_bernstein <- function(u, x, h, zero_at_origin = TRUE) {
  x <- prep_sample(x, positive = TRUE); u <- check_u(u)
  m   <- 1 / h
  K   <- floor(m)
  eq  <- eq_lattice((0:K) / m, x, zero_at_origin = zero_at_origin)
  deq <- diff(eq)
  B   <- outer(u, 0:(K - 1L), function(uu, k) stats::dbinom(k, K - 1L, uu))
  m * as.numeric(B %*% deq)
}

#' @rdname qdf_poisson
#' @export
# Exact leave-one-out Bernstein QDF estimates at the grid points u_i, i = 1..n:
# removing x_(i) shifts every order statistic above it by one position, so the
# lattice values of the leave-one-out empirical quantile function are
# x_(j + 1{j >= i}) with j the lattice index for a sample of size n - 1.
# Cost O(n K) per bandwidth, the same as the estimator itself.
qdf_bernstein_loo <- function(u, x, h, zero_at_origin = TRUE) {
  xs <- prep_sample(x, positive = TRUE); n <- length(xs)
  if (!is.numeric(u) || length(u) != n)
    stop("`u` must have one value per observation (u[i] is paired with the i-th ",
         "order statistic); got length ", length(u), " for ", n, " observations.", call. = FALSE)
  if (n < 2L) stop("`x` must contain at least 2 observations.", call. = FALSE)
  if (is.unsorted(u, strictly = TRUE) || any(u < 0 | u > 1))
    stop("`u` must be strictly increasing in [0, 1], e.g. (1:n)/(n + 1).", call. = FALSE)
  if (length(h) != 1L || !is.finite(h) || h <= 0 || floor(1 / h) < 2)
    stop("`h` must be a single positive number with floor(1/h) >= 2.", call. = FALSE)
  m <- 1 / h; K <- floor(m); t <- (0:K) / m
  j  <- pmin(floor((n - 1) * t) + 1L, n - 1L)
  B  <- outer(u, 0:(K - 1L), function(uu, k) stats::dbinom(k, K - 1L, uu))
  out <- numeric(n)
  for (i in seq_len(n)) {
    eq <- xs[j + (j >= i)]
    if (zero_at_origin) eq[t <= 0] <- 0
    eq[t >= 1] <- if (i == n) xs[n - 1L] else xs[n]
    out[i] <- m * sum(B[i, ] * diff(eq))
  }
  out
}

#' Indirect estimators of the quantile density
#'
#' @description
#' Estimators that exploit \eqn{q(u) = 1 / g(Q(u))} by first estimating the
#' density \eqn{g}.
#'
#' * `qdf_jones()` (Jones 1992; estimator 5) is
#'   \eqn{\hat q(u) = 1/\hat g_h(Q_n(u))} with a Gaussian kernel density
#'   estimator \eqn{\hat g_h}.
#' * `qdf_soni()` (Soni, Dewan and Jain 2012; estimator 6) double-smooths
#'   Jones' estimator on the probability scale. With `form = "integrated"`
#'   (the form printed in Doosti, Dewan and Talebian 2025),
#'   \deqn{\hat q(u) = \sum_{i=1}^n \frac{K((i/n - u)/H) - K(((i-1)/n - u)/H)}{\hat g_h(X_{(i)})},}
#'   where \eqn{K} is the kernel's distribution function; with
#'   `form = "riemann"` (the original implementation of Soni, Dewan and Jain)
#'   \deqn{\hat q(u) = \frac{1}{nH}\sum_{i=1}^n \frac{k((i/n - u)/H)}{\hat g_h(X_{(i)})}.}
#'   Here \eqn{h} is the density bandwidth on the data scale and \eqn{H} the
#'   smoothing bandwidth on the probability scale. For right-censored data
#'   supply `status`; the empirical scores \eqn{i/n} are then replaced by the
#'   Kaplan--Meier estimate \eqn{\hat F(X_{(i)})}, as in Soni et al.
#' * `qdf_indirect_poisson()` (Chaubey, Dewan and Li 2024; estimator 7)
#'   smooths the empirical *distribution* function with Poisson weights,
#'   inverts the result numerically to obtain a smooth quantile estimator, and
#'   takes the reciprocal of the implied density at that quantile. It requires
#'   positive data.
#'
#' @inheritParams qdf_kernel
#' @param h for `qdf_jones()` and `qdf_soni()` the kernel density bandwidth on
#'   the data scale; for `qdf_indirect_poisson()` the smoothing parameter
#'   \eqn{1/m} of the Poisson weights on the data scale.
#' @param H for `qdf_soni()`, the bandwidth of the outer Gaussian smoother on
#'   the probability scale.
#' @param tol root-finding tolerance passed to [stats::uniroot()] when
#'   inverting the smoothed distribution function.
#' @param kernel `"gaussian"` (default), `"triangular"`, `"epanechnikov"`, or a
#'   user-defined kernel from [make_kernel()].
#' @param form for `qdf_soni()`: `"integrated"` (default) or `"riemann"`.
#' @param status for `qdf_soni()`: optional 0/1 vector (1 = observed,
#'   0 = right-censored) aligned with `x`.
#' @return numeric vector of quantile density estimates at `u`.
#' @references
#' Soni, P., Dewan, I. and Jain, K. (2012). Nonparametric estimation of
#' quantile density function. *Computational Statistics & Data Analysis*, 56,
#' 3876--3886.
#' @examples
#' x <- rgamma(100, 5, 1)
#' u <- c(0.25, 0.5, 0.75)
#' qdf_jones(u, x, h = 0.5)
#' qdf_soni(u, x, h = 0.5, H = 0.05)
#' qdf_indirect_poisson(u, x, h = 0.05)
#' @export
qdf_jones <- function(u, x, h, kernel = "gaussian") {
  x <- prep_sample(x); u <- check_u(u); K <- kernel_fun(kernel)
  q0 <- eq_lattice(u, x, zero_at_origin = FALSE)
  g  <- rowMeans(K$d(outer(q0, x, "-") / h)) / h
  1 / g
}

#' @rdname qdf_jones
#' @export
qdf_soni <- function(u, x, h, H, kernel = "gaussian",
                     form = c("integrated", "riemann"), status = NULL) {
  form <- match.arg(form); u <- check_u(u); K <- kernel_fun(kernel)
  if (!is.null(status)) {
    if (length(status) != length(x)) stop("`status` must have the same length as `x`.", call. = FALSE)
    o <- order(x); x <- x[o]; status <- as.numeric(status[o])
    Fi <- km_cdf(x, status)
  } else {
    x <- prep_sample(x); Fi <- seq_len(length(x)) / length(x)
  }
  n  <- length(x)
  gx <- rowMeans(K$d(outer(x, x, "-") / h)) / h          # KDE at the order statistics
  if (form == "riemann") {
    Wi <- K$d(outer(Fi, u, "-") / H) / (n * H)             # n x |u|
  } else {
    s  <- c(0, Fi)
    Wi <- K$p(outer(s[-1], u, "-") / H) - K$p(outer(s[-(n + 1)], u, "-") / H)
  }
  as.numeric((1 / gx) %*% Wi)
}

#' @rdname qdf_jones
#' @param invert how the smoothed distribution function is inverted:
#'   `"newton"` tabulates it on a grid and refines by Newton steps (efficient
#'   for many evaluation points); `"uniroot"` reproduces the original
#'   implementation with [stats::uniroot()] (efficient for a few points;
#'   accuracy governed by `tol`); `"auto"` (default) picks by `length(u)`.
#' @export
qdf_indirect_poisson <- function(u, x, h, tol = .Machine$double.eps^0.25,
                                 invert = c("auto", "newton", "uniroot")) {
  invert <- match.arg(invert)
  x <- prep_sample(x, positive = TRUE); u <- check_u(u)
  sm <- ip_smoother(x, h)
  r  <- ip_invert(u, sm, invert, tol)
  1 / sm$gs(r)
}

# Poisson-smoothed distribution function and its density (Chaubey, Dewan & Li 2024)
# with the Poisson weights truncated to a window around their mean
ip_smoother <- function(x, h) {
  m <- 1 / h; M <- ceiling(max(x) * m)
  Sk <- esurv_at((0:M) / m, x)
  w  <- ecdf_at((1:M) / m, x) - ecdf_at((0:(M - 1L)) / m, x)
  win <- function(mu, kmax) {
    s <- 12 * sqrt(mu + 1) + 25
    c(min(max(0, floor(mu - s)), kmax), min(kmax, ceiling(mu + s)))
  }
  Gs <- function(t) vapply(t, function(tt) {
    mu <- tt * m; k <- win(mu, M); kk <- k[1]:k[2]
    1 - sum(stats::dpois(kk, mu) * Sk[kk + 1L])
  }, numeric(1))
  gs <- function(t) vapply(t, function(tt) {
    mu <- tt * m; k <- win(mu, M - 1L); kk <- k[1]:k[2]
    m * sum(stats::dpois(kk, mu) * w[kk + 1L])
  }, numeric(1))
  list(Gs = Gs, gs = gs, xmax = max(x), m = m, M = M)
}

ip_invert <- function(u, sm, invert = "auto", tol = .Machine$double.eps^0.25) {
  out <- numeric(length(u)); pos <- u > 0
  if (!any(pos)) return(out)
  if (invert == "auto") invert <- if (sum(pos) <= 24L) "uniroot" else "newton"
  if (invert == "uniroot") {
    out[pos] <- vapply(u[pos], function(ui)
      stats::uniroot(function(t) sm$Gs(t) - ui, c(0, sm$xmax), extendInt = "upX", tol = tol)$root,
      numeric(1))
    return(out)
  }
  ng <- 256L
  tg <- seq(0, sm$xmax * 1.5, length.out = ng)
  Gg <- sm$Gs(tg)
  Gg <- cummax(Gg) + seq_len(ng) * 1e-12            # strictly increasing for interpolation
  t0 <- stats::approx(Gg, tg, xout = u[pos], rule = 2)$y
  for (it in 1:3) {
    step <- (sm$Gs(t0) - u[pos]) / pmax(sm$gs(t0), 1e-12)
    step <- pmin(pmax(step, -sm$xmax / 4), sm$xmax / 4)
    t0 <- pmax(t0 - step, 0)
  }
  out[pos] <- t0
  out
}

#' Smoothed quantile function implied by the indirect Poisson estimator
#'
#' Returns the Poisson-smoothed quantile function \eqn{\tilde Q(u)} used by
#' [qdf_indirect_poisson()]; useful for plotting or for constructing
#' pseudo-values.
#' @inheritParams qdf_indirect_poisson
#' @return numeric vector of smoothed quantiles at `u`.
#' @keywords internal
#' @export
quantile_indirect_poisson <- function(u, x, h, tol = .Machine$double.eps^0.25,
                                      invert = c("auto", "newton", "uniroot")) {
  invert <- match.arg(invert)
  x <- prep_sample(x, positive = TRUE); u <- check_u(u)
  ip_invert(u, ip_smoother(x, h), invert, tol)
}
