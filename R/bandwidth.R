# Bandwidth selection ----------------------------------------------------------

# internal: evaluate an estimator on the CV grid with and without observation i
cv_fits <- function(x, h, fun, u, loo = TRUE) {
  n  <- length(x)
  q  <- fun(u, x, h)
  q_loo <- if (loo) {
    vapply(seq_len(n), function(i) fun(u[i], x[-i], h), numeric(1))
  } else {
    q * (n - 1) / n          # linear-smoother approximation (Doosti 2026b, Remark 2)
  }
  list(q = q, q_loo = q_loo)
}

#' Cross-validation criteria for quantile density estimators
#'
#' @description
#' Three data-driven criteria for the smoothing parameter of a quantile density
#' estimator, written as functions of \eqn{m = 1/h} so that they can be passed
#' to [stats::optimize()]. Each is evaluated on the grid \eqn{u_i = i/n} with
#' spacings \eqn{\Delta_i = X_{(i)} - X_{(i-1)}}, \eqn{X_{(0)} = 0}, so that
#' \eqn{\sum_i \Delta_i\, \phi(\hat q(u_i))} approximates \eqn{\int \phi\{1/\hat g(x)\}\,dx}.
#'
#' * `bcv_criterion()` is the biased (least-squares) cross-validation used by
#'   Chaubey, Dewan and Li (2024) and Doosti, Dewan and Talebian (2025):
#'   \deqn{\mathrm{BCV}(h) = \sum_i \Delta_i\, \hat q(u_i;h)^{-2} - \frac{2}{n}\sum_i \hat q_{-i}(u_i;h)^{-1},}
#'   i.e. the classical \eqn{\int \hat g^2 - 2\,\bar{\hat g}_{-i}} written in
#'   terms of \eqn{\hat q = 1/\hat g}. (Equation (11) of Doosti, Dewan and
#'   Talebian (2025) prints \eqn{\hat q^2} in the first term; the simulation
#'   code, reproduced here, uses \eqn{\hat q^{-2}}.)
#' * `rlcv_criterion()` is the robust likelihood cross-validation of Wu (2019),
#'   \eqn{-\{ \frac1n \sum_i \ln^*(\hat q_{-i}(u_i)^{-1}) - b^*(h) \}}, with
#'   \eqn{\ln^*} linearised below the threshold \eqn{a_n}. It is written so that
#'   it is *minimised*.
#' * `wbcv_criterion()` is the weighted criterion of Doosti (2026b), which
#'   targets \eqn{u\,q(u)} (the quantity entering the first-price auction
#'   inversion) rather than \eqn{q(u)}:
#'   \deqn{\mathrm{WBCV}(h) = \sum_i \Delta_i\{u_i \hat q(u_i)\}^2 - \frac{2}{n}\sum_i u_i / \hat q_{-i}(u_i).}
#'
#' @param m inverse bandwidth, \eqn{m = 1/h}.
#' @param x numeric sample (sorted internally).
#' @param fun an estimator with signature `fun(u, x, h)`, e.g. [qdf_kernel()].
#'   Extra arguments (such as `H` for [qdf_soni()]) should be bound with a
#'   closure before calling.
#' @param loo logical; if `TRUE` (default) leave-one-out fits are computed
#'   exactly by refitting without observation \eqn{i}. If `FALSE` the
#'   approximation \eqn{\hat q_{-i} \approx \hat q\,(n-1)/n} is used, which is
#'   exact up to \eqn{O(n^{-2})} for linear smoothers and much faster.
#' @param u_scheme grid on which the criterion is evaluated: `"i/n"` (the
#'   published convention) or `"i/(n+1)"`.
#' @return a single numeric value; smaller is better.
#' @references
#' Wu, X. (2019). Robust likelihood cross validation for kernel density
#' estimation. *Journal of Business & Economic Statistics*, 37, 761--770.
#' @seealso [select_bandwidth()], [qdf()]
#' @examples
#' x <- rgamma(60, 5, 1)
#' bcv_criterion(20, x, qdf_kernel)
#' rlcv_criterion(20, x, qdf_kernel)
#' optimize(bcv_criterion, c(10, 50), x = x, fun = qdf_bernstein)$minimum
#' @export
bcv_criterion <- function(m, x, fun, loo = TRUE, u_scheme = c("i/n", "i/(n+1)")) {
  x <- prep_sample(x); n <- length(x); h <- 1 / m
  u <- default_u(n, u_scheme)
  d <- diff(c(0, x))
  f <- cv_fits(x, h, fun, u, loo)
  sum(d / f$q^2) - 2 * mean(1 / f$q_loo)
}

#' @rdname bcv_criterion
#' @export
rlcv_criterion <- function(m, x, fun, loo = TRUE, u_scheme = c("i/n", "i/(n+1)")) {
  x <- prep_sample(x); n <- length(x); h <- 1 / m
  u <- default_u(n, u_scheme)
  d <- diff(c(0, x))
  f <- cv_fits(x, h, fun, u, loo)
  a <- a_n(x)
  g     <- 1 / f$q                    # implied density along the quantile grid
  g_loo <- 1 / f$q_loo
  ll    <- mean(log_star(g_loo, a))
  bstar <- sum(d * g * (g >= a)) + sum(d * g^2 * (g < a)) / (2 * a)
  -(ll - bstar)
}

#' @rdname bcv_criterion
#' @export
wbcv_criterion <- function(m, x, fun, loo = FALSE, u_scheme = c("i/(n+1)", "i/n")) {
  x <- prep_sample(x); n <- length(x); h <- 1 / m
  u <- default_u(n, u_scheme)
  d <- diff(c(0, x))
  f <- cv_fits(x, h, fun, u, loo)
  q <- pmax(f$q, 1e-10); q_loo <- pmax(f$q_loo, 1e-10)
  sum(d * (u * q)^2) - 2 * mean(u / q_loo)
}

#' Select the smoothing parameter of a quantile density estimator
#'
#' Minimises one of the criteria in [bcv_criterion()] over \eqn{m = 1/h} in
#' `m_range` with [stats::optimize()], optionally after a coarse grid search.
#'
#' @inheritParams bcv_criterion
#' @param selector `"bcv"`, `"rlcv"` or `"wbcv"`.
#' @param m_range numeric length-2 range for \eqn{m = 1/h}.
#' @param grid integer; if positive, the criterion is first evaluated on this
#'   many equally spaced points in `m_range` and [stats::optimize()] is then
#'   run in the bracket around the grid minimum. This guards against local
#'   minima of the (often non-convex) criteria. `0` runs `optimize()` directly
#'   over `m_range`, which reproduces the published simulations.
#' @return a list with `h`, `m`, `selector`, `criterion` (the minimised value)
#'   and `m_range`.
#' @examples
#' x <- rgamma(60, 5, 1)
#' select_bandwidth(x, qdf_poisson, "bcv", m_range = c(10, 50))
#' @export
select_bandwidth <- function(x, fun, selector = c("bcv", "rlcv", "wbcv"),
                             m_range = c(10, 50), loo = NULL, grid = 0L,
                             u_scheme = NULL) {
  selector <- match.arg(selector)
  crit <- switch(selector, bcv = bcv_criterion, rlcv = rlcv_criterion,
                 wbcv = wbcv_criterion)
  if (is.null(loo)) loo <- selector != "wbcv"
  if (is.null(u_scheme)) u_scheme <- if (selector == "wbcv") "i/(n+1)" else "i/n"
  obj <- function(m) {
    v <- tryCatch(crit(m, x, fun, loo = loo, u_scheme = u_scheme),
                  error = function(e) NA_real_)
    if (!is.finite(v)) Inf else v
  }
  lo <- m_range[1]; hi <- m_range[2]
  if (grid > 0L) {
    mg <- seq(lo, hi, length.out = grid)
    cv <- vapply(mg, obj, numeric(1))
    if (!any(is.finite(cv))) stop("Criterion is non-finite over the whole range.", call. = FALSE)
    k  <- which.min(cv)
    lo <- mg[max(1L, k - 1L)]; hi <- mg[min(grid, k + 1L)]
    if (hi <= lo) { lo <- m_range[1]; hi <- m_range[2] }
  }
  o <- stats::optimize(obj, c(lo, hi))
  if (!is.finite(o$objective))
    stop("Bandwidth selection failed: criterion is non-finite. ",
         "Try a different `m_range` or estimator.", call. = FALSE)
  list(h = 1 / o$minimum, m = o$minimum, selector = selector,
       criterion = o$objective, m_range = m_range)
}
