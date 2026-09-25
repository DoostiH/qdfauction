# Internal helpers -----------------------------------------------------------

#' Empirical quantile function at lattice points
#'
#' Left-continuous empirical quantile used by the series estimators:
#' `Q_n(u) = X_(floor(n u) + 1)` for `0 < u < 1`, `Q_n(0) = 0`, `Q_n(1) = X_(n)`.
#' This is the convention of Chaubey, Dewan and Li (2021, 2024) and of the
#' original simulation code; note that `Q_n(0) = 0` (not `X_(1)`) because the
#' series estimators are defined for non-negative data. Set
#' `zero_at_origin = FALSE` to get the usual `Q_n(0) = X_(1)`.
#'
#' @param u numeric vector in `[0, 1]`.
#' @param x sorted numeric sample.
#' @return numeric vector of the same length as `u`.
#' @keywords internal
#' @noRd
eq_lattice <- function(u, x, zero_at_origin = TRUE) {
  n <- length(x)
  out <- x[pmin(floor(n * u) + 1L, n)]
  if (zero_at_origin) out[u <= 0] <- 0
  out[u >= 1] <- x[n]
  out
}

#' Poisson probabilities p_k(a), k = 0, ..., m - 1
#' @keywords internal
#' @noRd
pois_w <- function(a, m) stats::dpois(0:(m - 1L), a)

#' Empirical distribution function `mean(x <= t)` at points `t`
#' @keywords internal
#' @noRd
ecdf_at <- function(t, x) {
  n <- length(x)
  findInterval(t, x) / n
}

#' Empirical survival function `mean(x > t)` at points `t`
#' @keywords internal
#' @noRd
esurv_at <- function(t, x) 1 - ecdf_at(t, x)

#' Robust logarithm of Wu (2019)
#'
#' `ln*(x) = log(x)` for `x >= a`, and the linear extension
#' `x / a - 1 + log(a)` below `a`.
#' @keywords internal
#' @noRd
log_star <- function(x, a) ifelse(x >= a, log(pmax(x, .Machine$double.xmin)),
                                  x / a - 1 + log(a))

#' Threshold `a_n` for robust likelihood cross-validation (Wu 2019)
#' @keywords internal
#' @noRd
a_n <- function(x) sqrt(log(length(x)) / 2) / (stats::sd(x) * length(x))

#' Check and sort a sample
#' @keywords internal
#' @noRd
prep_sample <- function(x, positive = FALSE) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  x <- x[is.finite(x)]
  if (length(x) < 5L) stop("`x` must contain at least 5 finite values.", call. = FALSE)
  if (positive && any(x <= 0))
    stop("This estimator requires strictly positive data ",
         "(it smooths on the lattice k/m of the positive half-line). ",
         "Shift the sample or choose another method.", call. = FALSE)
  sort(x)
}

#' Check an evaluation grid
#' @keywords internal
#' @noRd
check_u <- function(u) {
  if (!is.numeric(u) || any(!is.finite(u)) || any(u < 0) || any(u > 1))
    stop("`u` must be numeric with all values in [0, 1].", call. = FALSE)
  u
}

#' Default evaluation grid `i/n`
#' @keywords internal
#' @noRd
default_u <- function(n, scheme = c("i/n", "i/(n+1)")) {
  scheme <- match.arg(scheme)
  if (scheme == "i/n") seq_len(n) / n else seq_len(n) / (n + 1)
}

#' Kernel functions used by the kernel-type estimators
#' @return list with density `d`, distribution function `p` and name.
#' @keywords internal
#' @noRd
kernel_fun <- function(kernel = c("gaussian", "triangular", "epanechnikov")) {
  kernel <- match.arg(kernel)
  switch(kernel,
    gaussian = list(d = stats::dnorm, p = stats::pnorm, name = kernel),
    triangular = list(
      d = function(x) pmax(1 - abs(x), 0),
      p = function(x) ifelse(x <= -1, 0, ifelse(x >= 1, 1,
                       ifelse(x < 0, (1 + x)^2 / 2, 1 - (1 - x)^2 / 2))),
      name = kernel),
    epanechnikov = list(
      d = function(x) 0.75 * pmax(1 - x^2, 0),
      p = function(x) ifelse(x <= -1, 0, ifelse(x >= 1, 1, 0.5 + 0.75 * x - 0.25 * x^3)),
      name = kernel))
}

#' Kaplan-Meier estimate of F at the ordered observations (right censoring)
#'
#' `x` sorted, `status` 1 = event, 0 = censored (aligned with sorted `x`).
#' @keywords internal
#' @noRd
km_cdf <- function(x, status) {
  n <- length(x)
  at_risk <- n:1
  surv <- cumprod(1 - status / at_risk)
  1 - surv
}
