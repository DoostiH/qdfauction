# Periodised orthonormal wavelet basis on [0, 1] --------------------------------

#' Mother wavelet and periodised wavelet basis on the unit interval
#'
#' @description
#' `psi_jk()` evaluates the mother wavelet
#' \eqn{\psi_{j,k}(z) = 2^{j/2}\psi(2^j z - k)}, with \eqn{\psi} obtained from
#' the scaling function by the two-scale relation
#' \eqn{\psi(x) = \sqrt 2 \sum_k g_k \phi(2x - k)},
#' \eqn{g_k = (-1)^k h_{N-k}}. `phi_per()` and `psi_per()` evaluate the
#' periodised functions \eqn{\sum_l \phi_{j,k}(z + l)} and
#' \eqn{\sum_l \psi_{j,k}(z + l)}, which for \eqn{2^j \ge N} form an
#' orthonormal basis of \eqn{L_2([0,1])} (Daubechies 1992, Section 9.3); they
#' are used by the nonlinear estimators [qdf_wavelet_nl()] and
#' `fpa_density(method = "wavelet_theory")`.
#'
#' @param z numeric vector in the unit interval.
#' @param j resolution level.
#' @param k shift(s); for more than one, a matrix `length(z)` by `length(k)`.
#' @param filter scaling filter; see [wavelet_filter()].
#' @return numeric vector or matrix.
#' @examples
#' z <- seq(0, 1, length.out = 2001)
#' P <- psi_per(z, 4, 0:15, "db3")
#' round(crossprod(P)[1:3, 1:3] / 2000, 3)     # approximately the identity
#' @export
psi_jk <- function(z, j, k, filter = "db3") {
  h <- wavelet_filter(filter); N <- length(h) - 1L
  g <- (-1)^(0:N) * rev(h)
  ev <- function(x) {
    out <- numeric(length(x))
    for (m in 0:N) out <- out + sqrt(2) * g[m + 1L] * phi_jk(2 * x - m, 0, 0, filter)
    out
  }
  if (length(k) == 1L) return(2^(j / 2) * ev(2^j * z - k))
  2^(j / 2) * vapply(k, function(kk) ev(2^j * z - kk), numeric(length(z)))
}

.periodise <- function(fun, z, j, k, filter) {
  h <- wavelet_filter(filter); N <- length(h) - 1L
  lmax <- ceiling(N / 2^j) + 1L
  out <- 0
  for (l in -lmax:lmax) out <- out + fun(z + l, j, k, filter)
  out
}

#' @rdname psi_jk
#' @export
phi_per <- function(z, j, k, filter = "db3") .periodise(phi_jk, z, j, k, filter)

#' @rdname psi_jk
#' @export
psi_per <- function(z, j, k, filter = "db3") .periodise(psi_jk, z, j, k, filter)

# sup |psi'| (the constant K of Chesneau, Dewan and Doosti 2016, Theorem 3.2)
.psi_deriv_sup <- function(filter = "db3") {
  x <- seq(0, length(wavelet_filter(filter)) - 1, length.out = 20001)
  y <- psi_jk(x, 0, 0, filter)
  max(abs(diff(y) / diff(x)))
}
