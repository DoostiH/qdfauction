# Competitor density estimators for first-price auctions ----------------------

# R translation of kspdf_bc.m / evalkspdf*.m, (C) Brent R. Hickman and
# Timothy P. Hubbard; see inst/COPYRIGHTS and cite Hickman and Hubbard (2015).
# kernels of Hickman & Hubbard's implementation
hh_kernel <- function(u, kernel) {
  switch(kernel,
         gaussian     = stats::dnorm(u),
         epanechnikov = 0.75 * (1 - u^2) * (abs(u) <= 1),
         triweight    = 35 / 32 * (1 - u^2)^3 * (abs(u) <= 1),
         quartic      = 15 / 16 * (1 - u^2)^2 * (abs(u) <= 1),
         uniform      = 0.5 * (abs(u) <= 1),
         triangle     = (1 - abs(u)) * (abs(u) <= 1))
}
hh_hardle_c <- function(kernel)
  switch(kernel, gaussian = 1, epanechnikov = 2.214, triweight = 2.978,
         uniform = 1.740, triangle = 2.432, quartic = 2.623)

# b0 constant of Karunamuni & Zhang (2008) for the pilot bandwidth h0 = b0 h1
hh_b0 <- function(kernel) {
  num_t1   <- stats::integrate(function(u) u^2 * hh_kernel(u, kernel), -1, 1)$value
  num_t2   <- -(-144 / 5 + 432 / 4 - 468 / 3 + 216 / 2 - 36)
  denom_t1 <- -(-2 + 18 / 4 - 12 / 5)
  denom_t2 <- stats::integrate(function(u) hh_kernel(u, kernel)^2, -1, 1)$value
  ((num_t1^2 * num_t2) / (denom_t1^2 * denom_t2))^(1 / 5)
}

#' Boundary-corrected kernel density estimator of Hickman and Hubbard (2015)
#'
#' @description
#' Kernel density estimator with the transformation-based boundary correction
#' of Karunamuni and Zhang (2008) at the lower endpoint of the support, as
#' implemented by Hickman and Hubbard (2015) for first-price auctions (their
#' `kspdf_bc`; this function is an R translation of that code, whose
#' authors ask that Hickman and Hubbard (2015) be cited). The bandwidth is Silverman's rule scaled by Haerdle's (1991)
#' canonical-kernel constant. The estimate is set to zero outside the range
#' of the data. With `side = "both"` the same correction is applied at the
#' upper endpoint through the reflected sample.
#'
#' @param x evaluation points.
#' @param data observed sample (bids or pseudo-values).
#' @param kernel `"triweight"` (default), `"gaussian"`, `"epanechnikov"`,
#'   `"quartic"`, `"uniform"` or `"triangle"`.
#' @param side `"lower"` (the original implementation) or `"both"`.
#' @param hardle logical; use Haerdle's canonical bandwidth constant.
#' @param h optional bandwidth override.
#' @return a numeric vector with attribute `"h"` (the bandwidth used).
#' @references
#' Hickman, B. R. and Hubbard, T. P. (2015). Replacing sample trimming with
#' boundary correction in nonparametric estimation of first-price auctions.
#' *Journal of Applied Econometrics*, 30, 739--762. \doi{10.1002/jae.2385}
#'
#' Karunamuni, R. J. and Zhang, S. (2008). Some improvements on a boundary
#' corrected kernel density estimator. *Statistics & Probability Letters*, 78,
#' 499--507.
#' @examples
#' x <- rbeta(200, 2, 1)
#' g <- seq(0.01, 0.99, by = 0.01)
#' plot(g, kde_hh(g, x), type = "l"); lines(g, dbeta(g, 2, 1), lty = 2)
#' @export
kde_hh <- function(x, data, kernel = "triweight", side = c("lower", "both"),
                   hardle = TRUE, h = NULL) {
  side <- match.arg(side)
  kernel <- match.arg(kernel, c("triweight", "gaussian", "epanechnikov", "quartic", "uniform", "triangle"))
  if (side == "both") {
    f_lo <- kde_hh(x, data, kernel, "lower", hardle, h)
    f_hi <- kde_hh(-x, -data, kernel, "lower", hardle, h)
    lo <- min(data); hi <- max(data); mid <- (lo + hi) / 2
    out <- ifelse(x <= mid, f_lo, f_hi)
    attr(out, "h") <- attr(f_lo, "h")
    return(out)
  }
  n <- length(data)
  if (is.null(h)) {
    cc <- if (hardle) hh_hardle_c(kernel) else 1
    h <- cc * (4 / 3)^(1 / 5) * stats::sd(data) * n^(-1 / 5)
  }
  lowobs <- min(data); obs <- data - lowobs; xe <- x - lowobs
  h1 <- h * n^(-1 / 20); A <- 0.55
  h0 <- hh_b0(kernel) * h1
  fstar_h1 <- sum(hh_kernel((h1 - obs) / h1, kernel)) / (n * h1) + 1 / n^2
  u0 <- -obs / h0
  fstar_0 <- sum((6 + 18 * u0 + 12 * u0^2) * (u0 >= -1) * (u0 <= 0)) / (n * h0)
  f_0 <- max(fstar_0, 1 / n^2)
  dhat <- (log(fstar_h1) - log(f_0)) / h1
  ghat <- obs + dhat * obs^2 + A * dhat^2 * obs^3
  K1 <- hh_kernel(outer(xe, obs, "-") / h, kernel)
  K2 <- hh_kernel(outer(xe, ghat, "+") / h, kernel)
  fhat <- (rowSums(K1) + rowSums(K2)) / (n * h)
  fhat[xe > max(obs) | xe < min(obs)] <- 0
  attr(fhat, "h") <- h
  fhat
}

# density_marmer_shneyerov(): R translation of qbest4.m, (C) Vadim Marmer and
# Artyom Shneyerov; see inst/COPYRIGHTS.

#' Quantile-based density estimator of Marmer and Shneyerov (2012)
#'
#' Estimates the private-value density in a symmetric IPV first-price auction
#' directly from bids: the bid quantile function and density are estimated
#' with a triweight kernel, the preliminary value quantile
#' \eqn{Q(\tau) + \tau/\{(n-1)\hat g(q(\tau))\}} is monotonised around
#' \eqn{\tau_0 = 0.5}, inverted to give \eqn{\hat F(v)}, and the density is
#' \deqn{\hat f(v) = \Big[\tfrac{n}{n-1}\hat g^{-1} - \tfrac{1}{n-1}\hat F(v)\,\hat g'/\hat g^3\Big]^{-1}}
#' evaluated at \eqn{q(\hat F(v))}. Bandwidths follow the authors'
#' normal-reference rules (\eqn{1.06\,\hat\sigma\,n^{-1/5}} for \eqn{g},
#' \eqn{1.06\,\hat\sigma\,n^{-1/7}} for \eqn{g'}). This function is an R
#' translation of the authors' MATLAB implementation.
#'
#' @param x evaluation points on the value scale.
#' @param bids pooled bids.
#' @param n number of bidders per auction.
#' @param tau0 reference quantile for the monotonisation.
#' @param cap upper cap on the estimate (numerical safeguard).
#' @return numeric vector of density estimates at `x`.
#' @references
#' Marmer, V. and Shneyerov, A. (2012). Quantile-based nonparametric inference
#' for first-price auctions. *Journal of Econometrics*, 167(2), 345--357.
#' @export
density_marmer_shneyerov <- function(x, bids, n, tau0 = 0.5, cap = 1000) {
  b <- sort(as.numeric(bids)); N <- length(b)
  s <- stats::sd(b)
  h_pdf <- 1.06 * N^(-1 / 5) * s; h_prime <- 1.06 * N^(-1 / 7) * s
  tw  <- function(z) 35 / 32 * (1 - z^2)^3 * (abs(z) <= 1)
  twd <- function(z) 105 / 16 * (1 - z^2)^2 * z * (abs(z) <= 1)   # derivative of triweight (sign as in source)
  q <- b; t <- seq_len(N) / N
  gq <- colSums(tw(outer(b, q, "-") / h_pdf)) / (N * h_pdf)
  Q  <- q + t / gq / (n - 1)
  # monotonisation (eq. 13): sup on [tau0, tau], inf on [tau, tau0]
  k0 <- ceiling(N * tau0); Qm <- Q
  Qm[k0:N] <- cummax(Q[k0:N])
  k1 <- floor(N * tau0); Qm[1:k1] <- rev(cummin(rev(Q[1:k1])))
  Fx <- vapply(x, function(v) { i <- which(Qm <= v); if (length(i)) t[max(i)] else 0 }, numeric(1))
  im <- vapply(Fx, function(f) which.min(abs(t - f)), integer(1))
  qF <- q[im]
  gqF  <- colSums(tw(outer(b, qF, "-") / h_pdf)) / (N * h_pdf)
  g1qF <- colSums(twd(outer(b, qF, "-") / h_prime)) / (N * h_prime^2)
  tmp <- (n / (n - 1)) / gqF - (1 / (n - 1)) * Fx * g1qF / gqF^3
  f <- 1 / tmp
  f[!is.finite(f) | f < 0] <- 0
  pmin(f, cap)
}
