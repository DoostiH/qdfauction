# Inference for the GPV estimator ------------------------------------------------
#
# Implemented from Ma, Marmer and Shneyerov (2019), "Inference for first-price
# auctions with Guerre, Perrigne, and Vuong's estimator", Journal of
# Econometrics 211, 507-538, Sections 3-5 and 8. The authors' MATLAB code was
# used to validate this implementation.

# triweight kernel, its derivative, and the fourth-order triweight
.tw   <- function(x) (35 / 32) * (1 - x^2)^3 * (abs(x) <= 1)
.tw_d <- function(x) (105 / 32) * (1 - x^2)^2 * (-2 * x) * (abs(x) <= 1)
.tw4  <- function(x) (27 / 16) * (1 - 11 / 3 * x^2) * .tw(x)

# first stage: empirical cdf, fourth-order kernel density, pseudo-values, trimming
gpv_first_stage <- function(b, n, h_g, b_lo = min(b), b_hi = max(b), trim_estimator = FALSE) {
  N  <- length(b)
  Kg <- .tw4(outer(b, b, "-") / h_g)                 # symmetric N x N
  g  <- rowMeans(Kg) / h_g
  G  <- ecdf_at(b, sort(b))
  inner <- (b >= b_lo + h_g) & (b <= b_hi - h_g)     # GPV trimming region
  use <- g > 0 & (if (trim_estimator) inner else TRUE)
  V <- b + G / ((n - 1) * g)
  V[!use] <- NA_real_
  list(V = V, G = G, g = g, trim = use, inner = inner & g > 0, Kg = Kg)
}

# bootstrap first stage: a resample only reweights the original bids, so the
# kernel matrix of the original sample is reused (w = resampling counts)
gpv_first_stage_w <- function(w, b, n, h_g, Kg, ord, inner0, trim_estimator) {
  N <- length(b)
  g <- as.numeric(Kg %*% w) / (N * h_g)
  cw <- cumsum(w[ord]); G <- numeric(N)
  G[ord] <- cw[findInterval(b[ord], b[ord])] / N      # ties handled by findInterval
  use <- g > 0 & (if (trim_estimator) inner0 else TRUE)
  V <- b + G / ((n - 1) * g)
  V[!use] <- NA_real_
  list(V = V, use = use)
}

# Silverman's robust scale
.rscale <- function(x) { x <- x[is.finite(x)]; min(stats::sd(x), stats::IQR(x) / 1.349) }

# second stage on a grid (trimmed observations contribute zero; denominator N)
gpv_second_stage <- function(V, trim, grid, h_f) {
  Vt <- V[trim]
  colSums(.tw(outer(Vt, grid, "-") / h_f)) / (length(V) * h_f)
}

#' Variance estimator of Ma, Marmer and Shneyerov (2019)
#'
#' Computes the estimator (4.2) of the asymptotic variance of the GPV
#' density estimator,
#' \deqn{\hat V(v) = \frac{1}{N(N-1)^2 h_f^2 h_g}\,\frac{1}{n(n-1)(n-2)}
#'   \sum_i \sum_{j \ne i} \sum_{j' \ne i, j} \eta_{ij}(v)\,\eta_{ij'}(v),}
#' \eqn{\eta_{ij}(v) = T_j K_f'((\hat V_j - v)/h_f)\,\hat G(B_j)/\hat g(B_j)^2\,K_g((B_i - B_j)/h_g)},
#' written with separate first- and second-stage bandwidths. The triple sum
#' over distinct indices is computed exactly with matrix products.
#' @param fs output of the internal first stage.
#' @param grid evaluation points.
#' @param n bidders per auction.
#' @param h_g,h_f bandwidths.
#' @return numeric vector \eqn{\hat V(v)} on `grid`.
#' @keywords internal
#' @noRd
gpv_variance <- function(fs, grid, n, h_g, h_f) {
  N  <- length(fs$G)
  a0 <- ifelse(fs$trim, fs$G / fs$g^2, 0)
  Vt <- ifelse(fs$trim, fs$V, 0)
  A  <- .tw_d(outer(Vt, grid, "-") / h_f) * a0       # N x |grid|, rows j
  A[!fs$trim, ] <- 0
  k0 <- .tw4(0)
  C  <- fs$Kg %*% A                                  # C[i, v] = sum_j k_ij a_j(v)
  S0 <- colSums(C^2)                                 # all (j, j') for each i
  q  <- rowSums(fs$Kg^2)
  S2 <- colSums(A^2 * (q - k0^2))                    # j = j' != i
  S3 <- k0^2 * colSums(A^2)                          # j = j' = i
  S1 <- colSums(k0 * A * (C - k0 * A))               # j = i, j' != i
  U  <- (S0 - 2 * S1 - S2 - S3) / (N * (N - 1) * (N - 2))
  pmax(U, 0) / (n * (n - 1)^2 * h_f^2 * h_g)
}

#' Confidence intervals and uniform confidence bands for the private-value density
#'
#' @description
#' Inference for the two-step estimator of Guerre, Perrigne and Vuong (2000)
#' following Ma, Marmer and Shneyerov (2019). The first stage uses the
#' empirical distribution function of the bids and a fourth-order triweight
#' kernel density estimate with bandwidth \eqn{h_g = 3.72\,\hat\sigma_b\,n^{-1/5}};
#' bids within \eqn{h_g} of the sample extremes are trimmed. The second stage
#' is a triweight kernel density estimate of the pseudo-values with bandwidth
#' \eqn{h_f = 3.15\,\hat\sigma_v\,n_T^{-1/5}}, where \eqn{\hat\sigma_v} and
#' \eqn{n_T} are the standard deviation and number of the trimmed pseudo-values
#' (Section 8 of the paper). The trimmed pseudo-values are those of bids at
#' least \eqn{h_g} from the sample extremes. As scale estimate we use
#' Silverman's robust \eqn{\min(\hat\sigma, \mathrm{IQR}/1.349)} rather than
#' the standard deviation alone; the two coincide in the paper's designs, but
#' the robust version protects the bandwidths against the few very large
#' pseudo-values that arise where the estimated bid density is small.
#'
#' Standard errors use the estimator (4.2) of the asymptotic variance, which
#' accounts for the estimation error of the first stage. Pointwise intervals
#' are percentile bootstrap intervals (Theorem 4.2; the default) or normal
#' intervals \eqn{\hat f \pm z_{1-\alpha/2}\,\mathrm{se}}. The analytic
#' standard error contains only the leading, first-stage term of the variance;
#' in moderate samples, and where \eqn{F(v)} is small, the second-stage kernel
#' variance is not negligible and normal intervals can under-cover, whereas the
#' percentile intervals and the uniform band, whose critical value is
#' bootstrapped, are not affected. The uniform band over `grid` is
#' \eqn{\hat f \pm \zeta^*\,\mathrm{se}}, where \eqn{\zeta^*} is the
#' \eqn{1-\alpha} quantile of \eqn{\sup_v |\hat f^*(v) - \hat f(v)|/\mathrm{se}(v)}
#' over bootstrap samples of the pooled bids (Theorem 5.2 and Corollary 5.1 of
#' the paper).
#'
#' The theory requires `grid` to lie in the interior of the support of the
#' valuations. By default it spans the 5th to 95th percentiles of the
#' pseudo-values of the bids between the 10th and 90th bid percentiles, which
#' keeps it away from the long right tail of pseudo-values that the first
#' stage produces where the bid density is small; with real data, choose the
#' range with care, since the band is only valid where the density of values is
#' bounded away from zero. Observations whose estimated bid density is not positive
#' (possible with a fourth-order kernel near the boundary) are dropped.
#'
#' @param bids pooled bids (vector) or an `L` by `n` matrix.
#' @param n number of bidders per auction (from `ncol(bids)` for a matrix).
#' @param grid evaluation points; default 201 points between the 10th and
#'   90th percentiles of the trimmed pseudo-values.
#' @param level confidence level.
#' @section Reproducibility:
#' Fold assignment (or bootstrap resampling) uses R's random number generator;
#' call [set.seed()] beforehand to make results exactly reproducible.
#' @param B number of bootstrap replications.
#' @param pointwise `"normal"` or `"percentile"`.
#' @param h_g,h_f optional bandwidths overriding the rules above.
#' @param constants rule-of-thumb constants for \eqn{h_g} and \eqn{h_f}
#'   (default `c(3.72, 3.15)`, the Silverman constants for the fourth- and
#'   second-order triweight kernels, as in Section 8 of the paper).
#' @param sd_v standard deviation used in \eqn{h_f}: of the `"trimmed"`
#'   pseudo-values with \eqn{n_T} (default, as in the paper) or of `"all"`
#'   pseudo-values with \eqn{n}.
#' @param se `"augmented"` (default): the asymptotic variance (4.2) plus the
#'   second-stage kernel variance \eqn{\hat f(v)\int K^2/(n h_f)}, a term of
#'   smaller order that is asymptotically negligible but keeps the standard
#'   error positive and improves it in moderate samples, where the U-statistic
#'   (4.2) can be close to or below zero; `"asymptotic"`: the estimator (4.2)
#'   alone, as in the paper. Both give asymptotically valid uniform bands.
#' @param trim_estimator logical. The trimming of bids within \eqn{h_g} of the
#'   sample extremes is always used to compute \eqn{h_f}. With the default
#'   `FALSE` all pseudo-values then enter the density estimate and its
#'   variance, as in the authors' Monte Carlo implementation; with `TRUE` the
#'   trimmed bids are also excluded from the estimate (the trimming factor
#'   \eqn{T_{il}} of the paper), which with rule-of-thumb bandwidths removes a
#'   substantial part of the support in moderate samples.
#' @return An object of class `"fpa_density_ci"`: a list with `v` (grid),
#'   `f` (estimate), `se`, `lower`, `upper` (pointwise interval), `band_lower`,
#'   `band_upper` (uniform band), `crit` (\eqn{\zeta^*}), `level`, `h_g`, `h_f`,
#'   `values` (trimmed pseudo-values) and `B`.
#' @references
#' Ma, J., Marmer, V. and Shneyerov, A. (2019). Inference for first-price
#' auctions with Guerre, Perrigne, and Vuong's estimator. *Journal of
#' Econometrics*, 211, 507--538. \doi{10.1016/j.jeconom.2019.02.006}
#'
#' Guerre, E., Perrigne, I. and Vuong, Q. (2000). Optimal nonparametric
#' estimation of first-price auctions. *Econometrica*, 68, 525--574.
#' @examples
#' set.seed(1)
#' sim <- fpa_simulate(100, 5, Q = function(u) u)       # uniform values
#' ci  <- fpa_density_ci(sim$bids, grid = seq(0.2, 0.8, by = 0.01), B = 99)
#' ci
#' plot(ci, true = function(v) dunif(v))
#' @export
fpa_density_ci <- function(bids, n = NULL, grid = NULL, level = 0.95, B = 499,
                           pointwise = c("percentile", "normal"),
                           h_g = NULL, h_f = NULL, trim_estimator = FALSE,
                           constants = c(3.72, 3.15), sd_v = c("trimmed", "all"),
                           se = c("augmented", "asymptotic")) {
  sd_v <- match.arg(sd_v); se_type <- match.arg(se)
  pointwise <- match.arg(pointwise)
  bm <- as_bid_matrix(bids, n); n <- bm$n
  b <- as.numeric(bm$bids); N <- length(b)
  if (N < 20) stop("At least 20 bids are needed.", call. = FALSE)
  L <- N / n
  if (is.null(h_g)) h_g <- constants[1] * .rscale(b) * N^(-1 / 5)
  lo <- min(b); hi <- max(b)
  fs <- gpv_first_stage(b, n, h_g, lo, hi, trim_estimator)
  Vt <- fs$V[fs$inner]                               # trimmed pseudo-values
  if (length(Vt) < 10) stop("Too few bids remain after trimming.", call. = FALSE)
  if (is.null(h_f)) h_f <- if (sd_v == "trimmed") constants[2] * .rscale(Vt) * length(Vt)^(-1 / 5)
                           else constants[2] * .rscale(fs$V[fs$trim]) * N^(-1 / 5)
  if (is.null(grid)) {
    # pseudo-values of the central bids: the interior region where the theory applies
    qb <- stats::quantile(b, c(0.10, 0.90), names = FALSE)
    Vc <- fs$V[b >= qb[1] & b <= qb[2] & fs$trim]
    r  <- stats::quantile(Vc, c(0.05, 0.95), names = FALSE)
    grid <- seq(r[1], r[2], length.out = 201)
  }
  f  <- gpv_second_stage(fs$V, fs$trim, grid, h_f)
  Vh <- gpv_variance(fs, grid, n, h_g, h_f)
  se <- sqrt(Vh / (L * h_f^2 * h_g))
  if (se_type == "augmented")                         # + second-stage kernel variance
    se <- sqrt(se^2 + pmax(f, 0) * (350 / 429) / (N * h_f))
  se_safe <- ifelse(se > 0, se, NA_real_)

  boot <- matrix(NA_real_, B, length(grid))
  ord <- order(b); inner0 <- (b >= lo + h_g) & (b <= hi - h_g)
  for (r in seq_len(B)) {
    w  <- tabulate(sample.int(N, N, replace = TRUE), N)
    fb <- gpv_first_stage_w(w, b, n, h_g, fs$Kg, ord, inner0, trim_estimator)
    k  <- which(fb$use & w > 0)
    boot[r, ] <- colSums(w[k] * .tw(outer(fb$V[k], grid, "-") / h_f)) / (N * h_f)
  }
  Z <- sweep(boot, 2, f) / rep(se_safe, each = B)
  supZ <- apply(abs(Z), 1, function(z) if (all(is.na(z))) NA else max(z, na.rm = TRUE))
  alpha <- 1 - level
  crit <- stats::quantile(supZ, level, na.rm = TRUE, names = FALSE)
  if (pointwise == "normal") {
    z <- stats::qnorm(1 - alpha / 2)
    lower <- f - z * se; upper <- f + z * se
  } else {
    lower <- apply(boot, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE)
    upper <- apply(boot, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  }
  structure(list(v = grid, f = f, se = se, lower = lower, upper = upper,
                 band_lower = f - crit * se, band_upper = f + crit * se, crit = crit,
                 level = level, pointwise = pointwise, se_type = se_type, h_g = h_g, h_f = h_f,
                 values = Vt, n = n, N = N, B = B, call = match.call()),
            class = "fpa_density_ci")
}

#' @export
print.fpa_density_ci <- function(x, digits = 4, ...) {
  cat("GPV estimate of the private-value density with confidence bands\n")
  cat("  (Ma, Marmer and Shneyerov 2019)\n")
  cat("  bids        :", x$N, "pooled,", x$n, "bidders per auction\n")
  cat("  bandwidths  : h_g =", format(x$h_g, digits = digits), " h_f =", format(x$h_f, digits = digits), "\n")
  cat("  grid        :", format(min(x$v), digits = digits), "to", format(max(x$v), digits = digits),
      "(", length(x$v), "points )\n")
  cat("  level       :", x$level, "; pointwise intervals:", x$pointwise, "\n")
  cat("  uniform band: critical value", format(x$crit, digits = digits), "from", x$B,
      "bootstrap samples\n")
  cat("                max width", format(max(x$band_upper - x$band_lower, na.rm = TRUE), digits = digits), "\n")
  invisible(x)
}

#' Plot a GPV density estimate with confidence bands
#' @param x an `"fpa_density_ci"` object.
#' @param true optional true density function for simulated data.
#' @param ... passed to [graphics::plot()].
#' @return No return value, called for its side effect of drawing the
#'   density estimate with its pointwise confidence interval and uniform
#'   confidence band. The object `x` is returned invisibly.
#' @export
plot.fpa_density_ci <- function(x, true = NULL, ...) {
  yl <- range(c(x$band_lower, x$band_upper, x$lower, x$upper, 0), na.rm = TRUE)
  args <- list(x = x$v, y = x$f, type = "n", ylim = yl, xlab = "private value", ylab = "density",
               main = sprintf("GPV density with %g%% uniform band", 100 * x$level))
  args <- utils::modifyList(args, list(...))
  do.call(graphics::plot, args)
  graphics::polygon(c(x$v, rev(x$v)), c(x$band_lower, rev(x$band_upper)),
                    col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
  graphics::lines(x$v, x$lower, lty = 2); graphics::lines(x$v, x$upper, lty = 2)
  graphics::lines(x$v, x$f, lwd = 2)
  leg <- c("estimate", "pointwise interval", "uniform band"); lt <- c(1, 2, NA); cl <- c(1, 1, NA)
  fl <- c(NA, NA, grDevices::adjustcolor("steelblue", 0.25))
  if (!is.null(true)) { graphics::lines(x$v, true(x$v), col = 2, lwd = 2, lty = 3)
    leg <- c(leg, "true"); lt <- c(lt, 3); cl <- c(cl, 2); fl <- c(fl, NA) }
  graphics::legend("topright", leg, lty = lt, col = cl, fill = fl, border = NA, bty = "n")
  invisible(x)
}
