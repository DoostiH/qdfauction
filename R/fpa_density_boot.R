#' Bootstrap confidence intervals and bands for any private-value density estimator
#'
#' @description
#' Pointwise percentile intervals and a uniform confidence band for the
#' private-value density estimated by any [fpa_density()] method, for example
#' the quantile-based estimator of Marmer and Shneyerov (2012) or the wavelet
#' estimators. The whole two-step procedure (first stage, unit map, second
#' stage) is repeated on `B` bootstrap samples. The uniform band is
#' \eqn{\hat f(x) \pm c^*\,\hat\sigma^*(x)}, with \eqn{\hat\sigma^*} the
#' bootstrap standard deviation and \eqn{c^*} the \eqn{1-\alpha} quantile of
#' \eqn{\sup_x |\hat f^*(x) - \hat f(x)|/\hat\sigma^*(x)}. For the GPV estimator,
#' [fpa_density_ci()] gives the analytic studentisation of Ma, Marmer and
#' Shneyerov (2019).
#'
#' The percentile bootstrap reproduces the variability of the estimator but not
#' its smoothing bias, so at the usual bandwidths the intervals can under-cover
#' in moderate samples (in a small check with the Marmer--Shneyerov estimator,
#' 240 bids and nominal level 0.95, pointwise coverage was about 0.86).
#' Undersmoothing, i.e. a smaller bandwidth than is optimal for estimation,
#' is the usual remedy.
#'
#' @inheritParams fpa_density
#' @param B number of bootstrap replications.
#' @param level confidence level.
#' @param resample `"bids"` (pooled bids with replacement) or `"auctions"`
#'   (whole auctions; requires a bid matrix).
#' @param ... further arguments to [fpa_density()] (e.g. `j0`, `h`, `first_stage`).
#' @return an object of class `"fpa_density_boot"` with `x`, `f`, `se`,
#'   `lower`, `upper` (pointwise), `band_lower`, `band_upper`, `crit`, `level`,
#'   `method`, `B` and the matrix `boot` of bootstrap estimates.
#' @section Reproducibility:
#' Call [set.seed()] beforehand to make the bootstrap reproducible.
#' @examples
#' \donttest{
#' set.seed(1)
#' sim <- fpa_simulate(60, 4, Q = function(u) qbeta(u, 2, 2))
#' bb <- fpa_density_boot(sim$bids, method = "marmer_shneyerov",
#'                        x = seq(0.15, 0.85, by = 0.05), B = 49)
#' bb
#' plot(bb, true = function(v) dbeta(v, 2, 2))
#' }
#' @export
fpa_density_boot <- function(bids, n = NULL, method = "marmer_shneyerov", x = NULL, B = 199,
                             level = 0.95, resample = c("bids", "auctions"), ...) {
  resample <- match.arg(resample)
  bm <- as_bid_matrix(bids, n); n <- bm$n
  fit <- fpa_density(bm$bids, n = n, method = method, x = x, ...)
  x <- fit$x
  if (resample == "auctions" && !bm$is_mat)
    stop("`resample = \"auctions\"` needs the bids as a matrix.", call. = FALSE)
  boot <- matrix(NA_real_, B, length(x))
  for (r in seq_len(B)) {
    bs <- if (resample == "auctions") bm$bids[sample.int(nrow(bm$bids), replace = TRUE), , drop = FALSE]
          else sample(as.numeric(bm$bids), replace = TRUE)
    boot[r, ] <- tryCatch(fpa_density(bs, n = n, method = method, x = x, ...)$f,
                          error = function(e) rep(NA_real_, length(x)))
  }
  alpha <- 1 - level
  se <- apply(boot, 2, stats::sd, na.rm = TRUE)
  se_ok <- ifelse(se > 0, se, NA_real_)
  Z <- abs(sweep(boot, 2, fit$f)) / rep(se_ok, each = B)
  sup <- apply(Z, 1, function(z) if (all(is.na(z))) NA else max(z, na.rm = TRUE))
  crit <- stats::quantile(sup, level, na.rm = TRUE, names = FALSE)
  structure(list(x = x, f = fit$f, se = se,
                 lower = apply(boot, 2, stats::quantile, alpha / 2, na.rm = TRUE),
                 upper = apply(boot, 2, stats::quantile, 1 - alpha / 2, na.rm = TRUE),
                 band_lower = fit$f - crit * se, band_upper = fit$f + crit * se,
                 crit = crit, level = level, method = method, B = B, resample = resample,
                 boot = boot, fit = fit), class = "fpa_density_boot")
}

#' @export
print.fpa_density_boot <- function(x, digits = 4, ...) {
  cat("Bootstrap inference for the private-value density\n")
  cat("  method       :", x$method, "\n")
  cat("  bootstrap    :", x$B, "replications, resampling", x$resample, "\n")
  cat("  grid         :", format(min(x$x), digits = digits), "to", format(max(x$x), digits = digits),
      "(", length(x$x), "points )\n")
  cat("  uniform band : level", x$level, ", critical value", format(x$crit, digits = digits), "\n")
  invisible(x)
}

#' @rdname fpa_density_boot
#' @param true optional true density function (simulated data).
#' @export
plot.fpa_density_boot <- function(x, true = NULL, ...) {
  yl <- range(c(x$band_lower, x$band_upper, x$lower, x$upper, 0), na.rm = TRUE)
  args <- list(x = x$x, y = x$f, type = "n", ylim = yl, xlab = "private value", ylab = "density",
               main = sprintf("%s with %g%% bootstrap band", x$method, 100 * x$level))
  args <- utils::modifyList(args, list(...))
  do.call(graphics::plot, args)
  graphics::polygon(c(x$x, rev(x$x)), c(x$band_lower, rev(x$band_upper)),
                    col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
  graphics::lines(x$x, x$lower, lty = 2); graphics::lines(x$x, x$upper, lty = 2)
  graphics::lines(x$x, x$f, lwd = 2)
  if (!is.null(true)) graphics::lines(x$x, true(x$x), col = 2, lty = 3, lwd = 2)
  invisible(x)
}
