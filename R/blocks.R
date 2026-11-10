## =====================================================================
## R/blocks.R -- moving-block bootstrap over auctions for qdfauction 0.3.0
##
## Three exported functions plus one internal helper:
##
##   block_length_auto(B, mult = 2, warn_below = 15)
##        L* = ceiling(mult * L_pw), with L_pw the Politis-White (2004) /
##        Patton-Politis-White (2009) automatic block length computed on the
##        per-auction winning bid max_i B[t, i]. Returns an integer with
##        attributes L_pw and T_over_Lpw; warns when T / L_pw < warn_below.
##
##   resample_blocks(TT, L)
##        Row indices of a circular moving-block resample of TT auctions with
##        block length L. L = 1 is the i.i.d. auction bootstrap.
##
##   fpa_block_boot(B, stat, L = "auto", B_boot = 199, level = 0.90,
##                  method = "bernstein", bandwidth = "bcv", monotone = "spline",
##                  ...)
##        Generic block bootstrap of any statistic of an fpa_values() fit:
##        selects the bandwidth once on B, holds it fixed, refits on each block
##        resample, applies stat(fit) and returns point estimate, bootstrap
##        s.d., percentile interval, L used and the diagnostics.
##        Two convenience wrappers:
##          fpa_revenue_boot_blocks(B, r, v_s, ...)   stat = fpa_revenue(fit, r, v_s)
##          fpa_quantile_boot_blocks(B, u, ...)       stat = stats::quantile(fit$values, u)
##
## Integration into 0.3.0: add resample = c("auctions", "bids", "blocks") and
## L = "auto" to fpa_revenue_boot() and fpa_density_boot(); when
## resample == "blocks", obtain L via block_length_auto() (or the user's
## integer) and replace the row-sampling line with resample_blocks(TT, L).
## Everything else in those functions (bandwidth held fixed, refits,
## percentile intervals) is already in place.
##
## Dependencies: np (Suggests; required only for L = "auto").
## =====================================================================

#' Automatic block length for the auction bootstrap
#'
#' @param B numeric matrix of bids, auctions in rows in time order, bidders
#'   in columns.
#' @param mult multiplier on the Politis-White length; 2 is the value
#'   validated in the paper, results are insensitive over 1.5 to 3.
#' @param winner "max" for a selling auction (highest bid wins), "min" for a
#'   procurement auction (lowest bid wins); the winning-bid series is built
#'   accordingly. Pass procurement bids in their original orientation with
#'   winner = "min"; the fitting itself is done by proc_values(), which
#'   reflects them.
#' @param warn_below warn when T / L_pw is below this value: the sample is
#'   then too short for the dependence range it contains and no block length
#'   recovers the full sampling variance (the paper's Table 6).
#' @return integer block length with attributes "L_pw" and "T_over_Lpw".
#' @export
block_length_auto <- function(B, mult = 2, warn_below = 15, winner = c("max", "min")) {
  winner <- match.arg(winner)
  if (!requireNamespace("np", quietly = TRUE))
    stop("block_length_auto() needs package 'np' for the Politis-White length; ",
         "install.packages('np') or pass an integer L.", call. = FALSE)
  B <- as.matrix(B)
  TT <- nrow(B)
  if (TT < 10L) stop("at least 10 auctions are needed for a block length.", call. = FALSE)
  w <- apply(B, 1, winner, na.rm = TRUE)              # per-auction winning bid: max (selling) or min (procurement)
  bs <- tryCatch(np::b.star(w, round = TRUE), error = function(e) NULL)
  L_pw <- if (is.null(bs)) 1 else max(1, as.numeric(bs[1, "BstarCB"]))
  L <- as.integer(ceiling(mult * L_pw))
  L <- max(1L, min(L, TT %/% 2L))                     # never fewer than two blocks
  ratio <- TT / L_pw
  if (ratio < warn_below)
    warning(sprintf(paste0("T / L_pw = %.1f is below %g: the sample is short relative to the ",
                           "dependence range, and the block bootstrap interval is a lower bound ",
                           "on the sampling uncertainty."),
                    ratio, warn_below), call. = FALSE)
  structure(L, L_pw = L_pw, T_over_Lpw = ratio)
}

#' Circular moving-block resample of auction indices
#' @param TT number of auctions.
#' @param L block length; 1 reproduces the i.i.d. auction bootstrap.
#' @return integer vector of length TT of row indices.
#' @export
resample_blocks <- function(TT, L) {
  ## TT must be an integer too: with a double TT (e.g. resample_blocks(100, 10))
  ## the modulo below returns doubles and vapply(FUN.VALUE = integer(L)) fails.
  TT <- as.integer(TT)
  L  <- as.integer(L)
  if (is.na(TT) || TT < 1L) stop("`TT` must be a positive integer.", call. = FALSE)
  if (L <= 1L) return(sample.int(TT, TT, replace = TRUE))
  nb <- ceiling(TT / L)
  st <- sample.int(TT, nb, replace = TRUE)
  idx <- as.vector(vapply(st, function(s) ((s - 1L) + seq_len(L) - 1L) %% TT + 1L, integer(L)))
  idx[seq_len(TT)]
}

#' Block bootstrap of a statistic of the recovered values
#'
#' @param B bid matrix, auctions in rows in time order.
#' @param stat function of an fpa_values() fit returning a numeric vector.
#' @param L "auto" for block_length_auto(B), or an integer; 1 gives the
#'   i.i.d. auction bootstrap.
#' @param B_boot number of bootstrap resamples.
#' @param level nominal level of the percentile interval.
#' @param method,bandwidth,monotone passed to fpa_values(); the bandwidth is
#'   selected once on the full sample and held fixed in every refit.
#' @param winner which bid identifies the auction for the block-length rule:
#'   "max" for a sale auction, "min" for procurement.
#' @param v_s seller's valuation, passed on to the statistic where it applies.
#' @param ... further arguments to fpa_values().
#' @return list with point, boot_sd, lower, upper, L, L_pw, T_over_Lpw,
#'   boot (the B_boot x length(stat) matrix of bootstrap statistics).
#' @export
fpa_block_boot <- function(B, stat, L = "auto", B_boot = 199L, level = 0.90,
                           method = "bernstein", bandwidth = "bcv",
                           monotone = "spline", winner = "max", ...) {
  B <- as.matrix(B); TT <- nrow(B)
  ## bandwidth once, on the pooled sample
  h <- if (is.character(bandwidth)) qdf(as.numeric(B), method = method, bandwidth = bandwidth)$h
       else bandwidth
  fit0 <- fpa_values(B, method = method, bandwidth = h, monotone = monotone, ...)
  point <- as.numeric(stat(fit0))
  if (identical(L, "auto")) L <- block_length_auto(B, winner = winner)
  L_pw <- attr(L, "L_pw"); ratio <- attr(L, "T_over_Lpw"); L <- as.integer(L)
  boot <- matrix(NA_real_, B_boot, length(point))
  for (b in seq_len(B_boot)) {
    f <- tryCatch(fpa_values(B[resample_blocks(TT, L), , drop = FALSE], method = method,
                             bandwidth = h, monotone = monotone, ...),
                  error = function(e) NULL)
    if (!is.null(f)) boot[b, ] <- as.numeric(stat(f))
  }
  a <- (1 - level) / 2
  list(point   = point,
       boot_sd = apply(boot, 2, stats::sd, na.rm = TRUE),
       lower   = apply(boot, 2, stats::quantile, probs = a,     na.rm = TRUE, names = FALSE),
       upper   = apply(boot, 2, stats::quantile, probs = 1 - a, na.rm = TRUE, names = FALSE),
       L = L, L_pw = L_pw, T_over_Lpw = ratio, bandwidth = h, level = level,
       failures = sum(is.na(boot[, 1])), boot = boot)
}

#' @rdname fpa_block_boot
#' @param r reserve price(s); @param v_s seller valuation.
#' @export
fpa_revenue_boot_blocks <- function(B, r, v_s = 0, L = "auto", ...) {
  out <- fpa_block_boot(B, stat = function(f) fpa_revenue(f, r, v_s = v_s), L = L, ...)
  out$r <- r; out
}

#' @rdname fpa_block_boot
#' @param u probabilities at which the recovered value quantile is evaluated.
#' @export
fpa_quantile_boot_blocks <- function(B, u = c(0.25, 0.5, 0.75), L = "auto", ...) {
  out <- fpa_block_boot(B, stat = function(f) stats::quantile(f$values, u, type = 6, names = FALSE),
                        L = L, ...)
  out$u <- u; out
}

## ---- minimal check (not run on load) ---------------------------------
## library(qdfauction); set.seed(1)
## Q <- function(u) qgamma(u, 10, 2); q <- function(u) 1/dgamma(qgamma(u,10,2),10,2)
## bf <- fpa_bid_function(Q, q, 5); beta <- function(u) stats::approx(bf$U, bf$bid, xout = u, rule = 2)$y
## Z <- matrix(0, 200, 5); Z[1,] <- stats::rnorm(5); for (t in 2:200) Z[t,] <- 0.7*Z[t-1,] + sqrt(1-0.49)*stats::rnorm(5)
## B <- matrix(beta(as.numeric(stats::pnorm(Z))), 200, 5)
## block_length_auto(B)                                  # L_pw about 8, L* about 16
## fpa_revenue_boot_blocks(B, r = c(4, 5, 6), v_s = 2)[c("point", "boot_sd", "L")]
## fpa_revenue_boot_blocks(B, r = 4, v_s = 2, L = 1)$boot_sd   # i.i.d. for comparison: about 0.55 of the above
