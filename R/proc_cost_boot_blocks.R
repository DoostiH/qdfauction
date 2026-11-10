## =====================================================================
## R/proc_cost_boot_blocks.R -- block option for Project A's proc_cost_boot()
##
## Replaces proc_cost_boot() in Project A's R/procurement.R. The only change
## is the resampling line: with resample = "blocks" whole auctions are still
## resampled together (so the within-auction copula structure is kept), but
## in circular blocks of consecutive auctions, so dependence ACROSS auctions
## is kept too. Bandwidth and copula parameter stay at full-sample values,
## as in Project A. resample = "auctions" (default) reproduces Project A's
## function exactly, draw for draw, given the seed.
##
## Needs resample_blocks() and block_length_auto() from R/blocks.R.
## Rows of `bids` must be in time order for "blocks"; the package cannot
## check this, so it is the user's choice, never the default.
## =====================================================================

#' Bootstrap standard errors for expected procurement cost
#'
#' @inheritParams proc_values
#' @param r Price ceiling(s); \code{Inf} for no ceiling.
#' @param v0 Buyer's outside option.
#' @param B Bootstrap replications.
#' @param method First-stage quantile density estimator, passed to
#'   \code{\link{proc_values}}.
#' @param bandwidth First-stage bandwidth rule or numeric bandwidth. Whatever
#'   is given, the bandwidth is selected once on the full sample and held
#'   fixed in every bootstrap refit.
#' @param family Dependence specification: \code{"independence"},
#'   \code{"frank"} or \code{"nonparametric"} (see \code{\link{proc_values}}).
#' @param resample \code{"auctions"} (i.i.d. auction bootstrap, default) or
#'   \code{"blocks"} (moving-block bootstrap over consecutive auctions; use
#'   when the same firms bid repeatedly over time and the rows are in time
#'   order).
#' @param L Block length for \code{"blocks"}: \code{"auto"} (the rule of
#'   \code{\link{block_length_auto}}, applied with \code{winner = "min"}) or
#'   an integer. \code{L = 1} is the i.i.d. auction bootstrap.
#' @return List: \code{r}, \code{estimate}, \code{se}, \code{ci90},
#'   \code{boot}, \code{h}, \code{resample}, \code{L}.
#' @export
proc_cost_boot <- function(bids, r = Inf, v0, B = 499, method = "bernstein",
                           bandwidth = "sqrt", family = "independence",
                           resample = c("auctions", "blocks"), L = "auto", ...) {
  stopifnot(is.matrix(bids))
  resample <- match.arg(resample)
  TT <- nrow(bids)
  if (resample == "blocks") {
    if (identical(L, "auto")) L <- block_length_auto(bids, winner = "min")
    L <- as.integer(L)
    if (L < 1L || L > TT) stop("`L` must be between 1 and the number of auctions.", call. = FALSE)
  } else L <- 1L

  pv0 <- proc_values(bids, method = method, bandwidth = bandwidth, family = family, ...)
  h0  <- pv0$reflected$qdf$h; if (is.null(h0)) h0 <- bandwidth
  th0 <- pv0$reflected$copula$theta
  est <- proc_cost(pv0, r, v0)
  boot <- matrix(NA_real_, B, length(r))
  for (b in seq_len(B)) {
    idx <- if (resample == "blocks") resample_blocks(TT, L) else sample(TT, replace = TRUE)
    Bb  <- bids[idx, , drop = FALSE]
    pvb <- if (is.null(th0))
      proc_values(Bb, method = method, bandwidth = h0, family = family, ...)
    else proc_values(Bb, method = method, bandwidth = h0, family = family, theta = th0, ...)
    boot[b, ] <- proc_cost(pvb, r, v0)
  }
  list(r = r, estimate = est, se = apply(boot, 2, stats::sd),
       ci90 = apply(boot, 2, stats::quantile, c(.05, .95)), boot = boot, h = h0,
       resample = resample, L = L)
}
