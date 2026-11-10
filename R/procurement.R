# Procurement (lowest-bid-wins) auctions by reflection.
# Ported from the Project A pipeline (run against qdfauction 0.2.0, R 4.6,
# Nov 2026). See HANDOVER_procurement_to_0.3.0.md for evidence.

.proc_families <- c("independence", "frank", "nonparametric")

.proc_check_family <- function(dots) {
  fam <- dots$family
  if (is.null(fam)) return(invisible(TRUE))
  if (fam %in% c("clayton", "gumbel"))
    stop("Under reflection the bid copula becomes its survival copula; ",
         "for family = '", fam, "' the 180-degree rotated family is required, ",
         "which is not yet implemented. Use 'frank' (radially symmetric) or ",
         "'nonparametric'.", call. = FALSE)
  if (!fam %in% .proc_families)
    stop("family must be one of ", paste(.proc_families, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

#' Reflection constant for procurement bids
#'
#' Returns \eqn{M = \max b + \delta\,\mathrm{range}(b)}. With this choice the
#' reflected sample \eqn{M - b} starts at about zero, so the
#' \eqn{Q_n(0) = 0} anchor of the series quantile density estimators
#' (\code{"bernstein"}, \code{"poisson"}, \code{"indirect_poisson"}) is
#' correct to \eqn{O(\delta)}.
#'
#' @param bids Numeric vector or \eqn{L \times n} matrix of procurement bids.
#' @param delta Relative offset above the largest bid.
#' @return A single number.
#' @export
proc_M <- function(bids, delta = 1e-6) {
  b <- as.numeric(bids)
  stopifnot(all(is.finite(b)), delta > 0)
  max(b) + delta * diff(range(b))
}

#' Recover private costs from procurement (lowest-bid-wins) auctions
#'
#' Reflects the bids, \eqn{b' = M - b}, calls \code{\link{fpa_values}} on the
#' reflected sample, and maps the result back: costs \eqn{c = M - v'}, ranks
#' \eqn{u = 1 - u'}, quantile density \eqn{q(u) = q'(1-u)}. This implements
#' \eqn{c = b - (1-\eta)(1-u)q(u)/(n-1)} under independence and the copula
#' analogue otherwise.
#'
#' @param bids \eqn{L \times n} matrix of bids (one row per auction), or a
#'   vector with \code{n} given.
#' @param n Number of bidders; taken from \code{ncol(bids)} for a matrix.
#' @param M Reflection constant; default \code{proc_M(bids)}.
#' @param monotone Monotonisation of the pseudo-values; \code{"isotonic"} by
#'   default, because the spline replacement extrapolates into the lower
#'   (winning) tail of procurement bids.
#' @param ... Passed to \code{\link{fpa_values}} (\code{method},
#'   \code{bandwidth}, \code{family}, \code{theta}, \code{eta}, ...).
#'   \code{family} must be \code{"independence"}, \code{"frank"} or
#'   \code{"nonparametric"}.
#' @return An object of class \code{"proc_values"}: \code{bids} (ascending),
#'   \code{u}, \code{q}, \code{psi}, \code{markup} (\eqn{b - c}),
#'   \code{costs_raw}, \code{costs}, \code{costs_mat} (for matrix input),
#'   \code{M}, \code{n}, and \code{reflected}, the \code{fpa_values} fit on
#'   \eqn{M - b} used by the cost functionals.
#' @export
proc_values <- function(bids, n = NULL, M = NULL, monotone = "isotonic", ...) {
  .proc_check_family(list(...))
  if (is.null(M)) M <- proc_M(bids)
  if (is.matrix(bids) && is.null(n)) n <- ncol(bids)
  ## pass n only when given (see handover item 6)
  fit <- if (is.null(n)) fpa_values(M - bids, monotone = monotone, ...)
         else fpa_values(M - bids, n = n, monotone = monotone, ...)
  flip <- function(z) base::rev(as.numeric(z))
  if (length(fit$bids) != length(bids) || length(M) != 1 || !is.finite(M))
    stop("unexpected fpa_values() output", call. = FALSE)
  out <- list(
    M         = M,
    n         = n,
    bids      = M - flip(fit$bids),
    u         = 1 - flip(fit$u),
    q         = flip(fit$q),
    psi       = flip(fit$psi),
    markup    = flip(fit$markup),
    costs_raw = M - flip(fit$values_raw),
    costs     = M - flip(fit$values),
    reflected = fit)
  if (!is.null(fit$values_mat)) out$costs_mat <- M - fit$values_mat
  class(out) <- "proc_values"
  out
}

#' @export
print.proc_values <- function(x, ...) {
  rel <- (x$bids - x$costs) / x$bids
  cat("Procurement inversion by reflection (M = ", signif(x$M, 6), ")\n", sep = "")
  cat("bidders:", x$n, "  bids:", length(x$bids),
      "  negative raw costs:", sum(x$costs_raw < 0), "\n")
  cat("relative markup (b - c)/b at bid quartiles:\n")
  print(round(sapply(c(.25, .5, .75), function(p)
    mean(rel[abs(x$u - p) <= 0.025])), 4))
  invisible(x)
}

#' Expected procurement cost with a price ceiling
#'
#' \eqn{C(r) = M - R'(M - r; M - v_0)}, where \eqn{R'} is the
#' quantile-density revenue functional (\code{\link{fpa_revenue_qd}}) on the
#' reflected fit. The plug-in \code{\link{fpa_revenue}} is not used: in
#' procurement the unreliable tail of the inversion is the winners' tail.
#'
#' @param pv A \code{proc_values} object fitted to a bid matrix.
#' @param r Price ceiling(s) on the bid scale; \code{Inf} for no ceiling.
#' @param v0 Buyer's cost if no bid is accepted (outside option).
#' @param ... Passed to \code{fpa_revenue_qd}.
#' @return Numeric vector of expected costs, one per element of \code{r}.
#' @export
proc_cost <- function(pv, r = Inf, v0, ...) {
  stopifnot(inherits(pv, "proc_values"), !missing(v0))
  M <- pv$M
  M - fpa_revenue_qd(pv$reflected, r = M - r, v_s = M - v0, ...)
}

#' Optimal price ceiling for a procurement auction
#'
#' Mirrors \code{\link{fpa_optimal_reserve_qd}} on the reflected fit. An
#' interior optimum satisfies \eqn{r^* + F_c(r^*)/f_c(r^*) = v_0}.
#' The estimate depends on the lower tail of the recovered costs; check
#' \code{share_below}, the share of pseudo-costs below \eqn{r^*}, before
#' interpreting it.
#'
#' @inheritParams proc_cost
#' @return List: \code{r_star}, \code{p_star}, \code{cost}, \code{binding},
#'   \code{share_below}, and the curve \code{p}, \code{C}.
#' @export
proc_optimal_ceiling <- function(pv, v0, ...) {
  stopifnot(inherits(pv, "proc_values"))
  M <- pv$M
  o <- fpa_optimal_reserve_qd(pv$reflected, v_s = M - v0, ...)
  r_star <- M - o$r_star
  list(r_star = r_star, p_star = 1 - o$p_star, cost = M - o$revenue,
       binding = isTRUE(o$interior), share_below = mean(pv$costs < r_star),
       p = 1 - o$p, C = M - o$R)
}

## proc_cost_boot() lives in R/proc_cost_boot_blocks.R: the Project D version
## adds resample = c("auctions", "blocks") and is otherwise identical, with
## resample = "auctions" reproducing this project's function draw for draw.

#' Simulate equilibrium procurement bids
#'
#' Simulates the reflected standard auction with value quantile function
#' \eqn{Q'(u) = M - Q_c(1-u)} by \code{\link{fpa_simulate}} and maps back.
#' For costs uniform on [1, 2] this reproduces \eqn{b = c + (2-c)/n} to
#' about \eqn{10^{-6}}.
#'
#' @param L,n Auctions and bidders per auction.
#' @param Qc Cost quantile function.
#' @param qc Derivative of \code{Qc} (optional).
#' @param M Any constant above the cost support.
#' @param family,theta,eta Passed to \code{fpa_simulate}.
#' @return List: \code{costs}, \code{bids} (\eqn{L \times n}), \code{M},
#'   \code{reflected}.
#' @export
proc_simulate <- function(L, n, Qc, qc = NULL, M, family = "independence",
                          theta = NULL, eta = 0) {
  Qp <- function(u) M - Qc(1 - u)
  qp <- if (is.null(qc)) NULL else function(u) qc(1 - u)
  sim <- fpa_simulate(L, n, Q = Qp, q = qp, family = family, theta = theta, eta = eta)
  list(costs = M - sim$values, bids = M - sim$bids, M = M, reflected = sim)
}

#' Winning-bid check for the dependence specification
#'
#' Compares the model's expected winning bid without a ceiling with the
#' sample mean of the auction minima. Under the correct dependence
#' specification the two agree up to sampling error.
#'
#' @param pv A \code{proc_values} object fitted to a bid matrix.
#' @param bids The same bid matrix.
#' @return Named vector: \code{model}, \code{observed}, \code{difference}.
#' @export
proc_check_winning_bid <- function(pv, bids) {
  stopifnot(is.matrix(bids))
  m <- proc_cost(pv, Inf, v0 = max(bids))   # v0 irrelevant without a ceiling
  o <- mean(apply(bids, 1, min))
  c(model = m, observed = o, difference = m - o)
}

#' Between-auction share of residual bid variance
#'
#' One-way ANOVA intraclass correlation of (homogenised, log) bids by
#' auction. Large values indicate a component common to the bidders in an
#' auction, i.e. within-auction dependence.
#'
#' @param x Numeric vector, e.g. residuals of a log-bid regression.
#' @param auction Auction identifier, same length as \code{x}.
#' @return The ICC estimate.
#' @export
bid_icc <- function(x, auction) {
  stopifnot(length(x) == length(auction))
  a   <- stats::anova(stats::lm(x ~ factor(auction)))
  msb <- a$`Mean Sq`[1]; msw <- a$`Mean Sq`[2]
  n0  <- mean(table(auction))
  (msb - msw) / (msb + (n0 - 1) * msw)
}
