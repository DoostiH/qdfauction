# Evaluation and Monte Carlo ---------------------------------------------------

#' Evaluate recovered values or an estimated density against the truth
#'
#' @description
#' Performance criteria used in the papers behind this package.
#'
#' `fpa_evaluate_values()` scores recovered private values (Doosti, Dewan and
#' Talebian 2025): `mise`, the mean integrated squared error between a
#' Gaussian kernel density (Sheather--Jones bandwidth) of the values and the
#' true density; `loglik`, \eqn{\sum_i \log f(\hat v_i)}; `ks`, the
#' Kolmogorov--Smirnov distance between the empirical distribution of the
#' values and \eqn{F}; and `msep`, the mean squared error against the true
#' values when they are supplied (Doosti 2026a), with `msep_trim` excluding a
#' share `trim` at each end of the true-value distribution.
#'
#' `fpa_evaluate_density()` scores an [fpa_density()] object on its grid:
#' `mise` over the whole grid and over the upper `tail` share of the support,
#' and the expected log-likelihood \eqn{\int f \log \hat f} (Doosti, Dewan and
#' Sbai 2026).
#'
#' @param values numeric vector of recovered values, or an [fpa_values()] object.
#' @param pdf,cdf true density and distribution functions.
#' @param true optional true values, aligned with `values` (for an
#'   [fpa_values()] object, in the layout of the original bids).
#' @param trim trimming share for `msep_trim`.
#' @return a named numeric vector.
#' @examples
#' set.seed(1)
#' sim <- fpa_simulate(20, 5, Q = function(u) qgamma(u, 10, 2))
#' fv  <- fpa_values(sim$bids, method = "bernstein", bandwidth = "bcv", monotone = "spline")
#' fpa_evaluate_values(fv, pdf = function(x) dgamma(x, 10, 2),
#'                     cdf = function(x) pgamma(x, 10, 2), true = sim$values)
#' @export
fpa_evaluate_values <- function(values, pdf, cdf = NULL, true = NULL, trim = 0.05) {
  if (inherits(values, "fpa_values")) {
    v <- values$values
    if (!is.null(true)) true <- as.numeric(true)[values$order]
  } else v <- as.numeric(values)
  out <- c()
  d <- stats::density(v, bw = "SJ")
  out["mise"] <- mean((pdf(d$x) - d$y)^2) * diff(range(d$x))
  out["loglik"] <- sum(log(pdf(v)))
  if (!is.null(cdf)) out["ks"] <- max(abs(stats::ecdf(v)(sort(v)) - cdf(sort(v))),
                                     abs((seq_along(v) - 1) / length(v) - cdf(sort(v))))
  if (!is.null(true)) {
    out["msep"] <- mean((v - true)^2)
    qs <- stats::quantile(true, c(trim, 1 - trim))
    keep <- true >= qs[1] & true <= qs[2]
    out["msep_trim"] <- mean((v[keep] - true[keep])^2)
  }
  out
}

#' @rdname fpa_evaluate_values
#' @param density an [fpa_density()] object.
#' @param tail share of the upper support for the tail MISE.
#' @export
fpa_evaluate_density <- function(density, pdf, tail = 0.1) {
  x <- density$x; f <- density$f; ft <- pdf(x)
  k <- floor((1 - tail) * length(x)):length(x)
  c(mise = mean((f - ft)^2), mise_tail = mean((f[k] - ft[k])^2),
    eloglik = sum(diff(x) * (ft * log(pmax(f, 1e-12)))[-1]))
}

#' Monte Carlo comparison of estimators
#'
#' Repeats: simulate auctions with [fpa_simulate()], apply each estimator in
#' `estimators` to the bids, and score it with `score`. Designed to
#' regenerate the simulation tables of the papers behind this package.
#'
#' @param R number of replications.
#' @param L,n auctions and bidders per replication.
#' @param Q,q,family,theta private-value model, as in [fpa_simulate()].
#' @param estimators a named list of functions `function(bids, values)`
#'   returning an object accepted by `score` (typically [fpa_values()] or
#'   [fpa_density()] fits, or a numeric vector).
#' @param score a function `function(fit, bids, values)` returning a named
#'   numeric vector, e.g. a wrapper around [fpa_evaluate_values()].
#' @param seed optional seed.
#' @param verbose print progress every 10 replications.
#' @return a data frame with columns `rep`, `estimator` and one column per
#'   criterion; `summary()` on the result (via `aggregate`) gives the table.
#' @examples
#' \donttest{
#' set.seed(1)
#' pdf <- function(x) dgamma(x, 10, 2); cdf <- function(x) pgamma(x, 10, 2)
#' mc <- fpa_montecarlo(R = 5, L = 20, n = 5, Q = function(u) qgamma(u, 10, 2),
#'   estimators = list(
#'     B_B  = function(b, v) fpa_values(b, method = "bernstein", bandwidth = "bcv"),
#'     B_IP = function(b, v) fpa_values(b, method = "indirect_poisson", bandwidth = "bcv"),
#'     GPV  = function(b, v) fpa_values(b, method = "gpv")),
#'   score = function(fit, b, v) fpa_evaluate_values(fit, pdf, cdf, v))
#' aggregate(cbind(mise, loglik, ks, msep) ~ estimator, mc, mean)
#' }
#' @export
fpa_montecarlo <- function(R, L, n, Q, q = NULL, family = "independence", theta = NULL,
                           estimators, score, seed = NULL, verbose = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  bf <- fpa_bid_function(Q, q, n, family, theta)
  rows <- vector("list", R * length(estimators)); k <- 0L
  for (r in seq_len(R)) {
    sim <- fpa_simulate(L, n, Q, q, family, theta, bid_function = bf)
    for (nm in names(estimators)) {
      sc <- tryCatch(score(estimators[[nm]](sim$bids, sim$values), sim$bids, sim$values),
                     error = function(e) NULL)
      k <- k + 1L
      rows[[k]] <- if (is.null(sc)) data.frame(rep = r, estimator = nm, failed = TRUE)
                   else data.frame(rep = r, estimator = nm, failed = FALSE, as.list(sc))
    }
    if (verbose && r %% 10 == 0) message("replication ", r, "/", R)
  }
  out <- do.call(rbind, lapply(rows, function(d) { d[setdiff(unique(unlist(lapply(rows, names))), names(d))] <- NA; d }))
  out[, unique(unlist(lapply(rows, names)))]
}
