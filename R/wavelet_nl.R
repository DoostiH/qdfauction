# Nonlinear wavelet estimators with theoretical thresholds -----------------------

# keep/kill wavelet coefficients level by level
.wavelet_threshold <- function(d, levels, rule, n, lambda_fun, block_len, p_block, c_block) {
  s <- stats::median(abs(d[[length(d)]])) / 0.6745          # noise scale, finest level
  for (a in seq_along(levels)) {
    j <- levels[a]; dj <- d[[a]]
    if (rule == "block") {
      l <- max(1L, min(block_len, length(dj)))
      thr <- if (is.null(c_block)) 4.505^(p_block / 2) * s^p_block else c_block * n^(-p_block / 2)
      for (b in split(seq_along(dj), ceiling(seq_along(dj) / l)))
        if (mean(abs(dj[b])^p_block) <= thr) dj[b] <- 0
    } else {
      dj[abs(dj) < lambda_fun(j, s)] <- 0
    }
    d[[a]] <- dj
  }
  d
}

# synthesise sum_k c phi_{tau,k} + sum_j sum_k d_jk psi_{j,k} on a grid
.wavelet_synth <- function(u, cc, d, tau, levels, filter) {
  out <- as.numeric(phi_per(u, tau, 0:(2^tau - 1L), filter) %*% cc)
  for (a in seq_along(levels)) if (any(d[[a]] != 0)) {
    j <- levels[a]; k <- which(d[[a]] != 0) - 1L
    out <- out + as.numeric(as.matrix(psi_per(u, j, k, filter)) %*% d[[a]][k + 1L])
  }
  out
}

#' Nonlinear wavelet estimators of the quantile density with theoretical thresholds
#'
#' @description
#' Hard-thresholding and block-thresholding estimators of the quantile density
#' built on the wavelet coefficients of the quantile density itself,
#' \deqn{\hat c_{\tau,k} = \int \phi_{\tau,k}(\hat F(x))\,dx, \qquad
#'       \hat d_{j,k} = \int \psi_{j,k}(\hat F(x))\,dx,}
#' computed exactly for the empirical distribution function, on a periodised
#' orthonormal basis of \eqn{[0,1]}:
#' \deqn{\hat q(u) = \sum_k \hat c_{\tau,k}\phi_{\tau,k}(u) +
#'   \sum_{j=\tau}^{j_1} \sum_k \hat d_{j,k}\, \mathrm{keep}_{j,k}\, \psi_{j,k}(u).}
#' * `rule = "hard_general"`: hard thresholding with
#'   \eqn{\lambda_j = K 2^{3j/2}\sqrt{2 p \ln n / n}}, \eqn{K = \sup|\psi'|} and
#'   \eqn{2^{j_1} \approx (n/\ln n)^{1/4}} (Chesneau, Dewan and Doosti 2016,
#'   Theorem 3.2, no smoothness condition beyond the Besov ball);
#' * `rule = "hard_lipschitz"`: hard thresholding with
#'   \eqn{\lambda_j = \kappa\sqrt{\ln n / n}} and \eqn{2^{j_1} \approx \sqrt{n/\ln n}}
#'   (Theorem 3.4, under a Lipschitz-1/2 condition on \eqn{q});
#' * `rule = "block"`: block thresholding (Shirazi and Doosti 2022, eq. 6):
#'   at each level the coefficients are grouped in blocks of length
#'   \eqn{l = \lfloor\ln n\rfloor} and a block is kept when
#'   \eqn{l^{-1}\sum_{j \in B}|\hat d_{j,k}|^p > c\, n^{-p/2}}, with
#'   \eqn{j_1 = \lfloor\log_2(n / l^2)\rfloor}.
#'
#' The constants of the theorems are not available in practice. By default
#' \eqn{\kappa = \hat s\sqrt{2n}}, so that \eqn{\lambda = \hat s\sqrt{2\ln n}}
#' (the universal threshold), and \eqn{c = 4.505^{p/2}\hat s^p n^{p/2}}, so that a block
#' is kept when its mean \eqn{|\hat d|^p} exceeds \eqn{4.505^{p/2}\hat s^p}
#' (the BlockJS level of Cai 1999), where \eqn{\hat s} is the median absolute
#' coefficient at the finest level divided by 0.6745. With `smooth = TRUE` the
#' estimate is post-smoothed by local linear regression (Ramirez and Vidakovic
#' 2010), the "smoothed" versions of both papers.
#'
#' The periodised basis replaces the boundary-corrected basis of Cohen et al.
#' (1993) assumed by the theory; periodisation identifies the endpoints of
#' \eqn{[0,1]}, so near \eqn{u = 0} and \eqn{u = 1} these estimators should be
#' read with care, and the smoothed versions are usually preferable there.
#'
#' @inheritParams qdf_kernel
#' @param rule `"hard_general"`, `"hard_lipschitz"` or `"block"`.
#' @param tau coarsest level (default 3; must satisfy \eqn{2^\tau \ge} filter length).
#' @param j1 finest level; default from the rule as above.
#' @param const \eqn{\kappa} for `"hard_lipschitz"`, or \eqn{c} for `"block"`;
#'   `NULL` for the data-driven defaults above.
#' @param p risk exponent: enters \eqn{\lambda_j} for `"hard_general"` and the
#'   block statistic for `"block"` (default 2).
#' @param smooth logical; local linear post-smoothing.
#' @param h bandwidth of the post-smoother.
#' @param filter scaling filter; see [wavelet_filter()].
#' @return numeric vector of estimates at `u`, with attribute `"kept"`, the
#'   number of detail coefficients retained at each level.
#' @references
#' Chesneau, C., Dewan, I. and Doosti, H. (2016). *Computational Statistics &
#' Data Analysis*, 94, 161--174.
#'
#' Shirazi, E. and Doosti, H. (2022). *Communications in Statistics --
#' Simulation and Computation*, 51, 539--553.
#' @examples
#' set.seed(1)
#' x <- rbeta(500, 0.5, 0.5)
#' u <- seq(0.05, 0.95, by = 0.05)
#' q_true <- 1 / dbeta(qbeta(u, 0.5, 0.5), 0.5, 0.5)
#' cbind(q_true,
#'       lipschitz = qdf_wavelet_nl(u, x, "hard_lipschitz", smooth = TRUE),
#'       block     = qdf_wavelet_nl(u, x, "block", smooth = TRUE))
#' @export
qdf_wavelet_nl <- function(u, x, rule = c("hard_lipschitz", "hard_general", "block"),
                           tau = 3L, j1 = NULL, const = NULL, p = 2, smooth = FALSE,
                           h = 0.15, filter = "db3") {
  rule <- match.arg(rule)
  x <- prep_sample(x); u <- check_u(u); n <- length(x)
  if (2^tau < length(wavelet_filter(filter))) stop("`tau` too small for this filter.", call. = FALSE)
  lnn <- log(n); l <- max(1L, floor(lnn))
  if (is.null(j1)) j1 <- switch(rule,
    hard_general   = floor(log2((n / lnn)^(1 / 4))),
    hard_lipschitz = floor(log2(sqrt(n / lnn))),
    block          = floor(log2(n / l^2)))
  j1 <- max(j1, tau)
  levels <- tau:j1
  zi <- seq_len(n - 1) / n; dx <- diff(x)
  cc <- as.numeric(dx %*% phi_per(zi, tau, 0:(2^tau - 1L), filter))
  d  <- lapply(levels, function(j) as.numeric(dx %*% psi_per(zi, j, 0:(2^j - 1L), filter)))
  K  <- .psi_deriv_sup(filter)
  lam <- switch(rule,
    hard_general   = function(j, s) K * 2^(3 * j / 2) * sqrt(2 * p * lnn / n),
    hard_lipschitz = function(j, s) (if (is.null(const)) s * sqrt(2 * n) else const) * sqrt(lnn / n),
    block          = NULL)
  d <- .wavelet_threshold(d, levels, rule, n, lam, l, p, if (rule == "block") const else NULL)
  kept <- vapply(d, function(z) sum(z != 0), integer(1)); names(kept) <- levels
  if (smooth) {
    grid <- seq(0, 1, length.out = 2^(j1 + 1))
    est  <- .wavelet_synth(grid, cc, d, tau, levels, filter)
    out  <- loclin_smooth(u, grid, est, h)
  } else out <- .wavelet_synth(u, cc, d, tau, levels, filter)
  attr(out, "kept") <- kept
  out
}
