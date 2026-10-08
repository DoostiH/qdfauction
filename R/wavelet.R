# Wavelet quantile density estimators --------------------------------------

#' Wavelet estimators of the quantile density
#'
#' @description
#' Projection estimators of \eqn{q(u)} on a wavelet basis of \eqn{[0,1]}
#' (Chesneau, Dewan and Doosti 2016; Shirazi and Doosti 2022). The scaling
#' coefficients are estimated by the warped-basis formula
#' \deqn{\hat c_{j_0,k} = \int \phi_{j_0,k}(\hat F(x))\,dx
#'       = \sum_{i=1}^{n-1} (X_{(i+1)} - X_{(i)})\,\phi_{j_0,k}(i/n),}
#' which is exact for the empirical distribution function \eqn{\hat F}
#' (`integration = "exact"`); `"simpson"` reproduces the composite Simpson
#' rule of the original MATLAB implementation.
#'
#' Three estimators are available through `threshold`:
#' * `"none"`: the linear estimator \eqn{\hat g_L = \sum_k \hat c_{j_0,k}\phi_{j_0,k}};
#' * `"hard"`: term-by-term hard thresholding of a `levels`-level periodic DWT
#'   of the coefficient vector;
#' * `"block"`: Cai's BlockJS block thresholding (Shirazi and Doosti 2022).
#'
#' With `smooth = TRUE` the estimate is post-smoothed by the local-linear
#' smoother of Ramirez and Vidakovic (2010) with bandwidth `h`; this
#' "smoothed" version was the best performer in both papers' simulations.
#'
#' @inheritParams qdf_kernel
#' @param j0 coarsest resolution level; the estimate uses \eqn{2^{j_0}}
#'   scaling coefficients. Default 5.
#' @param h bandwidth of the local-linear post-smoother on the probability
#'   scale (default 0.15); ignored when `smooth = FALSE`.
#' @param threshold `"none"`, `"hard"` or `"block"`.
#' @param smooth logical; apply local-linear smoothing.
#' @param filter scaling filter for the basis; see [wavelet_filter()].
#'   Default `"db3"`.
#' @param dwt_filter filter for the secondary DWT used in thresholding
#'   (default `"coif1"`).
#' @param levels number of DWT levels for thresholding (default 3).
#' @param integration `"exact"` or `"simpson"`.
#' @param ... passed to [threshold_coefficients()] (`const`, `block_length`, `rule`).
#' @return numeric vector of quantile density estimates at `u`.
#' @references
#' Chesneau, C., Dewan, I. and Doosti, H. (2016). Nonparametric estimation of a
#' quantile density function by wavelet methods. *Computational Statistics &
#' Data Analysis*, 94, 161--174.
#'
#' Shirazi, E. and Doosti, H. (2022). Nonparametric estimation of a quantile
#' density function under \eqn{L_p} risk via block thresholding method.
#' *Communications in Statistics -- Simulation and Computation*, 51, 539--553.
#'
#' Ramirez, P. and Vidakovic, B. (2010). Wavelet density estimation for
#' stratified size-biased sample. *Journal of Statistical Planning and
#' Inference*, 140, 419--432.
#' @examples
#' set.seed(1)
#' x <- rbeta(500, 0.5, 0.5)
#' u <- seq(0.02, 0.98, by = 0.02)
#' q_true <- 1 / dbeta(qbeta(u, 0.5, 0.5), 0.5, 0.5)
#' q_lin  <- qdf_wavelet(u, x, j0 = 5, smooth = FALSE)
#' q_blk  <- qdf_wavelet(u, x, j0 = 5, threshold = "block", smooth = TRUE)
#' plot(u, q_true, type = "l"); lines(u, q_lin, col = 2); lines(u, q_blk, col = 4)
#' @export
qdf_wavelet <- function(u, x, j0 = 5L, h = 0.15,
                        threshold = c("none", "hard", "block"), smooth = TRUE,
                        filter = "db3", dwt_filter = "coif1", levels = 3L,
                        integration = c("exact", "simpson"), ...) {
  threshold <- match.arg(threshold); integration <- match.arg(integration)
  x <- prep_sample(x); u <- check_u(u); n <- length(x)
  K <- 2^j0
  if (K %% 2^levels != 0) stop("2^j0 must be divisible by 2^levels.", call. = FALSE)

  # scaling coefficients c_{j0,k}, k = 0..2^j0-1
  if (integration == "exact") {
    ug <- seq_len(n - 1) / n
    Phi <- phi_jk(ug, j0, 0:(K - 1L), filter)          # (n-1) x K
    chat <- as.numeric(diff(x) %*% Phi)
  } else {
    grid <- seq(x[1], x[n], length.out = n + 1)
    Fg <- ecdf_at(grid, x); Fg[grid <= x[1]] <- 0; Fg[grid >= x[n]] <- 1
    Phi <- phi_jk(Fg, j0, 0:(K - 1L), filter)
    w <- c(1, rep(c(4, 2), length.out = n - 1), 1)
    chat <- as.numeric((w * (x[n] - x[1]) / (3 * n)) %*% Phi)
  }

  # coefficient-vector representation used for thresholding (2^{j0/2} c_k
  # approximates g at the lattice k/2^j0, as in the MATLAB implementation)
  v <- 2^(j0 / 2) * chat
  if (threshold != "none") {
    w  <- dwt_periodic(v, levels, dwt_filter)
    w  <- threshold_coefficients(w, levels, n, type = threshold, ...)
    v  <- idwt_periodic(w, levels, dwt_filter)
  }
  lattice <- seq(0, 1, length.out = K)

  if (smooth) {
    return(loclin_smooth(u, lattice, v, h))
  }
  if (threshold == "none") {
    # proper synthesis on the basis
    Phi_u <- phi_jk(u, j0, 0:(K - 1L), filter)
    return(as.numeric(Phi_u %*% chat))
  }
  stats::approx(lattice, v, xout = u, rule = 2)$y
}

#' Cross-validated choice of (j0, h) for the wavelet estimator
#'
#' Chooses the resolution level and the post-smoothing bandwidth of
#' [qdf_wavelet()] by `folds`-fold cross-validation. The implied density of the
#' data is \eqn{\hat g(x) = 1/\hat q(\hat G_{\rm train}(x))}. With
#' `criterion = "lscv"` (default) the pair minimises the least-squares
#' criterion \eqn{\int \hat g^2 - 2 n_k^{-1} \sum \hat g(x_i)} over held-out
#' observations, using \eqn{\int \hat g(x)^2 dx = \int_0^1 1/\hat q(u)\,du};
#' this targets integrated squared error (Rudemo 1982; Bowman 1984). With
#' `criterion = "loglik"` it maximises the held-out log-likelihood. Pairs
#' with \eqn{h} below the lattice spacing \eqn{1/(2^{j_0}-1)} are not
#' considered, since the local-linear smoother is then ill-conditioned.
#'
#' @param x numeric sample.
#' @param j0_grid,h_grid candidate values.
#' @section Reproducibility:
#' Fold assignment (or bootstrap resampling) uses R's random number generator;
#' call [set.seed()] beforehand to make results exactly reproducible.
#' @param folds number of folds.
#' @param criterion `"lscv"` or `"loglik"`.
#' @param ... passed to [qdf_wavelet()] (e.g. `threshold`, `smooth`, `filter`).
#' @return a list with `j0`, `h`, and the matrix of criterion values
#'   (`ll`, rows `j0_grid`, columns `h_grid`).
#' @examples
#' set.seed(2)
#' x <- rbeta(300, 0.5, 0.5)
#' select_wavelet(x, j0_grid = 4:6, h_grid = c(0.05, 0.15), folds = 3)
#' @export
select_wavelet <- function(x, j0_grid = 3:6, h_grid = c(0.03, 0.05, 0.08, 0.12, 0.2, 0.3),
                           folds = 5L, criterion = c("lscv", "loglik"), ...) {
  criterion <- match.arg(criterion)
  x <- prep_sample(x); n <- length(x)
  ug <- seq(0.001, 0.999, length.out = 999)
  j0_grid <- j0_grid[2^j0_grid <= n]
  fold_id <- sample(rep_len(seq_len(folds), n))
  ll <- matrix(NA_real_, length(j0_grid), length(h_grid),
               dimnames = list(j0 = j0_grid, h = h_grid))
  smooth_used <- !isFALSE(list(...)$smooth)
  for (a in seq_along(j0_grid)) for (b in seq_along(h_grid)) {
    # local linear smoothing of 2^j0 lattice values is ill-conditioned when h is
    # below the lattice spacing 1/(2^j0 - 1): such pairs are not considered
    if (smooth_used && h_grid[b] < 1 / (2^j0_grid[a] - 1)) { ll[a, b] <- NA; next }
    tot <- 0
    for (f in seq_len(folds)) {
      tr <- sort(x[fold_id != f]); te <- x[fold_id == f]
      ut <- pmin(pmax(ecdf_at(te, tr), 0.5 / length(tr)), 1 - 0.5 / length(tr))
      q  <- tryCatch(qdf_wavelet(ut, tr, j0 = j0_grid[a], h = h_grid[b], ...),
                     error = function(e) NA_real_)
      q  <- pmax(q, 1e-8)
      if (criterion == "loglik") {
        tot <- tot + sum(-log(q))
      } else {                                         # -(int g^2 - 2 mean g(held out)), g = 1/q(G)
        qg <- tryCatch(pmax(qdf_wavelet(ug, tr, j0 = j0_grid[a], h = h_grid[b], ...), 1e-8),
                       error = function(e) NA_real_)
        tot <- tot - (mean(1 / qg) - 2 * mean(1 / q)) * length(te)
      }
    }
    ll[a, b] <- tot
  }
  best <- which(ll == max(ll, na.rm = TRUE), arr.ind = TRUE)[1, ]
  list(j0 = j0_grid[best[1]], h = h_grid[best[2]], ll = ll, criterion = criterion)
}
