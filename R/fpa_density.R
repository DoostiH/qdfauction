# Private-value density estimation ---------------------------------------------

.density_methods <- c("wavelet_ll", "wavelet_linear", "wavelet_hard", "wavelet_block",
                      "kde_hh", "kde", "empirical_quantile", "marmer_shneyerov")

# monotone map of bids and pseudo-values onto [0, 1] (unit_map in the MATLAB code)
unit_map <- function(x, mode = c("auto", "log", "minmax", "none"), bids = x) {
  mode <- match.arg(mode)
  # auto: identity (with pseudo-values clipped to [0,1]) when the bids live in the
  # unit interval, as in the simulation designs; otherwise the log map
  if (mode == "auto") mode <- if (all(bids >= 0 & bids <= 1)) "none" else "log"
  if (mode == "none")   # identity; wavelet methods project onto [0,1] themselves
    return(list(mode = "none", fwd = identity, inv = identity,
                jac = function(z) rep(1, length(z))))
  if (mode == "log") {
    if (any(x <= 0)) stop("`scale = \"log\"` requires positive bids and pseudo-values.", call. = FALSE)
    a <- log(min(x)); b <- log(max(x)); pad <- 0.02 * (b - a); a <- a - pad; b <- b + pad
    list(mode = "log", a = a, b = b,
         fwd = function(z) pmin(pmax((log(pmax(z, .Machine$double.xmin)) - a) / (b - a), 0), 1),
         inv = function(u) exp(a + u * (b - a)),
         jac = function(z) 1 / ((b - a) * z))                 # du/dx
  } else {
    a <- min(x); b <- max(x); pad <- 0.02 * (b - a); a <- a - pad; b <- b + pad
    list(mode = "minmax", a = a, b = b,
         fwd = function(z) pmin(pmax((z - a) / (b - a), 0), 1),
         inv = function(u) a + u * (b - a),
         jac = function(z) rep(1 / (b - a), length(z)))
  }
}

# first stage: pseudo-values with the empirical cdf and a kernel density of bids
first_stage <- function(b, n, type = c("gpv_hh", "gpv"), kernel = "triweight", side = "lower") {
  type <- match.arg(type)
  b <- sort(b); N <- length(b)
  G <- ecdf_at(b, b)
  g <- if (type == "gpv_hh") as.numeric(kde_hh(b, b, kernel, side)) else kde_at(b, b, stats::bw.SJ(b))
  g <- pmax(g, .Machine$double.eps)
  b + G / ((n - 1) * g)
}
# apply a first stage fitted on `tr` to new bids `te`
first_stage_apply <- function(tr, te, n, type = "gpv_hh", kernel = "triweight", side = "lower") {
  tr <- sort(tr); ntr <- length(tr)
  G <- pmin(pmax(ecdf_at(te, tr), 1 / (ntr + 1)), ntr / (ntr + 1))
  g <- if (type == "gpv_hh") as.numeric(kde_hh(te, tr, kernel, side)) else kde_at(te, tr, stats::bw.SJ(tr))
  te + G / ((n - 1) * pmax(g, .Machine$double.eps))
}

# wavelet lattice vector 2^{j0/2} * mean phi_{j0,k}(v) for v in [0,1]
wav_lattice <- function(v, j0, filter = "db3") {
  v <- pmin(pmax(v, 0), 1)             # projection Pi onto [0,1] (Algorithm 1)
  K <- 2^j0
  2^(j0 / 2) * colMeans(phi_jk(v, j0, 0:(K - 1L), filter))
}

# second-stage density on the unit scale, returned as a function of u
second_stage_unit <- function(v, method, j0 = 5L, h = 0.15, filter = "db3",
                              dwt_filter = "coif1", levels = 3L, kernel = "triweight",
                              side = "lower", ...) {
  N <- length(v)
  if (startsWith(method, "wavelet")) {
    lin <- wav_lattice(v, j0, filter)
    lattice <- seq(0, 1, length.out = 2^j0)
    if (method %in% c("wavelet_hard", "wavelet_block")) {
      w <- dwt_periodic(lin, levels, dwt_filter)
      w <- threshold_coefficients(w, levels, N, type = sub("wavelet_", "", method), ...)
      lin <- idwt_periodic(w, levels, dwt_filter)
    }
    if (method == "wavelet_ll") {
      pts <- seq_len(N) / (N + 1)
      sm  <- loclin_smooth(pts, lattice, lin, h)
      return(function(u) stats::approx(pts, sm, xout = u, rule = 2)$y)
    }
    return(function(u) stats::approx(lattice, lin, xout = u, rule = 2)$y)
  }
  switch(method,
         kde_hh = function(u) as.numeric(kde_hh(u, v, kernel, side)),
         kde    = { bw <- stats::bw.SJ(v); function(u) kde_at(u, v, bw) },
         empirical_quantile = {
           p  <- seq_len(N) / (N + 1)
           Qp <- stats::quantile(v, p, names = FALSE, type = 7)
           dQ <- c(diff(Qp[1:2]) / diff(p[1:2]),
                   (Qp[-(1:2)] - Qp[1:(N - 2)]) / (p[-(1:2)] - p[1:(N - 2)]),
                   diff(Qp[(N - 1):N]) / diff(p[(N - 1):N]))
           f  <- pmin(1 / pmax(dQ, 1e-10), 1e3)
           function(u) stats::approx(Qp, f, xout = u, rule = 2, ties = mean)$y   # density at v = Q(p)
         },
         stop("unknown method"))
}

#' Nonparametric estimation of the private-value density in first-price auctions
#'
#' @description
#' Two-step estimation of the density \eqn{f} of independent private values.
#' The first step recovers pseudo-values \eqn{\hat v_i} from the bids (by
#' default with the Hickman--Hubbard boundary-corrected GPV inversion, or
#' with any [fpa_values()] fit); the second step estimates the density of
#' the pseudo-values by
#' * `"wavelet_ll"`: linear wavelet estimator followed by local-linear
#'   smoothing (the best performer in Doosti, Dewan and Sbai 2026);
#' * `"wavelet_linear"`, `"wavelet_hard"`, `"wavelet_block"`: linear,
#'   hard-thresholded and block-thresholded wavelet estimators;
#' * `"kde_hh"`: the boundary-corrected kernel estimator of Hickman and
#'   Hubbard (2015), i.e. GPV with boundary correction at both stages;
#' * `"kde"`: Gaussian kernel with Sheather--Jones bandwidth;
#' * `"empirical_quantile"`: finite-difference density from the empirical
#'   quantile function of the pseudo-values;
#' * `"marmer_shneyerov"`: the quantile-based estimator of Marmer and
#'   Shneyerov (2012), computed directly from the bids (no first stage).
#'
#' Because the wavelet basis lives on \eqn{[0,1]}, bids and pseudo-values are
#' mapped to the unit interval by a monotone map (`scale`): with `"auto"` the
#' identity when the bids lie in \eqn{[0,1]} (pseudo-values outside are
#' projected onto \eqn{[0,1]} for the wavelet coefficients), otherwise \eqn{x \mapsto (\log x - a)/(b-a)}
#' (`"log"`, which handles the heavy right tail of the GPV ratio) or a linear
#' map (`"minmax"`). Densities are reported on the original scale through the
#' Jacobian.
#'
#' @param bids pooled bids (vector) or an `L` by `n` matrix.
#' @param n bidders per auction (from `ncol(bids)` for a matrix).
#' @param method second-stage estimator; see Description.
#' @param first_stage `"gpv_hh"` (default), `"gpv"`, an [fpa_values()] object,
#'   or the name of a [qdf()] method to use through [fpa_values()].
#' @param x evaluation grid on the value scale. Defaults to 200 points spanning
#'   the pseudo-values (excluding the extreme 0.5 percent at each end).
#' @param j0 resolution level, or `"cv"` to choose `j0` (and `h`) by
#'   cross-validated held-out log-likelihood ([select_fpa_density()]).
#' @param h local-linear bandwidth for `"wavelet_ll"` (ignored if `j0 = "cv"`).
#' @param scale `"auto"`, `"log"`, `"minmax"` or `"none"`.
#' @param filter,dwt_filter,levels wavelet settings; see [qdf_wavelet()].
#' @param kernel,side kernel and correction side for [kde_hh()] (used in the
#'   `"gpv_hh"` first stage and the `"kde_hh"` second stage).
#' @param ... for `first_stage` given as a [qdf()] method name, further
#'   arguments to [fpa_values()]; for `"wavelet_block"`, arguments to
#'   [threshold_coefficients()]; for `j0 = "cv"`, `j0_grid`, `h_grid`, `folds`,
#'   `criterion`.
#' @return An object of class `"fpa_density"`: a list with `x`, `f` (density
#'   on the value scale), `values` (pseudo-values), `values_u` (on the unit
#'   scale), `map`, `method`, `j0`, `h`, `first_stage`, `n`, and `fun`, a
#'   function returning the density at arbitrary values.
#' @references
#' Doosti, H., Dewan, I. and Sbai, E. (2026). Adaptive wavelet estimation of
#' private-value densities in first-price auctions. Working paper.
#'
#' Hickman, B. R. and Hubbard, T. P. (2015). *Journal of Applied
#' Econometrics*, 30, 739--762.
#'
#' Marmer, V. and Shneyerov, A. (2012). *Journal of Econometrics*, 167,
#' 345--357.
#' @examples
#' set.seed(1)
#' sim <- fpa_simulate(40, 5, Q = function(u) qbeta(u, 2, 2))
#' fd <- fpa_density(sim$bids, method = "wavelet_ll", j0 = 5, h = 0.1)
#' plot(fd, true = function(v) dbeta(v, 2, 2))
#' # timber auctions, three bidders, log scale
#' data(timber)
#' fd3 <- fpa_density(timber$n3, method = "wavelet_ll", j0 = 4, h = 0.15)
#' plot(fd3)
#' @export
fpa_density <- function(bids, n = NULL, method = c("wavelet_ll", "wavelet_linear", "wavelet_hard",
                                                 "wavelet_block", "kde_hh", "kde",
                                                 "empirical_quantile", "marmer_shneyerov"),
                        first_stage = "gpv_hh", x = NULL, j0 = 5L, h = 0.15,
                        scale = c("auto", "log", "minmax", "none"),
                        filter = "db3", dwt_filter = "coif1", levels = 3L,
                        kernel = "triweight", side = c("lower", "both"), ...) {
  method <- match.arg(method); scale <- match.arg(scale); side <- match.arg(side)
  dots <- list(...)
  if (is.matrix(bids)) n <- ncol(bids)
  if (is.null(n) || n < 2) stop("`n` (bidders per auction, >= 2) is required.", call. = FALSE)
  b <- sort(as.numeric(bids)); N <- length(b)

  # ---- first stage
  fs_label <- if (is.character(first_stage)) first_stage else "fpa_values"
  if (method == "marmer_shneyerov") {
    v <- b; fs_label <- "none"
  } else if (inherits(first_stage, "fpa_values")) {
    v <- sort(first_stage$values)
  } else if (first_stage %in% c("gpv_hh", "gpv")) {
    v <- first_stage(b, n, first_stage, kernel, side)
  } else {
    fv_args <- dots[intersect(names(dots), names(formals(fpa_values)))]
    dots <- dots[setdiff(names(dots), names(fv_args))]
    v <- sort(do.call(fpa_values, c(list(bids = b, n = n, method = first_stage), fv_args))$values)
  }

  # ---- evaluation grid
  if (is.null(x)) {
    r <- stats::quantile(v, c(0.005, 0.995), names = FALSE)
    x <- seq(r[1], r[2], length.out = 200)
  }

  if (method == "marmer_shneyerov") {
    f <- density_marmer_shneyerov(x, b, n)
    fun <- function(xx) density_marmer_shneyerov(xx, b, n)
    map <- unit_map(b, "none"); vu <- NULL
  } else {
    map <- unit_map(c(b, v), scale, bids = b)
    vu  <- map$fwd(v)
    if (identical(j0, "cv")) {
      cv_args <- dots[intersect(names(dots), c("j0_grid", "h_grid", "folds", "criterion", "trim"))]
      sel <- do.call(select_fpa_density, c(list(bids = b, n = n, method = method,
                                               first_stage = if (is.character(first_stage)) first_stage else "gpv_hh",
                                               scale = scale, filter = filter, dwt_filter = dwt_filter,
                                               levels = levels, kernel = kernel, side = side), cv_args))
      j0 <- sel$j0; h <- sel$h
    } else sel <- NULL
    thr_args <- dots[intersect(names(dots), c("const", "block_length", "rule"))]
    fu <- do.call(second_stage_unit, c(list(v = vu, method = method, j0 = j0, h = h, filter = filter,
                                            dwt_filter = dwt_filter, levels = levels, kernel = kernel,
                                            side = side), thr_args))
    fun <- function(xx) pmax(fu(map$fwd(xx)) * map$jac(xx), 0)
    f <- fun(x)
  }
  structure(list(x = x, f = f, fun = fun, values = v, values_u = vu, bids = b, n = n,
                 method = method, first_stage = fs_label, map = map,
                 j0 = if (startsWith(method, "wavelet")) j0 else NULL,
                 h = if (method == "wavelet_ll") h else NULL,
                 selection = if (exists("sel", inherits = FALSE)) sel else NULL,
                 call = match.call()),
            class = "fpa_density")
}

#' Cross-validated choice of the wavelet resolution level and bandwidth
#'
#' @description
#' Chooses `j0` (and `h` for `"wavelet_ll"`) by `folds`-fold cross-validation
#' on the bids: the first stage and the unit map are fitted on the training
#' bids and applied to the held-out bids, and the second-stage density
#' \eqn{\hat f_{-k}} is scored on the held-out pseudo-values by one of two
#' criteria:
#' * `"lscv"` (default): least-squares cross-validation (Rudemo 1982; Bowman
#'   1984), \eqn{\sum_k \{ \int \hat f_{-k}^2 - 2\,\bar{\hat f}_{-k}(\text{held-out}) \}},
#'   an unbiased estimate of the integrated squared error up to a constant, so
#'   that the selected pair targets MISE;
#' * `"loglik"`: the held-out log-likelihood \eqn{\sum_k \sum_i \log \hat f_{-k}(\hat v_i)},
#'   which targets Kullback--Leibler divergence. Because the held-out points are
#'   pseudo-values carrying first-stage error, this criterion tends to favour
#'   the smallest bandwidth on the grid; it is kept for comparison.
#'
#' Both criteria are evaluated on the interior \eqn{[\delta, 1-\delta]} of the
#' unit interval (`trim`): pseudo-values that fall outside \eqn{[0,1]} are
#' projected onto the boundary (Algorithm 1 of the paper), which creates a
#' spurious point mass there that a small bandwidth reproduces; excluding a
#' boundary strip from the score removes this artefact, in the spirit of the
#' boundary trimming of Guerre, Perrigne and Vuong (2000). For `"lscv"` the
#' trimmed criterion is \eqn{\int_\delta^{1-\delta}\hat f^2 - 2 n_k^{-1}\sum_{\hat v_i \in [\delta,1-\delta]} \hat f(\hat v_i)},
#' an unbiased estimate of the ISE over the interior up to a constant.
#'
#' The returned score is oriented so that larger is better in both cases
#' (`-LSCV` for `"lscv"`).
#'
#' @inheritParams fpa_density
#' @param j0_grid,h_grid candidate values.
#' @param folds number of folds.
#' @param criterion `"lscv"` or `"loglik"`.
#' @param n_int number of points for the numerical integral of \eqn{\hat f^2}.
#' @param trim boundary strip excluded from the score on the unit scale
#'   (default 0.05; 0 scores the whole interval).
#' @return a list with `j0`, `h`, `criterion` and the matrix of summed scores
#'   `ll` (rows `j0_grid`, columns `h_grid`; larger is better).
#' @references
#' Rudemo, M. (1982). Empirical choice of histograms and kernel density
#' estimators. *Scandinavian Journal of Statistics*, 9, 65--78.
#'
#' Bowman, A. W. (1984). An alternative method of cross-validation for the
#' smoothing of density estimates. *Biometrika*, 71, 353--360.
#' @export
select_fpa_density <- function(bids, n, method = "wavelet_ll", first_stage = "gpv_hh",
                               j0_grid = 3:6, h_grid = c(0.03, 0.05, 0.08, 0.12, 0.2, 0.3, 0.5),
                               folds = 5L, criterion = c("lscv", "loglik"), scale = "auto",
                               filter = "db3", dwt_filter = "coif1", levels = 3L,
                               kernel = "triweight", side = "lower", n_int = 2001L,
                               trim = 0.05) {
  criterion <- match.arg(criterion)
  ug <- seq(trim, 1 - trim, length.out = n_int)
  b <- sort(as.numeric(bids)); N <- length(b)
  if (method != "wavelet_ll") h_grid <- h_grid[1]
  j0_grid <- j0_grid[2^j0_grid <= N & 2^j0_grid %% 2^levels == 0]
  fold <- sample(rep_len(seq_len(folds), N))
  ll <- matrix(0, length(j0_grid), length(h_grid), dimnames = list(j0 = j0_grid, h = h_grid))
  for (k in seq_len(folds)) {
    tr <- b[fold != k]; te <- b[fold == k]
    if (first_stage %in% c("gpv_hh", "gpv")) {
      vtr <- first_stage(tr, n, first_stage, kernel, side)
      vte <- first_stage_apply(tr, te, n, first_stage, kernel, side)
    } else {
      ftr <- fpa_values(tr, n = n, method = first_stage)
      vtr <- sort(ftr$values)
      ute <- pmin(pmax(ecdf_at(te, tr), 1 / (length(tr) + 1)), length(tr) / (length(tr) + 1))
      vte <- te + stats::predict(ftr$qdf, ute) * ute / (n - 1)
    }
    mp <- unit_map(c(tr, vtr, vte), scale, bids = tr)
    utr <- mp$fwd(vtr); ute <- mp$fwd(vte)
    n_te <- length(ute)
    ute  <- ute[ute >= trim & ute <= 1 - trim]          # score on the interior only
    for (a in seq_along(j0_grid)) for (bb in seq_along(h_grid)) {
      fu <- tryCatch(second_stage_unit(utr, method, j0_grid[a], h_grid[bb], filter, dwt_filter,
                                       levels, kernel, side), error = function(e) NULL)
      sc <- if (is.null(fu)) -Inf else if (criterion == "loglik") {
        sum(log(pmax(fu(ute), 1e-8)))
      } else {
        f2 <- fu(ug)^2
        -(sum(diff(ug) * (f2[-1] + f2[-n_int]) / 2) - 2 * sum(fu(ute)) / n_te) * n_te
      }
      ll[a, bb] <- ll[a, bb] + sc
    }
  }
  best <- which(ll == max(ll), arr.ind = TRUE)[1, ]
  list(j0 = j0_grid[best[1]], h = if (method == "wavelet_ll") h_grid[best[2]] else NA_real_,
       ll = ll, criterion = criterion)
}

#' Compare density estimators by held-out log-likelihood
#'
#' Outer `folds`-fold cross-validation of several second-stage estimators
#' against a common target (the held-out pseudo-values from the first stage
#' fitted on the training fold), as in Table 3 of Doosti, Dewan and Sbai
#' (2026). Wavelet tuning parameters are fixed at `j0`, `h`.
#'
#' @inheritParams fpa_density
#' @param methods character vector of second-stage methods.
#' @param folds number of outer folds.
#' @return a named numeric vector of summed held-out log-likelihoods.
#' @export
fpa_density_cv <- function(bids, n = NULL, methods = c("wavelet_ll", "kde_hh", "marmer_shneyerov"),
                           folds = 5L, j0 = 5L, h = 0.15, first_stage = "gpv_hh",
                           scale = "auto", kernel = "triweight", side = "lower", ...) {
  if (is.matrix(bids)) n <- ncol(bids)
  b <- sort(as.numeric(bids)); N <- length(b)
  fold <- sample(rep_len(seq_len(folds), N))
  out <- stats::setNames(numeric(length(methods)), methods)
  for (k in seq_len(folds)) {
    tr <- b[fold != k]; te <- b[fold == k]
    vte <- first_stage_apply(tr, te, n, first_stage, kernel, side)
    for (m in methods) {
      fd <- fpa_density(tr, n, method = m, first_stage = first_stage, j0 = j0, h = h,
                        scale = scale, kernel = kernel, side = side, ...)
      out[m] <- out[m] + sum(log(pmax(fd$fun(vte), 1e-8)))
    }
  }
  out
}

#' @export
print.fpa_density <- function(x, digits = 4, ...) {
  cat("Private-value density, first-price auction\n")
  cat("  bids        :", length(x$bids), "pooled,", x$n, "bidders per auction\n")
  cat("  first stage :", x$first_stage, "\n")
  cat("  method      :", x$method,
      if (!is.null(x$j0)) sprintf("(j0 = %d%s)", x$j0, if (!is.null(x$h)) sprintf(", h = %s", format(x$h, digits = digits)) else ""), "\n")
  cat("  scale map   :", x$map$mode, "\n")
  if (!is.null(x$values_u))
    cat("  pseudo-values:", format(range(x$values), digits = digits), "; max v / max b =",
        format(max(x$values) / max(x$bids), digits = 3), "\n")
  invisible(x)
}

#' Plot an estimated private-value density
#' @param x an `"fpa_density"` object.
#' @param true optional true density function for simulated data.
#' @param rug logical; add a rug of the pseudo-values.
#' @param ... passed to [graphics::plot()].
#' @export
plot.fpa_density <- function(x, true = NULL, rug = TRUE, ...) {
  args <- list(x = x$x, y = x$f, type = "l", lwd = 2, xlab = "private value",
               ylab = "density", main = paste("Private-value density,", x$method))
  args <- utils::modifyList(args, list(...))
  do.call(graphics::plot, args)
  if (rug && !is.null(x$values)) graphics::rug(x$values, col = "grey60")
  if (!is.null(true)) {
    graphics::lines(x$x, true(x$x), lty = 2, col = 2, lwd = 2)
    graphics::legend("topright", c("estimate", "true"), lty = 1:2, col = 1:2, bty = "n")
  }
  invisible(x)
}
