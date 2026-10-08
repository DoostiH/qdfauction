# Wavelet machinery -----------------------------------------------------------
#
# R translation of the MATLAB routines of Brani Vidakovic (Phijk.m, dwtr.m,
# Idwtr.m, loc_lin.m; (C) B. Vidakovic, see inst/COPYRIGHTS) used in Ramirez and
# Vidakovic (2010), Chesneau, Dewan and Doosti (2016) and Shirazi and Doosti
# (2022): Daubechies-Lagarias evaluation
# of the scaling function, the periodic discrete wavelet transform, hard and
# block (Cai 1999) thresholding, and the local-linear post-smoother of
# Ramirez and Vidakovic (2010). No external wavelet package is needed.

# orthonormal scaling filters (low-pass, sum = sqrt(2))
.wavelet_filters <- list(
  haar  = c(1, 1) / sqrt(2),
  db2   = c(0.482962913144534, 0.836516303737808, 0.224143868042013, -0.129409522551260),
  db3   = c(0.332670552950083, 0.806891509311093, 0.459877502118492, -0.135011020010255,
            -0.085441273882027, 0.035226291885710),
  db4   = c(0.230377813308897, 0.714846570552916, 0.630880767929859, -0.027983769416860,
            -0.187034811719093, 0.030841381835561, 0.032883011666885, -0.010597401785069),
  sym4  = c(-0.075765714789341, -0.029635527645954, 0.497618667632458, 0.803738751805216,
            0.297857795605542, -0.099219543576935, -0.012603967262261, 0.032223100604071),
  coif1 = c(-0.015655728135465, -0.072732619512854, 0.384864846864203, 0.852572020212255,
            0.337897662457809, -0.072732619512854),
  coif2 = c(0.016387336463522, -0.041464936781979, -0.067372554721963, 0.386110066823092,
            0.812723635445542, 0.417005184421693, -0.076488599078306, -0.059434418646457,
            0.023680171946334, 0.005611434819394, -0.001823208870703, -0.000720549445364)
)

#' Orthonormal wavelet filters
#'
#' Returns the scaling (low-pass) filter of a compactly supported orthonormal
#' wavelet. Available: `"haar"`, `"db2"`, `"db3"`, `"db4"` (Daubechies with 2--4
#' vanishing moments), `"sym4"` (Symmlet 8-tap), `"coif1"`, `"coif2"`
#' (Coiflets). A numeric vector is returned unchanged after checking that it
#' sums to \eqn{\sqrt 2}.
#' @param filter character name or numeric filter.
#' @return numeric filter coefficients.
#' @examples
#' wavelet_filter("db3")
#' @export
wavelet_filter <- function(filter = "db3") {
  if (is.numeric(filter)) {
    if (abs(sum(filter) - sqrt(2)) > 1e-6) stop("A scaling filter must sum to sqrt(2).", call. = FALSE)
    return(filter)
  }
  f <- .wavelet_filters[[match.arg(filter, names(.wavelet_filters))]]
  f
}

# Daubechies-Lagarias matrices
.dl_matrices <- function(f) {
  n <- length(f); nn <- n - 1L
  T0 <- matrix(0, nn, nn); T1 <- matrix(0, nn, nn)
  for (i in seq_len(nn)) for (j in seq_len(nn)) {
    a <- 2 * i - j; b <- a + 1
    if (a > 0 && a <= n) T0[i, j] <- sqrt(2) * f[a]
    if (b > 0 && b <= n) T1[i, j] <- sqrt(2) * f[b]
  }
  list(T0 = T0, T1 = T1)
}

# cache of tabulated scaling functions
.phi_cache <- new.env(parent = emptyenv())

#' Tabulate the scaling function by the Daubechies--Lagarias algorithm
#'
#' Evaluates \eqn{\phi} on the dyadic grid \eqn{k/2^{L}} of its support
#' \eqn{[0, N-1]} (\eqn{N} = filter length) by the Daubechies--Lagarias
#' product formula (Vidakovic 1999, Section 3.4), and caches the result.
#' [phi_jk()] interpolates from this table. Translated from the MATLAB routine
#' `Phijk` of B. Vidakovic (Ramirez and Vidakovic 2010).
#' @param filter see [wavelet_filter()].
#' @param levels dyadic resolution of the table (default 10, i.e. 1024 points
#'   per unit interval).
#' @param extra number of additional (zero) binary digits used to converge the
#'   Daubechies--Lagarias product (default 20).
#' @return a list with `x` and `phi`.
#' @keywords internal
#' @export
scaling_function <- function(filter = "db3", levels = 10L, extra = 20L) {
  f <- wavelet_filter(filter)
  key <- paste(c(format(f, digits = 15), levels, extra), collapse = "_")
  if (!is.null(.phi_cache[[key]])) return(.phi_cache[[key]])
  N <- length(f) - 1L
  M <- .dl_matrices(f)
  # phi at the integers is the eigenvector for eigenvalue 1 of T0 (interior), 0 at ends
  # we build phi at dyadic points recursively: phi(x) for x = int + dec, dec = 0.b1 b2 ... bL
  dec <- (0:(2^levels - 1L)) / 2^levels
  vals <- matrix(0, N, length(dec))
  # product over the binary digits of dec, computed by a recursion on the digits
  # P(dec) = T_{b1} T_{b2} ... T_{bL}; mean of rows gives phi(int + dec) for int = 0..N-1
  prods <- vector("list", 2^levels)
  prods[[1]] <- diag(N)
  # iterate levels: at step l, prods for prefixes of length l
  cur <- list(diag(N))
  for (l in seq_len(levels)) {
    nxt <- vector("list", 2 * length(cur))
    for (i in seq_along(cur)) {
      nxt[[2 * i - 1]] <- cur[[i]] %*% M$T0
      nxt[[2 * i]]     <- cur[[i]] %*% M$T1
    }
    cur <- nxt
  }
  # trailing binary zeros (as in the MATLAB implementation with 25 digits):
  # multiply by T0^extra so that the product has converged to its rank-one limit
  T0p <- diag(N); for (i in seq_len(extra)) T0p <- T0p %*% M$T0
  for (i in seq_along(cur)) vals[, i] <- rowMeans(cur[[i]] %*% T0p)
  x   <- as.vector(outer(dec, 0:(N - 1L), "+"))
  phi <- as.vector(t(vals))
  o <- order(x); x <- x[o]; phi <- phi[o]
  out <- list(x = c(0, x[x > 0], N), phi = c(0, phi[x > 0], 0), N = N)
  # normalise so that the integral is exactly 1 (removes tabulation error)
  int <- sum(diff(out$x) * (out$phi[-1] + out$phi[-length(out$phi)]) / 2)
  out$phi <- out$phi / int
  .phi_cache[[key]] <- out
  out
}

#' Evaluate scaling functions \eqn{\phi_{j,k}(z) = 2^{j/2}\phi(2^j z - k)}
#'
#' @param z numeric vector of arguments.
#' @param j resolution level.
#' @param k integer vector of shifts. For a single `k` a vector is returned;
#'   otherwise a matrix `length(z)` by `length(k)`.
#' @inheritParams scaling_function
#' @return numeric vector or matrix of values.
#' @examples
#' z <- seq(0, 1, by = 0.01)
#' plot(z, phi_jk(z, 3, 2), type = "l")
#' @export
phi_jk <- function(z, j, k, filter = "db3", levels = 10L) {
  tab <- scaling_function(filter, levels)
  ev <- function(x) {
    out <- numeric(length(x))
    inside <- x > 0 & x < tab$N
    if (any(inside)) out[inside] <- stats::approx(tab$x, tab$phi, xout = x[inside])$y
    out
  }
  if (length(k) == 1L) return(2^(j / 2) * ev(2^j * z - k))
  2^(j / 2) * vapply(k, function(kk) ev(2^j * z - kk), numeric(length(z)))
}

#' Periodic discrete wavelet transform
#'
#' Forward (`dwt_periodic()`) and inverse (`idwt_periodic()`) periodised
#' orthonormal DWT with `L` levels, as in Vidakovic (1999). `length(x)` must be
#' divisible by `2^L`. Coefficients are returned as
#' `c(smooth, detail_coarsest, ..., detail_finest)`. Translated from the
#' MATLAB routines `dwtr` and `Idwtr` of B. Vidakovic.
#' @param x numeric vector.
#' @param L number of levels.
#' @param filter see [wavelet_filter()].
#' @return numeric vector of the same length as `x`.
#' @examples
#' x <- rnorm(32)
#' all.equal(idwt_periodic(dwt_periodic(x, 3, "coif1"), 3, "coif1"), x)
#' @export
dwt_periodic <- function(x, L, filter = "coif1") {
  h <- wavelet_filter(filter); n <- length(h)
  if (length(x) %% 2^L != 0) stop("length(x) must be divisible by 2^L.", call. = FALSE)
  C <- x; out <- numeric(0)
  H <- rev(h)
  G <- h; G[seq(1, n, by = 2)] <- -G[seq(1, n, by = 2)]
  for (j in seq_len(L)) {
    nn <- length(C)
    Cp <- c(C[((-(n - 1):-1) %% nn) + 1], C)          # periodic extension
    D  <- stats::convolve(Cp, rev(G), type = "open")
    D  <- D[seq(n, n + nn - 2, by = 2) + 1]
    Cn <- stats::convolve(Cp, rev(H), type = "open")
    C  <- Cn[seq(n, n + nn - 2, by = 2) + 1]
    out <- c(D, out)
  }
  c(C, out)
}

#' @rdname dwt_periodic
#' @param w wavelet coefficients from [dwt_periodic()].
#' @export
idwt_periodic <- function(w, L, filter = "coif1") {
  H <- wavelet_filter(filter); n <- length(H)
  G <- rev(H); G[seq(2, n, by = 2)] <- -G[seq(2, n, by = 2)]
  nn <- length(w); LL <- nn / 2^L
  C <- w[seq_len(LL)]
  for (j in seq_len(L)) {
    wrap <- ((0:(n / 2 - 1)) %% LL) + 1
    D  <- w[(LL + 1):(2 * LL)]
    Cu <- numeric(2 * LL + n); Du <- numeric(2 * LL + n)
    Cu[seq(1, 2 * LL + n, by = 2)] <- c(C, C[wrap])
    Du[seq(1, 2 * LL + n, by = 2)] <- c(D, D[wrap])
    Cn <- stats::convolve(Cu, rev(H), type = "open") + stats::convolve(Du, rev(G), type = "open")
    C  <- Cn[(n:(n + 2 * LL - 1)) - 1]
    LL <- 2 * LL
  }
  C
}

#' Threshold wavelet coefficients
#'
#' @description
#' Applies a thresholding rule to the detail coefficients of a periodic DWT.
#'
#' * `"hard"`: term-by-term hard thresholding, \eqn{\hat d\,1\{|\hat d| > \lambda\}},
#'   with \eqn{\lambda = c\,\hat\sigma\sqrt{\log(n)/n}} (Chesneau, Dewan and
#'   Doosti 2016; \eqn{c = 8} in their implementation).
#' * `"block"`: the BlockJS rule of Cai (1999) used by Shirazi and Doosti
#'   (2022): at each level the coefficients are grouped in blocks of length
#'   \eqn{L = \lfloor \log n \rfloor}; a block is kept if
#'   \eqn{\sum_{j \in B} \hat d_j^2 \ge \lambda^* L \hat\sigma^2}
#'   with \eqn{\lambda^* = 4.50524} (`rule = "truncate"`), or shrunk by
#'   \eqn{(1 - \lambda^* L \hat\sigma^2 / S_B^2)_+} (`rule = "james-stein"`).
#'
#' \eqn{\hat\sigma} is the MAD of the finest-level coefficients.
#'
#' @param w coefficient vector from [dwt_periodic()].
#' @param L number of DWT levels used to produce `w`.
#' @param n sample size (enters the threshold).
#' @param type `"hard"` or `"block"`.
#' @param const threshold constant: `c` for `"hard"`, \eqn{\lambda^*} for `"block"`.
#' @param block_length block length for `"block"`; default \eqn{\lfloor\log n\rfloor}.
#' @param rule for `"block"`: `"truncate"` (keep/kill) or `"james-stein"`.
#' @return thresholded coefficient vector.
#' @references
#' Cai, T. T. (1999). Adaptive wavelet estimation: a block thresholding and
#' oracle inequality approach. *Annals of Statistics*, 27, 898--924.
#' @export
threshold_coefficients <- function(w, L, n, type = c("hard", "block"),
                                   const = NULL, block_length = NULL,
                                   rule = c("truncate", "james-stein")) {
  type <- match.arg(type); rule <- match.arg(rule)
  nn <- length(w); finest <- w[(nn / 2 + 1):nn]
  sigma <- 1.4826 * stats::median(abs(finest - stats::median(finest)))
  smooth_n <- nn / 2^L
  if (type == "hard") {
    if (is.null(const)) const <- 8
    lambda <- const * sigma * sqrt(log(n) / n)
    idx <- (smooth_n + 1):nn
    w[idx] <- w[idx] * (abs(w[idx]) > lambda)
    return(w)
  }
  if (is.null(const)) const <- 4.50524
  if (is.null(block_length)) block_length <- max(1L, floor(log(n)))
  pos <- smooth_n
  for (lev in seq_len(L)) {
    len <- smooth_n * 2^(lev - 1)
    idx <- (pos + 1):(pos + len)
    Lb  <- min(block_length, len)
    blocks <- split(idx, ceiling(seq_along(idx) / Lb))
    for (b in blocks) {
      S2 <- sum(w[b]^2)
      thr <- const * length(b) * sigma^2
      w[b] <- if (rule == "truncate") w[b] * (S2 >= thr) else w[b] * max(0, 1 - thr / S2)
    }
    pos <- pos + len
  }
  w
}

#' Local-linear smoother with Gaussian kernel
#'
#' Post-smoothing of a wavelet estimate on a grid (Ramirez and Vidakovic 2010;
#' Fan 1992): for each `x` solves the weighted least-squares line through
#' `(grid, y)` with weights \eqn{K((grid - x)/h)} and returns its value at `x`.
#' Translated from the MATLAB routine `loc_lin` of B. Vidakovic.
#' @param x evaluation points.
#' @param grid,y tabulated function.
#' @param h bandwidth.
#' @return numeric vector.
#' @export
loclin_smooth <- function(x, grid, y, h) {
  vapply(x, function(x0) {
    d <- grid - x0
    k <- stats::dnorm(d / h)
    s1 <- sum(k * d); s2 <- sum(k * d^2)
    wgt <- k * (s2 - d * s1)
    sum(wgt * y) / sum(wgt)
  }, numeric(1))
}
