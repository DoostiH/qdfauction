# User-defined kernels ---------------------------------------------------------

#' Construct and check a kernel for the kernel-type estimators
#'
#' @description
#' `make_kernel()` turns a kernel density function into the object used by
#' [qdf_kernel()], [qdf_kernel_corrected()], [qdf_jones()], [qdf_soni()] and
#' [qdf()] (argument `kernel`). The distribution function needed by the
#' boundary-corrected and Soni estimators is computed numerically.
#' `check_kernel()` reports the properties that matter for these estimators:
#' total mass, symmetry, support, order (the first non-vanishing moment),
#' \eqn{\int k^2} and whether the kernel takes negative values; it warns when a
#' property required by the theory fails.
#'
#' @param k a vectorised function, the kernel density.
#' @param name a label.
#' @param support half-width of an interval \eqn{[-s, s]} containing the
#'   support (or, for kernels with unbounded support, essentially all of the
#'   mass).
#' @param n_grid number of grid points for the numerical integration.
#' @param kernel a kernel name, a function, or an object from `make_kernel()`.
#' @return `make_kernel()`: an object of class `"qdf_kernel"` (a list with the
#'   density `d`, the distribution function `p`, `name` and the `check`
#'   results). `check_kernel()`: a list of properties, printed invisibly.
#' @examples
#' # the biweight (quartic) kernel
#' biweight <- make_kernel(function(u) 15 / 16 * (1 - u^2)^2 * (abs(u) <= 1), "biweight")
#' biweight
#' x <- rgamma(100, 5)
#' qdf_kernel(c(0.25, 0.5, 0.75), x, h = 0.05, kernel = biweight)
#' check_kernel("epanechnikov")
#' @export
make_kernel <- function(k, name = "custom", support = 1, n_grid = 20001L) {
  if (!is.function(k)) stop("`k` must be a function.", call. = FALSE)
  s <- 1.5 * support
  xg <- seq(-s, s, length.out = n_grid)
  kg <- k(xg)
  if (length(kg) != n_grid || any(!is.finite(kg)))
    stop("`k` must be vectorised and finite on [-", s, ", ", s, "].", call. = FALSE)
  cg <- c(0, cumsum((kg[-1] + kg[-n_grid]) / 2 * diff(xg)))
  cg <- cg / cg[n_grid]
  keep_dim <- function(y, x) { if (!is.null(dim(x))) dim(y) <- dim(x); y }
  p <- function(x) keep_dim(stats::approx(xg, cg, xout = as.numeric(x), yleft = 0, yright = 1)$y, x)
  d <- function(x) keep_dim(k(x), x)
  out <- structure(list(d = d, p = p, name = name, support = support), class = "qdf_kernel")
  out$check <- check_kernel(out, verbose = FALSE)
  out
}

#' @rdname make_kernel
#' @param verbose print the report.
#' @export
check_kernel <- function(kernel, verbose = TRUE) {
  K <- kernel_fun(kernel)
  s <- if (!is.null(K$support)) 1.5 * K$support else 8
  xg <- seq(-s, s, length.out = 4001); kg <- K$d(xg)
  sup <- range(xg[abs(kg) > 1e-12])
  lo <- sup[1] - 2 * s / 4000; hi <- sup[2] + 2 * s / 4000   # integrate over the support only
  f <- function(g) stats::integrate(g, lo, hi, subdivisions = 500L, rel.tol = 1e-8)$value
  mass <- f(K$d)
  xs <- seq(0, s, length.out = 2001)
  asym <- max(abs(K$d(xs) - K$d(-xs)))
  mom <- vapply(1:4, function(j) f(function(u) u^j * K$d(u)), numeric(1))
  ord <- which(abs(mom) > 1e-6)[1]; if (is.na(ord)) ord <- NA_integer_
  res <- list(name = K$name, mass = mass, max_asymmetry = asym, order = ord,
              roughness = f(function(u) K$d(u)^2), negative = any(kg < -1e-12),
              support = sup)
  if (abs(mass - 1) > 1e-3) warning("Kernel integrates to ", signif(mass, 4), ", not 1.", call. = FALSE)
  if (asym > 1e-6) warning("Kernel is not symmetric.", call. = FALSE)
  if (verbose) {
    cat("Kernel:", res$name, "\n")
    cat(sprintf("  integral       : %.6f\n  symmetric      : %s\n  order          : %s\n  int k^2        : %.4f\n  negative values: %s\n  support        : [%.3g, %.3g]\n",
                mass, asym <= 1e-6, if (is.na(ord)) "> 4" else ord, res$roughness, res$negative, sup[1], sup[2]))
  }
  invisible(res)
}

#' @export
print.qdf_kernel <- function(x, ...) { check_kernel(x, verbose = TRUE); invisible(x) }
