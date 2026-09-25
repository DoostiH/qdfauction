# Main interface ---------------------------------------------------------------

.qdf_methods <- c("indirect_poisson", "kernel", "kernel_corrected", "poisson",
                  "bernstein", "jones", "soni", "wavelet", "wavelet_hard",
                  "wavelet_block")

# default search ranges for m = 1/h, as used in the published simulations
.default_m_range <- function(method) {
  switch(method,
         jones = , soni = c(1, 50),
         c(10, 50))
}

#' Nonparametric estimation of a quantile density function
#'
#' @description
#' Estimates the quantile density function \eqn{q(u) = Q'(u) = 1/g(Q(u))}
#' (Parzen 1979; also called the sparsity function, Tukey 1965) from an i.i.d.
#' sample, with the smoothing parameter chosen by cross-validation. Seven
#' estimators are available; see [qdf_kernel()], [qdf_poisson()] and
#' [qdf_jones()] for their definitions and the simulation evidence in Doosti,
#' Dewan and Talebian (2025), where the indirect Poisson estimator with BCV
#' bandwidth performed best for first-price auction inversion.
#'
#' @param x numeric sample. The series estimators (`"poisson"`,
#'   `"bernstein"`, `"indirect_poisson"`) require strictly positive data.
#' @param method one of `"indirect_poisson"` (default), `"kernel"`,
#'   `"kernel_corrected"`, `"poisson"`, `"bernstein"`, `"jones"`, `"soni"`, or
#'   the wavelet estimators `"wavelet"` (linear), `"wavelet_hard"`,
#'   `"wavelet_block"` of [qdf_wavelet()], all with local-linear smoothing.
#' @param bandwidth either a selector name (`"bcv"`, `"rlcv"`, `"wbcv"`; for
#'   the wavelet methods `"cv"`, see [select_wavelet()]) or a numeric
#'   bandwidth \eqn{h} to use directly.
#' @param m_range search range for \eqn{m = 1/h} when a selector is used.
#'   Defaults to `c(10, 50)` (`c(1, 50)` for `"jones"` and `"soni"`), the
#'   ranges used in the published simulations.
#' @param u evaluation grid in \eqn{[0, 1]}. Defaults to \eqn{i/n}.
#' @param H second bandwidth for `method = "soni"` (the outer smoother on the
#'   probability scale). Defaults to the Sheather--Jones bandwidth of the
#'   uniform scores, `bw.SJ(rank(x)/n)`.
#' @param kernel kernel for the kernel-type estimators (`"gaussian"`,
#'   `"triangular"`, `"epanechnikov"`); see [qdf_kernel()].
#' @param j0 resolution level for the wavelet methods (default 5 when a
#'   numeric `h` is given; chosen by cross-validation with `bandwidth = "cv"`).
#' @param loo,grid,u_scheme passed to [select_bandwidth()].
#' @param ... currently unused.
#'
#' @return An object of class `"qdf"`: a list with components
#' \describe{
#'   \item{x}{the sorted sample}
#'   \item{u, q}{the evaluation grid and the estimated quantile density on it}
#'   \item{Q}{the empirical quantile function on `u` (for convenience)}
#'   \item{method, h, m, H}{the estimator and its smoothing parameters}
#'   \item{selection}{the output of [select_bandwidth()], or `NULL` when a
#'     numeric bandwidth was supplied}
#' }
#' with methods [predict.qdf()], `print()` and [plot.qdf()].
#'
#' @references
#' Doosti, H., Dewan, I. and Talebian, M. (2025). Nonparametric estimation of
#' private value distributions in first-price auctions: Evaluating quantile
#' density function approaches. *Economics Letters*, 257, 112670.
#' \doi{10.1016/j.econlet.2025.112670}
#'
#' Parzen, E. (1979). Nonparametric statistical data modeling. *Journal of the
#' American Statistical Association*, 74, 105--131.
#' @examples
#' set.seed(1)
#' x <- rgamma(80, 5, 1)
#' fit <- qdf(x, method = "bernstein", bandwidth = "bcv")
#' fit
#' predict(fit, u = c(0.1, 0.5, 0.9))
#' # compare with the truth 1 / dgamma(qgamma(u, 5), 5)
#' plot(fit, true = function(u) 1 / dgamma(qgamma(u, 5, 1), 5, 1))
#' @export
qdf <- function(x, method = c("indirect_poisson", "kernel", "kernel_corrected",
                              "poisson", "bernstein", "jones", "soni",
                              "wavelet", "wavelet_hard", "wavelet_block"),
                bandwidth = c("bcv", "rlcv", "wbcv", "cv"), m_range = NULL,
                u = NULL, H = NULL, kernel = "gaussian", j0 = NULL, loo = NULL,
                grid = 0L, u_scheme = NULL, ...) {
  method <- match.arg(method)
  cl <- match.call()
  if (startsWith(method, "wavelet")) return(qdf_wavelet_fit(x, method, bandwidth, u, j0, cl, ...))
  positive <- method %in% c("poisson", "bernstein", "indirect_poisson")
  x <- prep_sample(x, positive = positive)
  n <- length(x)
  if (is.null(u)) u <- default_u(n) else u <- check_u(u)

  if (method == "soni" && is.null(H)) {
    H <- stats::bw.SJ(seq_len(n) / n)
  }
  fun <- estimator_fun(method, H, kernel)

  if (is.numeric(bandwidth)) {
    if (length(bandwidth) != 1L || bandwidth <= 0)
      stop("`bandwidth` must be a single positive number or a selector name.", call. = FALSE)
    h <- bandwidth; sel <- NULL
  } else {
    bandwidth <- match.arg(bandwidth)
    if (is.null(m_range)) m_range <- .default_m_range(method)
    sel <- select_bandwidth(x, fun, bandwidth, m_range = m_range, loo = loo,
                            grid = grid, u_scheme = u_scheme)
    h <- sel$h
  }
  q <- fun(u, x, h)
  structure(list(call = cl, x = x, n = n, u = u, q = q,
                 Q = stats::quantile(x, u, type = 1, names = FALSE),
                 method = method, h = h, m = 1 / h, H = H, kernel = kernel,
                 selection = sel),
            class = "qdf")
}

# wavelet branch of qdf()
qdf_wavelet_fit <- function(x, method, bandwidth, u, j0, cl, ...) {
  x <- prep_sample(x); n <- length(x)
  if (is.null(u)) u <- default_u(n) else u <- check_u(u)
  thr <- switch(method, wavelet = "none", wavelet_hard = "hard", wavelet_block = "block")
  if (is.numeric(bandwidth)) {
    h <- bandwidth; if (is.null(j0)) j0 <- 5L; sel <- NULL
  } else {
    bandwidth <- match.arg(bandwidth, c("cv", "bcv", "rlcv", "wbcv"))
    if (bandwidth != "cv") stop("Wavelet methods use `bandwidth = \"cv\"` or a numeric h.", call. = FALSE)
    sel <- select_wavelet(x, threshold = thr, smooth = TRUE, ...)
    if (!is.null(j0)) sel$j0 <- j0
    h <- sel$h; j0 <- sel$j0; sel$selector <- "cv"; sel$criterion <- max(sel$ll, na.rm = TRUE)
    sel$m_range <- c(NA, NA)
  }
  q <- qdf_wavelet(u, x, j0 = j0, h = h, threshold = thr, smooth = TRUE, ...)
  structure(list(call = cl, x = x, n = n, u = u, q = q,
                 Q = stats::quantile(x, u, type = 1, names = FALSE),
                 method = method, h = h, m = 1 / h, H = NULL, j0 = j0,
                 threshold = thr, selection = sel, dots = list(...)),
            class = "qdf")
}

# build a (u, x, h) closure for a method
estimator_fun <- function(method, H = NULL, kernel = "gaussian") {
  switch(method,
         kernel            = function(u, x, h) qdf_kernel(u, x, h, kernel = kernel),
         kernel_corrected  = function(u, x, h) qdf_kernel_corrected(u, x, h, kernel = kernel),
         poisson           = qdf_poisson,
         bernstein         = qdf_bernstein,
         jones             = function(u, x, h) qdf_jones(u, x, h, kernel = kernel),
         soni              = function(u, x, h) qdf_soni(u, x, h, H = H, kernel = kernel),
         indirect_poisson  = qdf_indirect_poisson)
}

#' Predict from a fitted quantile density estimate
#'
#' @param object a `"qdf"` object.
#' @param u evaluation points in \eqn{[0, 1]}; defaults to the grid stored in
#'   `object`.
#' @param ... unused.
#' @return numeric vector of \eqn{\hat q(u)}.
#' @export
predict.qdf <- function(object, u = NULL, ...) {
  if (is.null(u)) return(object$q)
  u <- check_u(u)
  if (startsWith(object$method, "wavelet"))
    return(do.call(qdf_wavelet, c(list(u, object$x, j0 = object$j0, h = object$h,
                                       threshold = object$threshold, smooth = TRUE),
                                  object$dots)))
  estimator_fun(object$method, object$H, if (is.null(object$kernel)) "gaussian" else object$kernel)(u, object$x, object$h)
}

#' @export
print.qdf <- function(x, digits = 4, ...) {
  cat("Quantile density estimate\n")
  cat("  method    :", x$method, "\n")
  cat("  n         :", x$n, "\n")
  cat("  h (1/m)   :", format(x$h, digits = digits),
      sprintf(" (m = %s)", format(x$m, digits = digits)), "\n")
  if (!is.null(x$H)) cat("  H         :", format(x$H, digits = digits), "\n")
  if (!is.null(x$kernel) && x$method %in% c("kernel", "kernel_corrected", "jones", "soni"))
    cat("  kernel    :", x$kernel, "\n")
  if (!is.null(x$j0)) cat("  j0        :", x$j0, "\n")
  if (!is.null(x$selection))
    cat("  selector  :", x$selection$selector,
        sprintf(" (criterion = %s, m in [%s, %s])",
                format(x$selection$criterion, digits = digits),
                x$selection$m_range[1], x$selection$m_range[2]), "\n")
  else cat("  selector  : fixed bandwidth\n")
  cat("  q(u) on", length(x$u), "grid points; range",
      paste(format(range(x$q), digits = digits), collapse = " to "), "\n")
  invisible(x)
}

#' Plot a quantile density estimate
#'
#' @param x a `"qdf"` object.
#' @param true optional function of `u` giving the true quantile density, drawn
#'   for comparison.
#' @param ... further arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.qdf <- function(x, true = NULL, ...) {
  args <- list(x = x$u, y = x$q, type = "l", xlab = "u",
               ylab = expression(hat(q)(u)),
               main = paste("Quantile density,", x$method))
  args <- utils::modifyList(args, list(...))
  do.call(graphics::plot, args)
  if (!is.null(true)) {
    ug <- seq(min(x$u), max(x$u), length.out = 200)
    graphics::lines(ug, true(ug), lty = 2, col = 2)
    graphics::legend("topleft", c("estimate", "true"), lty = 1:2, col = 1:2, bty = "n")
  }
  invisible(x)
}
