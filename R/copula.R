# Copula module for affiliated private values -------------------------------
#
# Port of copqdf_core.R (Doosti 2026a, corrected version). Exchangeable
# n-dimensional Archimedean copulas: Clayton (theta > 0), Gumbel (theta >= 1),
# Frank (theta > 0), plus independence. All closed forms nest independence.

.copula_families <- c("independence", "clayton", "gumbel", "frank")

#' Copula multiplier and diagonal functions of exchangeable Archimedean copulas
#'
#' @description
#' For an exchangeable \eqn{n}-dimensional Archimedean copula \eqn{C}:
#' * `copula_psi()` is the copula multiplier of Doosti (2026a),
#'   \eqn{\psi(u) = C_2(u,\ldots,u)/c(u,\ldots,u)}, which enters the inverse bid
#'   function \eqn{v = b + \psi(u)\,q(u)}. Under independence
#'   \eqn{\psi(u) = u/(n-1)} and the GPV formula is recovered.
#' * `copula_delta()` is the diagonal \eqn{\delta(u) = C(u,\ldots,u)}, the
#'   distribution function of the maximum of the \eqn{n} uniform scores.
#' * `copula_delta_prime()` is its derivative.
#'
#' @param u numeric vector in \eqn{(0, 1)}.
#' @param family `"independence"`, `"clayton"`, `"gumbel"` or `"frank"`.
#' @param theta copula parameter (ignored for independence).
#' @param n number of bidders.
#' @return numeric vector.
#' @references
#' Doosti, H. (2026a). A copula-quantile density approach to affiliated
#' private value auctions. *Economics Letters*, 261, 112837.
#' \doi{10.1016/j.econlet.2026.112837}
#' @examples
#' u <- c(0.2, 0.5, 0.8)
#' copula_psi(u, "clayton", theta = 2, n = 3)
#' copula_psi(u, "independence", n = 3)      # u / (n - 1)
#' copula_delta(u, "gumbel", theta = 1.5, n = 3)
#' @export
copula_psi <- function(u, family = "independence", theta = NULL, n) {
  family <- match.arg(family, .copula_families)
  check_theta(family, theta)
  switch(family,
         independence = u / (n - 1),
         clayton = (n * u - (n - 1) * u^(theta + 1)) / ((n - 1) * (theta + 1)),
         gumbel  = {
           lu <- -log(u)
           n * u * lu / ((n - 1) * (theta - 1 + n^(1 / theta) * lu))
         },
         frank   = {
           a <- expm1(-theta * u); D <- expm1(-theta)
           y <- 1 + a^n / D^(n - 1)
           -y * a / ((n - 1) * theta * exp(-theta * u))
         })
}

#' @rdname copula_psi
#' @export
copula_delta <- function(u, family = "independence", theta = NULL, n) {
  family <- match.arg(family, .copula_families)
  check_theta(family, theta)
  switch(family,
         independence = u^n,
         clayton = (n * u^(-theta) - (n - 1))^(-1 / theta),
         gumbel  = u^(n^(1 / theta)),
         frank   = {
           a <- expm1(-theta * u); D <- expm1(-theta)
           -log1p(a^n / D^(n - 1)) / theta
         })
}

#' @rdname copula_psi
#' @export
copula_delta_prime <- function(u, family = "independence", theta = NULL, n) {
  family <- match.arg(family, .copula_families)
  check_theta(family, theta)
  switch(family,
         independence = n * u^(n - 1),
         clayton = {
           A <- n * u^(-theta) - (n - 1)
           n * A^(-1 / theta - 1) * u^(-theta - 1)
         },
         gumbel  = { m <- n^(1 / theta); m * u^(m - 1) },
         frank   = {
           a <- expm1(-theta * u); D <- expm1(-theta)
           y <- 1 + a^n / D^(n - 1)
           n * a^(n - 1) * exp(-theta * u) / (D^(n - 1) * y)
         })
}

check_theta <- function(family, theta) {
  if (family == "independence") return(invisible())
  if (is.null(theta) || !is.finite(theta)) stop("`theta` is required for family '", family, "'.", call. = FALSE)
  bad <- if (family == "gumbel") theta < 1 else theta <= 0
  if (bad) stop("`theta` out of range for family '", family, "'.", call. = FALSE)
  invisible()
}

# bivariate log copula densities (pseudo-ML)
biv_logc <- function(u, v, family, theta) {
  switch(family,
         clayton = log(1 + theta) - (1 + theta) * (log(u) + log(v)) -
           (2 + 1 / theta) * log(u^(-theta) + v^(-theta) - 1),
         gumbel  = {
           x <- (-log(u))^theta; y <- (-log(v))^theta; S <- x + y
           -S^(1 / theta) + log(S^(1 / theta) + theta - 1) +
             (1 / theta - 2) * log(S) +
             (theta - 1) * (log(-log(u)) + log(-log(v))) - log(u) - log(v)
         },
         frank   = {
           Dm <- -expm1(-theta)
           num <- log(theta) + log(Dm) - theta * (u + v)
           den <- Dm - (-expm1(-theta * u)) * (-expm1(-theta * v))
           num - 2 * log(abs(den))
         })
}

#' Pseudo-maximum-likelihood estimation of an Archimedean copula from bids
#'
#' Transforms bids to pooled Hazen scores \eqn{\hat u = \mathrm{rank}/(N+1)},
#' forms all within-auction pairs, and maximises the bivariate copula
#' log-likelihood over \eqn{\theta}. Only pairs within the same auction are
#' used (pairs straddling auctions are independent and would bias
#' \eqn{\hat\theta} toward independence).
#'
#' @param bids numeric matrix, one row per auction and one column per bidder.
#' @param family `"clayton"`, `"gumbel"` or `"frank"`.
#' @param bounds optional search interval for \eqn{\theta}.
#' @return a list with `theta`, `tau` (implied Kendall's tau), `family`,
#'   `loglik` and `n_pairs`.
#' @examples
#' set.seed(1)
#' U <- rcopula_archimedean(100, 3, "clayton", theta = 2)
#' copula_fit(U, "clayton")
#' @export
copula_fit <- function(bids, family = c("clayton", "gumbel", "frank"), bounds = NULL) {
  family <- match.arg(family)
  if (!is.matrix(bids) || ncol(bids) < 2) stop("`bids` must be an L x n matrix with n >= 2.", call. = FALSE)
  L <- nrow(bids); n <- ncol(bids); N <- length(bids)
  uh <- matrix(rank(bids) / (N + 1), nrow = L)
  cb <- utils::combn(n, 2)
  prs <- do.call(rbind, lapply(seq_len(ncol(cb)), function(k) cbind(uh[, cb[1, k]], uh[, cb[2, k]])))
  if (is.null(bounds))
    bounds <- switch(family, clayton = c(1e-3, 30), gumbel = c(1 + 1e-4, 15), frank = c(1e-3, 40))
  nll <- function(th) -sum(biv_logc(prs[, 1], prs[, 2], family, th))
  o <- stats::optimize(nll, bounds)
  list(theta = o$minimum, tau = copula_theta_to_tau(o$minimum, family), family = family,
       loglik = -o$objective, n_pairs = nrow(prs))
}

#' Convert between Kendall's tau and the copula parameter
#' @param tau,theta scalar values.
#' @param family `"clayton"`, `"gumbel"` or `"frank"`.
#' @return numeric scalar.
#' @examples
#' copula_tau_to_theta(0.5, "clayton")   # 2
#' copula_theta_to_tau(2, "gumbel")      # 0.5
#' @export
copula_tau_to_theta <- function(tau, family = c("clayton", "gumbel", "frank")) {
  family <- match.arg(family)
  if (tau <= 0) return(switch(family, clayton = 1e-2, gumbel = 1 + 1e-2, frank = 1e-2))
  switch(family,
         clayton = 2 * tau / (1 - tau),
         gumbel  = 1 / (1 - tau),
         frank   = stats::uniroot(function(th) copula_theta_to_tau(th, "frank") - tau,
                                  c(1e-4, 100), tol = 1e-8)$root)
}

#' @rdname copula_tau_to_theta
#' @export
copula_theta_to_tau <- function(theta, family = c("clayton", "gumbel", "frank")) {
  family <- match.arg(family)
  switch(family,
         clayton = theta / (theta + 2),
         gumbel  = 1 - 1 / theta,
         frank   = {
           D1 <- stats::integrate(function(t) t / expm1(t), 0, theta, rel.tol = 1e-10)$value / theta
           1 + 4 * (D1 - 1) / theta
         })
}

# positive stable, LT = exp(-s^alpha)
rpostable <- function(m, alpha) {
  Th <- stats::runif(m, 0, pi); W <- stats::rexp(m)
  (sin(alpha * Th) / sin(Th)^(1 / alpha)) * (sin((1 - alpha) * Th) / W)^((1 - alpha) / alpha)
}

# logarithmic series (Kemp's 2nd algorithm, Devroye 1986)
rlogseries <- function(m, p) {
  out <- integer(m); lg <- log1p(-p)
  for (i in seq_len(m)) {
    u <- stats::runif(1)
    if (u >= p) { out[i] <- 1L; next }
    q <- -expm1(stats::runif(1) * lg)
    out[i] <- if (u < q^2) as.integer(floor(1 + log(u) / log(q))) else if (u > q) 1L else 2L
  }
  out
}

#' Simulate exchangeable Archimedean copula scores
#'
#' Marshall--Olkin frailty samplers in base R.
#' @param L number of auctions (rows).
#' @param n number of bidders (columns).
#' @inheritParams copula_psi
#' @return an `L` by `n` matrix with uniform margins.
#' @examples
#' U <- rcopula_archimedean(500, 2, "frank", theta = 5)
#' cor(U, method = "kendall")[1, 2]        # about 0.46
#' @export
rcopula_archimedean <- function(L, n, family = "independence", theta = NULL) {
  family <- match.arg(family, .copula_families)
  if (family == "independence" ||
      (family == "clayton" && theta < 5e-3) || (family == "frank" && theta < 5e-3) ||
      (family == "gumbel" && theta < 1 + 5e-3))
    return(matrix(stats::runif(L * n), L, n))
  E <- matrix(stats::rexp(L * n), L, n)
  switch(family,
         clayton = { V <- stats::rgamma(L, shape = 1 / theta, rate = 1); (1 + E / V)^(-1 / theta) },
         gumbel  = { V <- rpostable(L, 1 / theta); exp(-(E / V)^(1 / theta)) },
         frank   = { V <- rlogseries(L, -expm1(-theta)); -log1p(expm1(-theta) * exp(-E / V)) / theta })
}
