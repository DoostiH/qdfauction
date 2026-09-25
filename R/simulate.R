# Equilibrium bids and simulation ---------------------------------------------

#' Quantile grid used for the equilibrium computations
#'
#' Three segments with refined tails, because \eqn{1/\psi(u)} is stiff as
#' \eqn{u \to 0} for every family and as \eqn{u \to 1} for Gumbel.
#' @param n_interior number of interior points.
#' @return sorted numeric vector in \eqn{(0, 1)}.
#' @keywords internal
#' @export
fpa_ugrid <- function(n_interior = 2001L) {
  sort(unique(c(
    exp(seq(log(1e-6), log(0.02), length.out = 250)),
    seq(0.02, 0.90, length.out = n_interior),
    1 - exp(seq(log(0.10), log(1e-6), length.out = 1200)))))
}

trapz <- function(x, y) sum(0.5 * diff(x) * (y[-1] + y[-length(y)]))

# Lambda(u) = int 1/psi, up to an additive constant, on grid U
make_lambda <- function(U, psiv) {
  g <- 1 / psiv
  c(0, cumsum(0.5 * diff(U) * (g[-1] + g[-length(g)])))
}

# markup m(u) solving m' = q - m/psi, m(u_{j0}) = 0, via a stable one-pass recursion
markup_grid <- function(U, qv, lam, j0 = 1L) {
  NU <- length(U); m <- numeric(NU)
  if (j0 < NU) for (j in (j0 + 1):NU) {
    K <- exp(lam[j - 1] - lam[j])
    m[j] <- K * m[j - 1] + 0.5 * (U[j] - U[j - 1]) * (K * qv[j - 1] + qv[j])
  }
  m
}

#' Equilibrium bid function of a symmetric first-price auction
#'
#' @description
#' Computes the symmetric equilibrium bid function in quantile space,
#' \eqn{\beta(u) = Q(u) - m(u)} with the markup solving
#' \eqn{m'(u) = q(u) - m(u)/\psi(u)}, \eqn{m(u_0) = 0} (Milgrom and Weber 1982;
#' Doosti 2026a). Under independence \eqn{\psi(u) = u/(n-1)} and this is the
#' Guerre--Perrigne--Vuong bid function
#' \eqn{\beta(v) = v - \int_0^v F(t)^{n-1}dt / F(v)^{n-1}}. A reserve price at
#' quantile `p0` sets \eqn{m(p_0) = 0}. For CRRA bidders with utility
#' \eqn{U(t) = t^{1-\eta}} the first-order condition is
#' \eqn{v = b + (1-\eta)\psi(u) q(u)} (Guerre, Perrigne and Vuong 2009), so the
#' markup solves the same ODE with \eqn{\psi} replaced by \eqn{(1-\eta)\psi};
#' under independence this gives Zincenko's (2024) risk-averse bid function
#' \eqn{\beta(v) = e_1 \int_0^v x F(x)^{e_2} f(x)\,dx / F(v)^{e_1}},
#' \eqn{e_1 = (n-1)/(1-\eta)}, \eqn{e_2 = (n+\eta-2)/(1-\eta)}.
#'
#' @param Q quantile function of private values (a function of `u`).
#' @param q quantile density of private values (a function of `u`); if `NULL`
#'   it is obtained by numerical differentiation of `Q`.
#' @param n number of bidders.
#' @inheritParams copula_psi
#' @param p0 screening quantile of the reserve price (default 0, no reserve).
#' @param eta CRRA risk-aversion parameter in \eqn{[0, 1)} (0: risk neutral).
#' @param U quantile grid; see [fpa_ugrid()].
#' @return a list with `U`, `bid` (\eqn{\beta(U)}), `value` (\eqn{Q(U)}),
#'   `markup` and a function `beta(v)` interpolating the bid at value `v`.
#' @examples
#' bf <- fpa_bid_function(Q = function(u) qgamma(u, 2, 10),
#'                        q = function(u) 1 / dgamma(qgamma(u, 2, 10), 2, 10), n = 5)
#' bf$beta(c(0.1, 0.2, 0.5))
#' # compare with the direct integral form for IPV
#' v <- 0.2; v - integrate(function(t) pgamma(t, 2, 10)^4, 0, v)$value / pgamma(v, 2, 10)^4
#' @export
fpa_bid_function <- function(Q, q = NULL, n, family = "independence", theta = NULL,
                             p0 = 0, eta = 0, U = fpa_ugrid()) {
  if (eta < 0 || eta >= 1) stop("`eta` must be in [0, 1).", call. = FALSE)
  if (is.null(q)) {
    eps <- 1e-6
    q <- function(u) (Q(pmin(u + eps, 1 - 1e-9)) - Q(pmax(u - eps, 1e-9))) /
      (pmin(u + eps, 1 - 1e-9) - pmax(u - eps, 1e-9))
  }
  psiv <- (1 - eta) * copula_psi(U, family, theta, n)
  lam  <- make_lambda(U, psiv)
  qv   <- q(U); Qv <- Q(U)
  j0   <- if (p0 <= 0) 1L else which.min(abs(U - p0))
  m    <- markup_grid(U, qv, lam, j0)
  bid  <- Qv - m
  beta <- function(v) stats::approx(Qv, bid, xout = v, rule = 2)$y
  list(U = U, bid = bid, value = Qv, markup = m, beta = beta, n = n,
       family = family, theta = theta, p0 = p0, eta = eta)
}

#' Simulate first-price sealed-bid auctions
#'
#' Draws private values with a given marginal quantile function and an
#' exchangeable Archimedean copula (independent private values by default),
#' and computes exact equilibrium bids with [fpa_bid_function()].
#'
#' @param L number of auctions.
#' @param n number of bidders per auction.
#' @inheritParams fpa_bid_function
#' @param bid_function optional precomputed output of [fpa_bid_function()].
#' @return a list with matrices `values` and `bids` (`L` by `n`) and the bid
#'   function used.
#' @examples
#' set.seed(1)
#' sim <- fpa_simulate(20, 5, Q = function(u) qlnorm(u, 0, 0.5))
#' str(sim$bids)
#' @export
fpa_simulate <- function(L, n, Q, q = NULL, family = "independence", theta = NULL,
                         p0 = 0, eta = 0, bid_function = NULL) {
  if (is.null(bid_function)) bid_function <- fpa_bid_function(Q, q, n, family, theta, p0, eta)
  U <- bid_function$U
  Uv <- rcopula_archimedean(L, n, family, theta)
  Uc <- pmin(pmax(Uv, U[1]), U[length(U)])
  vals <- matrix(Q(Uc), nrow = L)
  bids <- matrix(stats::approx(U, bid_function$bid, xout = Uc)$y, nrow = L)
  list(values = vals, bids = bids, bid_function = bid_function, n = n, L = L,
       eta = bid_function$eta)
}
