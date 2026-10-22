# Computational cost of the main functions of qdfauction as the number of
# pooled bids grows. Run with: source(system.file("benchmarks", "scaling.R",
# package = "qdfauction")). Times are elapsed seconds on one core.
library(qdfauction)
sizes <- c(200, 1000, 3000)
n_bidders <- 5
cases <- list(
  `qdf, Bernstein, BCV`                 = function(b, bm) qdf(b, "bernstein", "bcv"),
  `qdf, indirect Poisson, BCV`          = function(b, bm) qdf(b, "indirect_poisson", "bcv"),
  `qdf, block wavelet, fixed (j0, h)`   = function(b, bm) qdf(b, "wavelet_block", bandwidth = 0.15),
  `fpa_values, indirect Poisson, fixed h` = function(b, bm) fpa_values(bm, method = "indirect_poisson", bandwidth = 0.05),
  `fpa_density, local linear wavelet`   = function(b, bm) fpa_density(bm, method = "wavelet_ll"),
  `fpa_density, LSCV tuning (j0, h)`    = function(b, bm) fpa_density(bm, method = "wavelet_ll", j0 = "cv", j0_grid = 3:6, h_grid = c(0.05, 0.1, 0.2)),
  `fpa_density, Hickman-Hubbard`        = function(b, bm) fpa_density(bm, method = "kde_hh"),
  `fpa_density_ci, B = 199`             = function(b, bm) fpa_density_ci(bm, B = 199),
  `fpa_revenue_boot, Bernstein, B = 50` = function(b, bm) fpa_revenue_boot(bm, r = c(0.3, 0.5), B = 50, method = "bernstein", bandwidth = 0.05)
)
res <- matrix(NA_real_, length(cases), length(sizes), dimnames = list(names(cases), paste0("n=", sizes)))
for (s in seq_along(sizes)) {
  set.seed(1)
  sim <- fpa_simulate(sizes[s] / n_bidders, n_bidders, function(u) qbeta(u, 2, 2))
  for (k in seq_along(cases)) {
    t0 <- proc.time()[["elapsed"]]
    ok <- tryCatch({ cases[[k]](as.numeric(sim$bids), sim$bids); TRUE }, error = function(e) FALSE)
    res[k, s] <- if (ok) proc.time()[["elapsed"]] - t0 else NA
  }
}
print(round(res, 2))

# Functions for affiliated values, timed on values drawn from a Gumbel copula
# (theta = 1.5). The density bootstrap with the Marmer-Shneyerov estimator is
# the slowest case (several minutes at 3000 bids).
cases_apv <- list(
  `qdf_wavelet_nl, block`               = function(bm, fv) qdf_wavelet_nl(seq(0.01, 0.99, by = 0.01), as.numeric(bm), "block"),
  `fpa_values, nonparametric copula`    = function(bm, fv) fpa_values(bm, method = "poisson", bandwidth = 0.05, family = "nonparametric"),
  `copula_spec_test, B = 199`           = function(bm, fv) copula_spec_test(bm, "gumbel", B = 199),
  `fpa_density_boot, MS, B = 199`       = function(bm, fv) fpa_density_boot(bm, method = "marmer_shneyerov", B = 199),
  `fpa_optimal_reserve_qd, fixed h`     = function(bm, fv) fpa_optimal_reserve_qd(fv, v_s = 0, bandwidth = 0.05)
)
res_apv <- matrix(NA_real_, length(cases_apv), length(sizes),
                  dimnames = list(names(cases_apv), paste0("n=", sizes)))
for (s in seq_along(sizes)) {
  set.seed(1)
  sim <- fpa_simulate(sizes[s] / n_bidders, n_bidders, function(u) qbeta(u, 2, 2),
                      family = "gumbel", theta = 1.5)
  fv <- fpa_values(sim$bids, method = "poisson", bandwidth = 0.05, family = "gumbel")
  for (k in seq_along(cases_apv)) {
    t0 <- proc.time()[["elapsed"]]
    ok <- tryCatch({ cases_apv[[k]](sim$bids, fv); TRUE }, error = function(e) FALSE)
    res_apv[k, s] <- if (ok) proc.time()[["elapsed"]] - t0 else NA
  }
}
print(round(res_apv, 2))
