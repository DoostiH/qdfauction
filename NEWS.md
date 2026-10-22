# qdfauction 0.2.0

This release combines the changes made for the revisions of Doosti (2026b,
2026c) with fixes from the wavelet branch.

## Bug fixes (these change numerical results)

* Frank copula objects (`copula_psi()`, `copula_delta()`,
  `copula_delta_prime()`) are evaluated without the catastrophic
  cancellation in `1 + a^n / D^(n-1)` at large `theta` (Kendall's tau above
  about 0.8); previously `psi` was exactly 0 near `u = 1` at tau = 0.9.
* The tail integral in the revenue engine of `fpa_revenue_qd()` and
  `fpa_optimal_reserve_qd()` is accumulated from the top; the previous forward
  difference lost precision for Clayton at `theta >= 2` (up to 0.05 in revenue
  at the 90th reserve percentile for tau = 0.75).
* `fpa_optimal_reserve_model()` searches from `p_lo = 1e-4` (was 0.02), so
  interior optima close to the lower support are no longer reported as corners.
* `select_bandwidth()` and `qdf()` use exact leave-one-out by default for
  every selector, including `"wbcv"`. The scaling approximation `q (n-1)/n`, invalid
  for estimators built from order statistics, is used only with `loo = FALSE`.
* `bcv_criterion()`, `rlcv_criterion()` and `wbcv_criterion()` can be called
  directly with their default `u_scheme` (previously an error), and
  `wbcv_criterion()` now uses exact leave-one-out by default (`loo = TRUE`), as
  `select_bandwidth()` does.
* `qdf()` with a wavelet method and no `bandwidth` argument now uses
  cross-validation (`bandwidth = "cv"`) instead of failing.
* The `print()` method of `fpa_values()` no longer shows `theta = NULL` for
  the nonparametric copula multiplier.

## Documentation and data

* The documentation of the `timber` dataset now states that the bids are
  homogenised, not raw: within each stratum the log bids were regressed on the
  log appraised per-unit value and the log timber volume and rescaled to the
  average lot (R-squared 0.671 for n = 2 and 0.554 for n = 3). The data are
  unchanged.
* `data-raw/timber.R` now rebuilds the dataset from the Journal of Applied
  Econometrics Data Archive files of Lu and Perrigne (2008)
  (doi:10.15456/jae.2022319.0720356598), applying that regression, instead of
  storing the bids as literal values.

## New features

* `fpa_virtual_value()` and `fpa_virtual_value_model()`: the affiliated
  virtual value `omega(u) = Q(u) - q(u) H(u)`.
* `fpa_corner_test()` (with a `print()` method): a crossing test for whether
  any binding reserve price is optimal (Doosti 2026c). The decision is the
  point rule on the full sample, `min over u in [0.10, 0.98] of omega_hat(u)
  - v_s >= 0`; a bootstrap over auctions gives the frequency of that decision.
  The copula `family` is a required argument, because the result depends on it.
* `qdf(bandwidth = "sqrt")` (and `fpa_values()`, the revenue functions and
  `fpa_revenue_boot()` through it): the deterministic rule
  `m = ceiling(sqrt(N))` for the Bernstein estimator, used in Doosti (2026b).
  It lies inside the range covered by that paper's asymptotic theory.
* `qdf_bernstein()` gains `zero_at_origin`; `qdf_bernstein_loo()` computes
  exact leave-one-out Bernstein estimates by an index shift in O(nK) and is
  used automatically by the cross-validation criteria.
* `select_bandwidth(grid = "exhaustive")` returns the global minimiser over all
  integer `m` in `m_range` (not the default).
* `fpa_montecarlo()` gains `p0` and `eta`, passed to the simulator.
* Stage 3 of `fpa_revenue_qd()`, `fpa_optimal_reserve_qd()`,
  `fpa_virtual_value()` and `fpa_corner_test()`: if the data-driven bandwidth
  for the pseudo-values is minimised at an end of `m_range`, the quantile
  density is refitted with the published range `m_range = c(10, 50)` and
  `stats::optimize()` (Doosti 2026c), with a warning; if that fit is also at an
  end, the function stops and asks for a numeric bandwidth.
* Argument checks: `select_bandwidth()` and `qdf()` accept only
  `grid = "exhaustive"` or a non-negative integer (a mistyped value used to fall
  back silently to `optimize()`); `qdf_bernstein_loo()` checks that `u` has one
  increasing value per observation, that the data are positive and that
  `floor(1/h) >= 2`.
* Documentation: `qdf()` lists all ten estimators and gives the correct pages
  for Parzen (1979); the default grid of `fpa_density_ci()` is described
  correctly; `qdf_wavelet_nl()` explains when post-smoothing helps.
* `inst/benchmarks/scaling.R` also times the functions for affiliated values,
  the specification test and the density bootstrap.

## Changed defaults

* `fpa_revenue_boot()` resamples whole auctions by default
  (`resample = "auctions"`; previously `"bids"`). The revenue estimator uses the
  highest pseudo-value in each auction, and resampling pooled bids breaks that
  structure: in the design of Zincenko (2024), with the Bernstein estimator and
  `m = ceiling(sqrt(N))`, pooled resampling gave pointwise coverage of about
  0.65-0.80 at nominal 0.90, against 0.85-0.89 with auction resampling.
  `"bids"` is kept for comparison. `fpa_corner_test()` already resamples
  auctions; `fpa_density_boot()` keeps `"bids"`, because the density
  estimators do not use the auction structure.
* `bandwidth = "wbcv"`, `select_bandwidth(selector = "wbcv")` and
  `wbcv_criterion()` are kept only for backward compatibility, are not
  recommended, and now warn. The two terms of the WBCV criterion have
  different units, so the selected bandwidth depends on the units of the data,
  and on data of order one or larger the criterion tends to select the
  smoothest fit in `m_range`. The criterion is being revised (Doosti 2026b);
  it is kept so that existing code runs. Examples, the README and the vignette
  now use BCV.
* `select_bandwidth()` and `qdf()` keep `stats::optimize()` over `m_range` as
  the default search for every selector, as in the published simulations;
  `grid = "exhaustive"` (every integer `m` in `m_range`, the global minimiser)
  is available as an option. A user-supplied `optimizer` takes precedence over
  the exhaustive search.
* Exact leave-one-out is the default for `"wbcv"` in `select_bandwidth()` and
  `qdf()` (previously the scaling approximation); see Bug fixes.
* `fpa_revenue_boot()` and the bootstrap of `fpa_corner_test()` select the
  bandwidth once, on the full sample, and hold it fixed in every bootstrap
  replication, as the theory assumes (Doosti 2026b, Theorems 6-7); previously
  the bandwidth was re-selected in each replicate. Results with a data-driven
  `bandwidth` change; results with a numeric `bandwidth` do not. The selected
  `m` and whether it lies at a boundary are returned and printed.
* `select_bandwidth()` detects boundary solutions for every selector and search
  method: if the minimiser lies at either end of `m_range`, it warns and
  returns `at_boundary = TRUE` (and `boundary = "lower"` or `"upper"`), and
  `print()` of a `qdf()` fit reports it. A boundary minimum means that the
  criterion did not select a bandwidth in that range.

# qdfauction 0.1.0

First CRAN release.

## Quantile density estimation and auction methods

* Quantile density estimation: ten estimators (kernel, corrected kernel,
  Poisson, Bernstein, Jones, Soni, indirect Poisson, linear / hard- /
  block-thresholded wavelets) with BCV, RLCV, WBCV and cross-validated
  wavelet tuning (`qdf()`).
* Private-value recovery in first-price auctions for independent or
  Archimedean-copula affiliated values and risk-neutral or CRRA bidders
  (`fpa_values()`), with GPV, boundary-corrected and Zincenko kernel
  benchmarks and monotonicity options.
* Private-value density estimation with wavelet, Hickman-Hubbard,
  Marmer-Shneyerov and empirical-quantile second stages
  (`fpa_density()`), tuned by interior least-squares or log-likelihood
  cross-validation (`select_fpa_density()`).
* Inference for the GPV estimator of the private-value density: analytic
  standard errors, pointwise confidence intervals and bootstrap uniform
  confidence bands (Ma, Marmer and Shneyerov 2019) (`fpa_density_ci()`).
* Expected revenue, optimal reserve price and bootstrap inference
  (`fpa_revenue()`, `fpa_optimal_reserve()`, `fpa_revenue_boot()`), exact
  equilibrium bids and simulation (`fpa_bid_function()`, `fpa_simulate()`),
  Monte Carlo and evaluation tools (`fpa_montecarlo()`,
  `fpa_evaluate_values()`, `fpa_evaluate_density()`).
* Extensibility: user-defined kernels via `make_kernel()`, with
  `check_kernel()` reporting mass, symmetry, order and roughness; a
  user-supplied optimiser for bandwidth selection (`optimizer` in
  `select_bandwidth()` and `qdf()`).
* Bid matrices may be supplied as data frames; incomplete auctions, missing
  values and univariate/multivariate mix-ups give informative errors.
* `inst/benchmarks/scaling.R` times the main functions as the sample grows.
* USFS timber auction data (`timber`) and a vignette.

## Further estimators and inference

* Expected revenue and optimal reserve under affiliated private values by the
  copula--quantile density functional (`fpa_revenue_qd()`,
  `fpa_optimal_reserve_qd()`; Doosti 2026c).
* Nonlinear wavelet estimators of the quantile density built on genuine
  wavelet coefficients, with the theoretical thresholds of Chesneau, Dewan and
  Doosti (2016, Theorems 3.2 and 3.4) and the block thresholding of Shirazi
  and Doosti (2022) (`qdf_wavelet_nl()`), on a periodised orthonormal basis
  (`phi_per()`, `psi_per()`, `psi_jk()`).
* The hard-thresholding estimator of Theorem 2 of Doosti, Dewan and Sbai
  (2026) with level-dependent thresholds (`fpa_density(method = "wavelet_theory")`).
* Flexible, nonparametric copula multiplier (`copula_psi_np()`,
  `fpa_values(family = "nonparametric")`) and a bootstrap specification test
  for parametric copulas (`copula_spec_test()`) (Doosti 2026a). The
  semiparametric estimator of Hubbard, Li and Paarsch (2012) is documented as
  `fpa_values(method = "gpv", family = ...)`.
* Bootstrap pointwise intervals and uniform bands for any private-value
  density estimator, including Marmer and Shneyerov (2012)
  (`fpa_density_boot()`).
* `qdf()` is now scale-equivariant for every estimator: it estimates on
  `x / sd(x)` and rescales (`rescale = TRUE`). This matters for the indirect
  Poisson, Jones and Soni estimators, whose smoothing parameter lives on the
  data scale; `rescale = FALSE` reproduces the earlier behaviour. The default
  evaluation grid is now `i/(n+1)`.
* `select_wavelet()` (and `qdf(..., bandwidth = "cv")` for the wavelet
  estimators) uses least-squares cross-validation by default
  (`criterion = "lscv"`) and skips bandwidths below the lattice spacing;
  cross-validation arguments are no longer passed on to the estimator.
