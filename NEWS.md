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
