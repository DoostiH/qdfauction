# qdfauction 1.0.0

First public release.

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
