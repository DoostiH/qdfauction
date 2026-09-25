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
* Expected revenue, optimal reserve price and bootstrap inference
  (`fpa_revenue()`, `fpa_optimal_reserve()`, `fpa_revenue_boot()`), exact
  equilibrium bids and simulation (`fpa_bid_function()`, `fpa_simulate()`),
  Monte Carlo and evaluation tools (`fpa_montecarlo()`,
  `fpa_evaluate_values()`, `fpa_evaluate_density()`).
* USFS timber auction data (`timber`) and a vignette.
