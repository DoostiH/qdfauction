# qdfauction

Nonparametric quantile density estimation and structural analysis of
first-price auctions.

**Version 1.0.0.** See the vignette
`vignette("timber", package = "qdfauction")` for an end-to-end analysis of
the USFS timber auctions. Quantile density
estimation, private-value recovery, private-value density estimation,
expected revenue with bootstrap inference, and a simulation/evaluation
toolkit that regenerates the Monte Carlo tables of the underlying papers.

## Installation

```r
# from the source tarball
install.packages("qdfauction_1.0.0.tar.gz", repos = NULL, type = "source")
```

## What it does

`qdf()` estimates the quantile density function q(u) = Q'(u) = 1/g(Q(u))
with one of seven estimators and a data-driven smoothing parameter:

| `method`             | Estimator                                     | Reference |
|----------------------|-----------------------------------------------|-----------|
| `kernel`             | direct kernel smoother of Q_n                 | Jones (1992) |
| `kernel_corrected`   | boundary-corrected kernel                     | Chaubey, Dewan & Li (2024) |
| `poisson`            | Poisson-weight smoother of Q_n                | Chaubey, Dewan & Li (2021) |
| `bernstein`          | Bernstein-polynomial smoother of Q_n          | Chaubey, Dewan & Li (2021) |
| `jones`              | indirect, 1/ĝ(Q_n(u))                         | Jones (1992) |
| `soni`               | double-smoothed indirect (integrated or Riemann form; optional right-censoring via Kaplan–Meier) | Soni, Dewan & Jain (2012) |
| `indirect_poisson`   | inverted Poisson-smoothed CDF (default)       | Chaubey, Dewan & Li (2024) |
| `wavelet`            | linear wavelet projection + local-linear smoothing | Chesneau, Dewan & Doosti (2016) |
| `wavelet_hard`       | hard-thresholded wavelet + smoothing           | Chesneau, Dewan & Doosti (2016) |
| `wavelet_block`      | block-thresholded (BlockJS) wavelet + smoothing | Shirazi & Doosti (2022) |

The kernel-type estimators accept `kernel = "gaussian"`, `"triangular"` or
`"epanechnikov"`.

Bandwidth selectors: `"bcv"` (biased CV), `"rlcv"` (robust likelihood CV,
Wu 2019) and `"wbcv"` (weighted BCV targeting u·q(u), Doosti 2026). The
wavelet estimators choose the resolution level `j0` and smoothing bandwidth
`h` by cross-validated held-out log-likelihood (`bandwidth = "cv"`).

The wavelet machinery (`phi_jk()`, `dwt_periodic()`,
`threshold_coefficients()`, `loclin_smooth()`) is a dependency-free R port
of the MATLAB/WaveLab routines used in the papers, and reproduces Table 1 of
Chesneau, Dewan and Doosti (2016) within Monte Carlo error.

```r
library(qdfauction)
set.seed(1)
x   <- rgamma(100, 5, 1)
fit <- qdf(x, method = "indirect_poisson", bandwidth = "bcv")
fit
predict(fit, u = c(0.1, 0.5, 0.9))
plot(fit, true = function(u) 1 / dgamma(qgamma(u, 5, 1), 5, 1))
```

## First-price auctions

`fpa_values()` inverts the equilibrium bid function through the quantile
density of bids, v = b + λ⁻¹{ψ(u) q(u)}:

* independent private values (ψ = u/(n−1), the GPV formula), or affiliated
  values with a Clayton, Gumbel or Frank copula (Doosti 2026a), with the
  copula parameter estimated by pseudo-ML on within-auction pairs;
* risk-neutral or CRRA bidders (`eta`);
* the conventional kernel benchmarks `method = "gpv"` / `"gpv_bc"`;
* monotonicity via the natural-spline fix of the 2025 paper or isotonic
  regression.

`fpa_simulate()` draws values from any marginal and copula, for risk-neutral
or CRRA bidders (`eta`), and computes exact equilibrium bids by the markup
recursion in quantile space (`fpa_bid_function()`).

```r
set.seed(1)
sim <- fpa_simulate(L = 50, n = 3, Q = function(u) (1 - 8 * u / 9)^(-1/2),
                    family = "clayton", theta = 2)
fit <- fpa_values(sim$bids, method = "indirect_poisson", family = "clayton")
fit; plot(fit, true = sim$values)
```

`fpa_density()` estimates the private-value density in two steps: pseudo-values
(Hickman–Hubbard boundary-corrected GPV by default, or any `fpa_values()`
fit) and then a second-stage density estimator — the local-linear, linear,
hard- or block-thresholded wavelet estimators of Doosti, Dewan and Sbaï
(2026), or the competitors `kde_hh` (Hickman–Hubbard 2015), `kde`,
`empirical_quantile` and `marmer_shneyerov` (2012). Dollar-scale data are
mapped to [0,1] by a log map and densities returned through the Jacobian.
`select_fpa_density()` chooses (j0, h) by held-out log-likelihood and
`fpa_density_cv()` compares estimators the same way. The USFS timber data of
Lu and Perrigne (2008) are included as `timber`.

```r
data(timber)
fd <- fpa_density(timber$n3, method = "wavelet_ll", j0 = "cv")
plot(fd)
fpa_density_cv(timber$n3, methods = c("wavelet_ll", "kde_hh", "marmer_shneyerov"))
```

## Expected revenue and counterfactuals

`fpa_revenue()` is the plug-in estimator of the seller's expected revenue at
reserve price r from recovered values (Zincenko 2024; Doosti 2026b), for
risk-neutral or CRRA bidders; `fpa_optimal_reserve()` maximises it.
`fpa_revenue_boot()` gives bootstrap standard errors, percentile-t
pointwise intervals and a uniform band. `fpa_revenue_model()` and
`fpa_optimal_reserve_model()` compute the true revenue and optimal reserve
of a known model (any marginal, any Archimedean copula), for use as the
benchmark in simulations. `method = "zincenko"` in `fpa_values()` provides
Zincenko's kernel-BC pseudo-values.

```r
sim <- fpa_simulate(50, 3, Q = function(u) 10 * u)
rb  <- fpa_revenue_boot(sim$bids, r = 3:8, B = 500, v_s = 2,
                        method = "bernstein", bandwidth = "wbcv")
plot(rb); fpa_revenue_model(function(u) 10 * u, n = 3, r = 3:8, v_s = 2)
```

## Simulation and evaluation

`fpa_montecarlo()` repeats simulate → estimate → score for a list of
estimators; `fpa_evaluate_values()` (MISE, log-likelihood, KS, MSEP) and
`fpa_evaluate_density()` (MISE, tail MISE, expected log-likelihood) are the
criteria used in the papers. Together with `fpa_simulate()` these regenerate
the simulation tables of Doosti, Dewan and Talebian (2025), Doosti (2026a,
2026b) and Doosti, Dewan and Sbaï (2026).

The low-level estimators (`qdf_kernel()`, `qdf_bernstein()`, ...) and
criteria (`bcv_criterion()`, ...) are exported for use with a fixed
bandwidth or custom selection.

## Reproducibility

The test suite (`tests/testthat/test-legacy-agreement.R`) checks every
estimator and criterion against the original simulation code of Doosti,
Dewan and Talebian (2025), kept verbatim in `helper-legacy.R`. `test-soni-original.R` does the same against the original 2012 code of
Soni, Dewan and Jain. Two corrections to the 2025 code are applied and documented there: `n` is defined
before use in the corrected-kernel estimator, and the RLCV criteria for the
Soni and indirect-Poisson estimators evaluate the full grid rather than a
stale loop index.

## References

Chaubey, Y. P., Dewan, I., Li, J. (2021). Sankhya B 83, S114–S139.
Chaubey, Y. P., Dewan, I., Li, J. (2024). Comm. Statist. Theory Methods 53, 5553–5573.
Chesneau, C., Dewan, I., Doosti, H. (2016). Comput. Statist. Data Anal. 94, 161–174.
Doosti, H. (2026a). A copula-quantile density approach to affiliated private value auctions. Economics Letters 261, 112837.
Doosti, H. (2026b). Estimating seller's expected revenue in first-price sealed-bid auctions via quantile density functions. Working paper, SSRN 7071380.
Doosti, H., Dewan, I., Sbaï, E. (2026). Adaptive wavelet estimation of private-value densities in first-price auctions. Working paper.
Doosti, H., Dewan, I., Talebian, M. (2025). Economics Letters 257, 112670.
Guerre, E., Perrigne, I., Vuong, Q. (2000). Econometrica 68, 525–574.
Hickman, B. R., Hubbard, T. P. (2015). J. Appl. Econometrics 30, 739–762.
Lu, J., Perrigne, I. (2008). J. Appl. Econometrics 23, 871–896.
Marmer, V., Shneyerov, A. (2012). J. Econometrics 167, 345–357.
Zincenko, F. (2024). J. Econometrics 241, 105734.
Jones, M. C. (1992). Ann. Inst. Statist. Math. 44, 721–727.
Shirazi, E., Doosti, H. (2022). Comm. Statist. Simulation Comput. 51, 539–553.
Soni, P., Dewan, I., Jain, K. (2012). Comput. Statist. Data Anal. 56, 3876–3886.
Wu, X. (2019). J. Bus. Econom. Statist. 37, 761–770.
