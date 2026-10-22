## Update

This is an update of qdfauction (0.1.0 -> 0.2.0). It fixes five problems in
the version currently on CRAN:

* `qdf()` fails for the wavelet methods when no `bandwidth` argument is given;
* the Frank copula objects (`copula_psi()`, `copula_delta()`,
  `copula_delta_prime()`) lose all precision to catastrophic cancellation at
  large `theta`, returning exactly 0 near u = 1 for Kendall's tau = 0.9;
* the tail integral of the expected-revenue estimator loses precision for
  strong dependence;
* `fpa_optimal_reserve_model()` reports interior optima below the 2nd
  percentile as corner solutions;
* the leave-one-out shortcut in the cross-validation criteria is invalid for
  estimators built from order statistics, so the criteria select the roughest
  bandwidth in the search range.

It also adds functions for the affiliated virtual value and a corner test for
the optimal reserve, and warns when a selected bandwidth lies at the end of its
search range or when the scale-dependent WBCV criterion is used. See NEWS.md.

## Test environments

* local Windows 11, R 4.6.1
* win-builder, R-release and R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

* Days since last update: 2. The update fixes the bugs listed above, which
  change results for users of the affected functions; the first makes a
  documented call fail outright.

The two further notes on the local machine ("Skipping checking math rendering:
package 'V8' unavailable" and a 'lastMiKTeXException' file left in the temp
directory) are local to that machine.