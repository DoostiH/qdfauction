## Submission of qdfauction 0.3.0

This is a feature release with three bug fixes. It adds procurement
(lowest-bid-wins) auctions, recovered by reflecting the bids and reusing the
existing inversion, and a moving-block bootstrap for auctions observed over
time. The Bernstein quantile density estimator gains a `zero_at_origin`
argument, and the package now ships `inst/CITATION`.

All defaults are unchanged, so results from version 0.2.0 reproduce exactly,
apart from the three bug fixes listed in NEWS.md. Those affect only the cases
described there: a numerically negligible fitted dependence parameter, a
vector of seller valuations, and a non-integer number of auctions passed to
the new block resampler.

## Test environments

* win-builder, R 4.6.2 (release) -- Status: OK
* win-builder, R Under development (2026-11-07 r90850) -- Status: OK
* local, Windows x86_64, R 4.6.1

## R CMD check results

0 errors | 0 warnings | 0 notes, on both win-builder checks.

The local check additionally reports two notes that do not arise on CRAN's
machines:

* "checking HTML version of manual ... Skipping checking math rendering:
  package 'V8' unavailable" -- V8 is not installed locally.

* "checking for detritus in the temp directory ... Found the following
  files/directories: 'lastMiKTeXException'" -- a MiKTeX temporary file left
  behind by the local build of the PDF manual.

The test suite is 541 expectations and runs in about 210 seconds on
win-builder. One test is skipped: it reproduces the pseudo-values of a
published paper against that paper's original reference code, which is not
distributed with the package.

## Note on the Author field

The Author and Maintainer fields are given explicitly rather than left for
R CMD build to generate, so that they match the form derived by R 4.6.2 and
later, which labels the ORCID and affiliation entries of a person() comment.
The package is built under R 4.6.1, whose generated fields omit those labels
and therefore differ from the derivation.

## Reverse dependencies

None. No package on CRAN depends on, imports or suggests qdfauction.
