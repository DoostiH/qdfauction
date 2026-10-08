This is a resubmission. In response to the review of the first
submission (version 1.0.0), I have:

* Written every reference in the Description field in the form
  authors (year) <doi:...>, or authors (year) <https:...> where no
  DOI exists, so that all references are auto-linked.
* Added a \value section to every exported method, including
  plot.fpa_density.Rd, plot.fpa_revenue_boot.Rd and plot.fpa_values.Rd,
  explaining that they are called for their side effect of drawing a
  plot and return the object invisibly.

As the package is not yet on CRAN, the version is now 0.1.0 and
includes the current development version (see NEWS.md).

R CMD check results (local Windows R 4.6.1, win-builder R-release and
R-devel): 0 errors | 0 warnings | 1 note (new submission).

The words flagged as possibly misspelled in DESCRIPTION (Chaubey,
Chesneau, Dewan, Doosti, Marmer, Shirazi, Shneyerov, Soni, Talebian,
Zincenko) are author names in the cited references.