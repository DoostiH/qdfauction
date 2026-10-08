## Resubmission

This is a resubmission. In response to the review of the first submission
(version 1.0.0), I have:

* Written every reference in the Description field of DESCRIPTION in the form
  authors (year) <doi:...>, or authors (year) <https:...> where no DOI exists
  (Jones 1992), so that all references are auto-linked.
* Added a \value section to the documentation of every exported method,
  including plot.fpa_density.Rd, plot.fpa_revenue_boot.Rd and
  plot.fpa_values.Rd, explaining that these methods are called for their
  side effect of drawing a plot and return the object invisibly. The
  \value sections of plot.qdf.Rd and plot.fpa_density_ci.Rd were expanded in
  the same way.

The version number is now 0.1.0. Because the package is not yet on CRAN, I
have taken the opportunity to submit the current development version, which
adds further estimators and inference tools (listed in NEWS.md) to the code
reviewed in the first submission.

## Test environments

* local Windows 11, R 4.6.1
* win-builder, R-release and R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

  The words flagged as possibly misspelled in DESCRIPTION (Chaubey,
  Chesneau, Dewan, Doosti, Hickman, Jain, Marmer, Shirazi, Shneyerov, Soni,
  Talebian, Zincenko) are author names in the cited references.

## Comments

Parts of the package are R translations or adaptations of code by other
authors (Chaubey, Dewan and Li; Hickman and Hubbard; Vidakovic; Marmer and
Shneyerov; Zincenko). They are listed as contributors and copyright holders
in Authors@R, and the origin of each file is documented in inst/COPYRIGHTS.
The authors have been informed, and those who could be reached have agreed
to the inclusion of their code.
