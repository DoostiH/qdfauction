#' USFS timber sealed-bid auctions (homogenised bids)
#'
#' Homogenised per-unit bids (dollars per unit of timber) in the sealed-bid
#' subsample of the U.S. Forest Service timber auction data assembled by Lu
#' and Perrigne (2008) and deposited in the Journal of Applied Econometrics
#' Data Archive, split by the number of bidders. The bids are **not** the raw
#' bids: observed lot heterogeneity has been removed as described in Details.
#' The announced reserve price is treated as non-binding in the literature
#' using these data.
#'
#' @format A list with two numeric matrices, one row per auction and one
#'   column per bidder, in the order in which the archive files record them:
#' \describe{
#'   \item{n2}{107 auctions with 2 bidders (columns `bid1`, `bid2`)}
#'   \item{n3}{108 auctions with 3 bidders (columns `bid1`, `bid2`, `bid3`)}
#' }
#' @details
#' Each auction in the archive reports the per-unit bid of every bidder
#' together with two lot characteristics: the appraised per-unit value and the
#' timber volume. Observed heterogeneity in these characteristics is removed
#' before the structural analysis, as is standard for these data (Lu and
#' Perrigne 2008; Haile, Hong and Shum 2003), so that the methods of this
#' package can be applied to bids drawn from a common distribution within each
#' stratum.
#'
#' Within each stratum (\eqn{n = 2} and \eqn{n = 3} separately), the log bid is
#' regressed by ordinary least squares on
#' \eqn{X = (1, \log(\mathrm{appraisal}), \log(\mathrm{volume}))}, which is
#' constant within an auction, and every bid is rescaled to the average lot,
#' \deqn{\tilde b = b \exp(\bar X' \hat\beta - X' \hat\beta).}
#' The intercept cancels in the difference, so only the slope coefficients
#' matter. The regressions fit the log bids with \eqn{R^2 = 0.671} for
#' \eqn{n = 2} (slopes 0.642 and 0.043) and \eqn{R^2 = 0.554} for
#' \eqn{n = 3} (slopes 0.553 and 0.070). The stored values are rounded to six
#' decimal places.
#'
#' The script that rebuilds the object from the archive files is
#' `data-raw/timber.R` of the package sources; the archive files themselves
#' are not redistributed here.
#' @source Journal of Applied Econometrics Data Archive, replication files
#'   of Lu, J. and Perrigne, I. (2008). Estimating risk aversion from
#'   ascending and sealed-bid auctions: the case of timber auction data.
#'   *Journal of Applied Econometrics*, 23, 871--896. \doi{10.1002/jae.1032};
#'   data archive entry \doi{10.15456/jae.2022319.0720356598}.
#'   Please cite that paper when using the data.
#' @references
#' Haile, P. A., Hong, H. and Shum, M. (2003). Nonparametric tests for common
#' values at first-price sealed-bid auctions. Cowles Foundation Discussion
#' Paper 1445.
#' @examples
#' data(timber)
#' fit <- fpa_values(timber$n3, method = "bernstein", bandwidth = "bcv")
#' summary(fit$values - fit$bids)     # implied markups
"timber"
