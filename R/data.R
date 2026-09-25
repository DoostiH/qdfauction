#' USFS timber sealed-bid auctions
#'
#' Per-unit bids (dollars per unit of timber) in the sealed-bid subsample of
#' the U.S. Forest Service timber auction data assembled by Lu and Perrigne
#' (2008) and deposited in the Journal of Applied Econometrics Data Archive,
#' split by the number of bidders. The announced reserve price is treated as
#' non-binding in the literature using these data.
#'
#' @format A list with two numeric matrices, one row per auction:
#' \describe{
#'   \item{n2}{107 auctions with 2 bidders (columns `bid1`, `bid2`)}
#'   \item{n3}{108 auctions with 3 bidders (columns `bid1`, `bid2`, `bid3`)}
#' }
#' @details The bids were extracted from the replication data of Lu and
#'   Perrigne (2008), which record U.S. Forest Service timber sales; only the
#'   per-unit bids of the sealed-bid auctions with two and three bidders are
#'   retained, pooled by the number of bidders. The script that builds the
#'   object is in `data-raw/timber.R` of the package sources.
#' @source Journal of Applied Econometrics Data Archive, replication files
#'   of Lu, J. and Perrigne, I. (2008). Estimating risk aversion from
#'   ascending and sealed-bid auctions: the case of timber auction data.
#'   *Journal of Applied Econometrics*, 23, 871--896. \doi{10.1002/jae.1032}.
#'   Please cite that paper when using the data.
#' @examples
#' data(timber)
#' fit <- fpa_values(timber$n3, method = "bernstein", bandwidth = "bcv")
#' summary(fit$values - fit$bids)     # implied markups
"timber"
