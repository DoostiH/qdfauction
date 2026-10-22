## =====================================================================
## data-raw/timber.R -- rebuild data/timber.rda
##
## The `timber` dataset holds the HOMOGENISED per-unit bids of the sealed-bid
## subsample of the U.S. Forest Service timber auctions (1979, western states)
## of Lu and Perrigne (2008): 107 auctions with two bidders and 108 with three.
##
## Source files (not redistributed with this package). Download the
## replication files of Lu, J. and Perrigne, I. (2008), "Estimating risk
## aversion from ascending and sealed-bid auctions: the case of timber auction
## data", Journal of Applied Econometrics 23, 871-896, from the JAE Data
## Archive (doi:10.15456/jae.2022319.0720356598), unpack lp-data.zip and put
## these two files in this directory:
##
##   SealedBid_2Bidders.txt   107 rows, 5 columns
##   SealedBid_3Bidders.txt   108 rows, 6 columns

## These two files are the copyright of their authors and are NOT part of this
## package: keep them out of the repository (data-raw/.gitignore lists them)
## and out of any redistributed source. Cite Lu and Perrigne (2008) when you
## use the data.
##
## The files have no header row and are whitespace separated. Columns, as
## documented in the archive's readme:
##
##   2 bidders:  volume  bid1  bid2        appraisal  n_bidders
##   3 bidders:  volume  bid1  bid2  bid3  appraisal  n_bidders
##
## Bids and the appraisal value are per unit of timber; volume is the total
## volume of the lot. The two ascending-auction files of the archive are not
## used here.
##
## Homogenisation (Lu and Perrigne 2008; Haile, Hong and Shum 2003). Within
## each stratum n, with X = (1, log appraisal, log volume) constant within an
## auction, ordinary least squares of the log bid on X over all bids of the
## stratum gives beta-hat, and every bid is rescaled to the average lot:
##
##     b-tilde = b * exp(Xbar' beta-hat - X' beta-hat).
##
## The intercept cancels in the difference, so only the slopes matter.
## Expected fit: R^2 = 0.6714 (n = 2) and 0.5540 (n = 3). The rebuilt bids
## agree with data/timber.rda to the six decimals at which the stored values
## are rounded (one value of the three-bidder matrix differs in the last bits
## of its double representation).
##
## Usage, with the working directory set to the package root:
##     source("data-raw/timber.R")
##
## The script prints the fitted slopes and R^2 of each stratum, checks that the
## homogenised bids no longer depend on the covariates, compares the result
## with the version currently in data/timber.rda, and then rewrites that file.
## =====================================================================

RAW_DIR   <- "data-raw"
SPEC <- list(
  n2 = list(file = "SealedBid_2Bidders.txt", n = 2L, auctions = 107L),
  n3 = list(file = "SealedBid_3Bidders.txt", n = 3L, auctions = 108L)
)
ROUND_DP <- 6L             # the stored values are rounded to six decimals

## ---- helpers --------------------------------------------------------------

## bids as an `auctions` x n matrix, and the auction-level covariates
read_stratum <- function(path, n, auctions) {
  raw <- utils::read.table(path, header = FALSE)
  if (ncol(raw) != n + 3L)
    stop(basename(path), ": expected ", n + 3L, " columns (volume, ", n,
         " bids, appraisal, number of bidders), found ", ncol(raw), ".",
         call. = FALSE)
  if (nrow(raw) != auctions)
    stop(basename(path), ": expected ", auctions, " auctions, found ",
         nrow(raw), ".", call. = FALSE)
  if (!all(raw[[n + 3L]] == n))
    stop(basename(path), ": the last column should be the number of bidders (",
         n, ").", call. = FALSE)

  bids <- as.matrix(raw[, 1L + seq_len(n), drop = FALSE])
  storage.mode(bids) <- "double"
  dimnames(bids) <- list(NULL, paste0("bid", seq_len(n)))
  volume <- as.numeric(raw[[1L]])
  appraisal <- as.numeric(raw[[n + 2L]])

  if (anyNA(bids) || anyNA(volume) || anyNA(appraisal))
    stop(basename(path), ": missing values.", call. = FALSE)
  if (any(bids <= 0) || any(volume <= 0) || any(appraisal <= 0))
    stop(basename(path), ": bids, volume and the appraisal value must be positive.",
         call. = FALSE)

  list(bids = bids,
       X = cbind(log_appraisal = log(appraisal), log_volume = log(volume)))
}

## rescale every bid of a stratum to its average lot
homogenise <- function(bids, X) {
  n  <- ncol(bids)
  y  <- log(as.vector(t(bids)))                    # stacked auction by auction
  Xb <- X[rep(seq_len(nrow(X)), each = n), , drop = FALSE]
  fit   <- stats::lm(y ~ Xb)
  beta  <- stats::coef(fit)[-1L]                   # slopes; the intercept cancels
  names(beta) <- colnames(Xb)
  shift <- as.vector(sweep(Xb, 2L, colMeans(Xb)) %*% beta)    # (X - Xbar)' beta
  out <- matrix(exp(y - shift), ncol = n, byrow = TRUE,
                dimnames = list(NULL, colnames(bids)))
  list(bids = out, r2 = summary(fit)$r.squared, beta = beta, Xb = Xb)
}

## ---- build ----------------------------------------------------------------

timber <- list()
for (nm in names(SPEC)) {
  s    <- SPEC[[nm]]
  path <- file.path(RAW_DIR, s$file)
  if (!file.exists(path))
    stop("Missing ", path, ". See the header of this script for where to get it.",
         call. = FALSE)

  raw <- read_stratum(path, s$n, s$auctions)
  hom <- homogenise(raw$bids, raw$X)
  message(sprintf("%s: %d auctions x %d bidders; R^2 = %.4f; slopes: %s",
                  nm, nrow(hom$bids), ncol(hom$bids), hom$r2,
                  paste(sprintf("%s = %.4f", names(hom$beta), hom$beta),
                        collapse = ", ")))

  ## the homogenised bids should no longer depend on the covariates
  left <- summary(stats::lm(log(as.vector(t(hom$bids))) ~ hom$Xb))$r.squared
  if (left > 1e-8)
    warning(nm, ": the homogenised bids still depend on the covariates (R^2 = ",
            signif(left, 3), ").", call. = FALSE)

  timber[[nm]] <- round(hom$bids, ROUND_DP)
}

## ---- compare with the version currently shipped ---------------------------

if (file.exists("data/timber.rda")) {
  shipped <- local({ e <- new.env(); load("data/timber.rda", envir = e); e$timber })
  for (nm in names(timber)) {
    if (!identical(dim(timber[[nm]]), dim(shipped[[nm]]))) {
      warning(nm, ": dimensions differ from the shipped data.", call. = FALSE)
      next
    }
    worst <- max(abs(timber[[nm]] - shipped[[nm]]))
    message(sprintf("%s: largest deviation from the shipped data: %.2e", nm, worst))
    if (worst > 10^(-ROUND_DP))
      warning(nm, ": the rebuilt bids differ from the shipped data by more than ",
              "the rounding of the stored values.", call. = FALSE)
  }
}

save(timber, file = "data/timber.rda", compress = "xz")
message("Wrote data/timber.rda")
