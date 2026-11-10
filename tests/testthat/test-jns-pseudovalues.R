# The public API must reproduce the Bernstein pseudo-values of the revenue
# paper (Doosti, Journal of Nonparametric Statistics, under review) exactly.
#
# The paper's estimator, as implemented in its replication script
# (R_code/JNS/jns_montecarlo.R, function bb_values()), is: sort the pooled bids,
# score them u_i = i / (N + 1), take q-hat from qdf_bernstein() with h = 1/m and
# zero_at_origin = FALSE at m = ceiling(sqrt(N)), and set
#     V_i = b_(i) + (1 - eta) / (n - 1) * u_i * q-hat(u_i).
#
# The reference below is a verbatim transcription of that function. If this test
# ever fails, the package no longer computes the published estimator through its
# documented arguments, and the paper's numbers would have to be regenerated.

bb_values_reference <- function(bids, eta, m) {
  n <- ncol(bids)
  x <- as.numeric(bids)
  o <- order(x)
  xs <- x[o]
  N <- length(xs)
  u <- seq_len(N) / (N + 1)
  qh <- qdf_bernstein(u, xs, 1 / m, zero_at_origin = FALSE)
  v <- numeric(N)
  v[o] <- xs + (1 - eta) / (n - 1) * u * qh
  matrix(v, nrow = nrow(bids))
}

test_that("fpa_values() reproduces the published Bernstein pseudo-values", {
  for (group in c("n2", "n3")) {
    b <- timber[[group]]
    m <- ceiling(sqrt(length(b)))
    for (eta in c(0.25, 0.60)) {
      ref <- bb_values_reference(b, eta, m)
      fit <- fpa_values(b, method = "bernstein", bandwidth = "sqrt",
                        zero_at_origin = FALSE, u_scheme = "i/(n+1)",
                        monotone = "none", trim = 0, eta = eta)

      # bandwidth = "sqrt" must be the m the paper fixes
      expect_equal(fit$qdf$m, m, info = paste(group, eta))

      # the recovered values themselves
      expect_equal(sort(as.numeric(fit$values)), sort(as.numeric(ref)),
                   tolerance = 1e-8, info = paste(group, eta))

      # and the auction layout, which fpa_revenue() depends on through the
      # per-auction maximum: equal revenue implies the values sit in the right
      # rows, not merely that the multiset agrees
      r_grid <- stats::quantile(b, c(0.25, 0.5, 0.75), names = FALSE)
      expect_equal(fpa_revenue(fit$values_mat, r = r_grid, eta = eta, v_s = 0),
                   fpa_revenue(ref, r = r_grid, eta = eta, v_s = 0),
                   tolerance = 1e-8, info = paste(group, eta))
    }
  }
})

test_that("the winning-bid check of Table 4 reproduces", {
  # Remark 3: at a reserve below the bid support, estimated revenue must equal
  # the expected winning bid. The published gaps are +15.81 (n = 2) and +16.60
  # (n = 3) at eta = 0.60, against mean winning bids of 87.63 and 98.47.
  targets <- list(n2 = c(mean_win = 87.63, gap = 15.81),
                  n3 = c(mean_win = 98.47, gap = 16.60))
  for (group in c("n2", "n3")) {
    b <- timber[[group]]
    m <- ceiling(sqrt(length(b)))
    mean_win <- mean(apply(b, 1, max))
    rev0 <- fpa_revenue(bb_values_reference(b, 0.60, m),
                        r = min(b) - 1, eta = 0.60, v_s = 0)
    expect_equal(mean_win, targets[[group]][["mean_win"]], tolerance = 1e-3)
    expect_equal(rev0 - mean_win, targets[[group]][["gap"]], tolerance = 1e-2)
  }
})
