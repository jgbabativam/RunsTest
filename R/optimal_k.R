# =============================================================================
# Calculation of the optimal trimming parameter k* for the Ck statistic
#
# CRITERION: k* = max( ceil(c / (1-q)),  ceil(gamma * n) )
#
#   Component 1 - STABILITY: ceil(c / (1-q))
#     Under worst-case misspecification m0 = Q_q, n1k ~ Bin(k, 1-q).
#     Requires E[n1k] >= c for sufficient support in the conditional
#     runs distribution (c=5, q=0.70 gives bound = 17).
#
#   Component 2 - POWER: ceil(gamma * n)
#     Retains at least fraction gamma of the sample.
#     Simulation with n in {20,30,50,100} shows gamma=0.60 reproduces
#     the empirically optimal k* for n >= 30.
#     For n < 30 the stability bound dominates (k*=17).
# =============================================================================


#' Optimal trimming parameter k* for the Ck statistic
#'
#' Computes the optimal number of observations \eqn{k^*} to retain when
#' computing the \eqn{C_k} conditional runs test for symmetry.
#' The criterion balances stability of the conditional distribution
#' (under worst-case median misspecification) against power retention.
#'
#' @param n      Sample size (positive integer).
#' @param q      Upper quantile of maximum misspecification considered
#'   (default 0.70). Must satisfy \code{0 < q < 1}.
#' @param c      Minimum required value for \eqn{E[n_{1k}]} under
#'   misspecification (default 5).
#' @param gamma  Minimum fraction of the sample to retain (default 0.60).
#'   Must satisfy \code{0 < gamma <= 1}.
#' @param verbose Logical. If \code{TRUE}, prints a detailed diagnostic
#'   table around \eqn{k^*}. Default \code{FALSE}.
#'
#' @return An integer giving \eqn{k^*}.
#'
#' @details
#' The criterion is
#' \deqn{k^* = \max\!\left(\lceil c/(1-q) \rceil,\; \lceil \gamma n \rceil\right)}
#' capped at \eqn{n}. The first component ensures that, when the true median
#' equals the \eqn{q}-th percentile, the expected number of observations above
#' \eqn{m_0} in the trimmed subsample is at least \eqn{c}. The second
#' component guarantees that at least a fraction \eqn{\gamma} of the sample
#' contributes to the test statistic.
#'
#' @references
#' Baklizi, A. (2003). A conditional distribution free runs test for
#' symmetry. \emph{Journal of Nonparametric Statistics}, \bold{15}(6),
#' 713--718. \doi{10.1080/10485250310001634737}
#'
#' @examples
#' cat("k* for n=30:", optimal_k(30), "\n")
#' cat("k* for n=50:", optimal_k(50), "\n")
#'
#' # Detailed diagnostic for study sample sizes
#' optimal_k(30,  verbose = TRUE)
#' optimal_k(50,  verbose = TRUE)
#' optimal_k(100, verbose = TRUE)
#'
#' # Table for a range of sample sizes
#' print(optimal_k_table(c(20, 25, 30, 40, 50, 75, 100, 150, 200)),
#'       row.names = FALSE)
#' @export
optimal_k <- function(n, q = 0.70, c = 5, gamma = 0.60, verbose = FALSE) {

  p_above <- 1 - q                        # P(X > m0) under incorrect H0
  k_stab  <- ceiling(c / p_above)         # stability bound
  k_power <- ceiling(gamma * n)           # power bound
  k_star  <- min(max(k_stab, k_power), n) # combined criterion

  if (verbose) {
    cat(sprintf("\n=== Optimal k* for n = %d ===\n", n))
    cat(sprintf("  q     = %.2f  (maximum misspecification percentile)\n", q))
    cat(sprintf("  p     = %.2f  (P(X > m0) under incorrect H0)\n", p_above))
    cat(sprintf("  c     = %.0f   (minimum required E[n1k])\n", c))
    cat(sprintf("  gamma = %.2f  (minimum sample fraction)\n", gamma))
    cat(sprintf("  k_stability = ceil(%.0f / %.2f) = %d\n", c, p_above, k_stab))
    cat(sprintf("  k_power     = ceil(%.2f * %d) = %d\n", gamma, n, k_power))
    cat(sprintf("  k* = max(%d, %d) = %d  (%.0f%% of n)\n\n",
                k_stab, k_power, k_star, 100 * k_star / n))

    cat(sprintf("%5s  %6s  %8s  %12s  %12s\n",
                "k", "k/n", "E[n1k]", "P(n1k<=2)", "P(n1k<=3)"))
    cat(strrep("-", 52), "\n")
    k_lo <- max(1L, k_star - 5L)
    k_hi <- min(n,  k_star + 5L)
    for (k in k_lo:k_hi) {
      marker <- if (k == k_star) " <-- k*" else ""
      cat(sprintf("%5d  %5.0f%%  %8.2f  %12.4f  %12.4f%s\n",
                  k, 100 * k / n, k * p_above,
                  stats::pbinom(2, k, p_above),
                  stats::pbinom(3, k, p_above),
                  marker))
    }
    cat("\n")
  }

  return(k_star)
}


#' Table of optimal k* values for a vector of sample sizes
#'
#' Applies \code{\link{optimal_k}} to each element of \code{n_vec} and
#' returns a summary \code{data.frame}.
#'
#' @param n_vec  Integer vector of sample sizes.
#' @param q      Upper quantile of maximum misspecification (default 0.70).
#' @param c      Minimum required \eqn{E[n_{1k}]} (default 5).
#' @param gamma  Minimum sample fraction to retain (default 0.60).
#'
#' @return A \code{data.frame} with columns \code{n}, \code{k_star},
#'   \code{prop} (k*/n as a percentage string), \code{E_n1k}, and
#'   \code{P_n1k_le2}.
#'
#' @examples
#' print(optimal_k_table(c(20, 30, 50, 100)), row.names = FALSE)
#' @export
optimal_k_table <- function(n_vec, q = 0.70, c = 5, gamma = 0.60) {
  p_above <- 1 - q
  results <- lapply(n_vec, function(n) {
    k <- optimal_k(n, q = q, c = c, gamma = gamma, verbose = FALSE)
    data.frame(
      n         = n,
      k_star    = k,
      prop      = paste0(round(100 * k / n), "%"),
      E_n1k     = round(k * p_above, 2),
      P_n1k_le2 = round(stats::pbinom(2, k, p_above), 4)
    )
  })
  do.call(rbind, results)
}
