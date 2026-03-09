#' @import dplyr
#' @importFrom randtests druns
#' @export
#'
#' @title
#' Perform symmetry test
#' @description
#' This function contains a large number of tests for symmetry with known median. To see the available tests, see
#' details.
#' @return
#' A data frame with statistical values and p.values
#' @details
#' Perform symmetry tests on numeric vectors or objects, choose test with vector stat.
#' Jk: Corzo, J. and Babativa, G. (2013),
#' R: McWiliams, P. (1990),
#' Rs: Baklizi, A. (2003)
#' Mp: Modarres and Gastwirth (1996).
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param x vector of numeric information.
#' @param statis Test statistic to be used. By default \code{stat = c("Bk", "Jk")}
#' @param Bk Parameter of \eqn{B_k} statistic. By default \code{Bk = 5}
#' @param Jk Parameter of \eqn{J_k} statistic. By default \code{Jk = 6}
#' @param median The centre parameter around which to test symmetry. By default \code{median=0}
#' @param type Type of test. When the test is with know median, select \code{type = "k"} else \code{type = "u"} for unknown median.
#' @references
#' Corzo, J., & Babativa, G. (2013). A modified runs test for symmetry. Journal of Statistical Computation and Simulation, 83(5), 984-991.
#' @examples
#' x <- rnorm(20)
#' #--- All test
#' (test <- symmetry_test(x))
#' #--- Choose any test
#' (Jk_test <- symmetry_test(x, statis = "Jk", Jk = 6))
#' #--- Choose severals tests
#' (MyTest1 <- symmetry_testOP(x, statis = c("Bk", "Jk"),  Bk = c(5, 7), Jk = c(6, 8)))
#' (MyTest2 <- symmetry_testOP(x, statis = c("Mp", "Jk"), Jk = 6, Mp = c(10, 20, 25)))


symmetry_testOP <- function(x, statis = c("Bk", "Jk", "R", "Rs", "Mp", "Bkc"), type = "k",
                            Bk = 5, Jk = 6, Bkc = 11, Mp = c(10, 20, 25), median = 0) {
  
  if (type == "u") median <- median(x)
  
  n <- length(x)
  
  # --- Precálculos ---
  xc   <- x - median
  xabs <- abs(xc)
  si   <- (sign(xc) + 1) / 2
  ord  <- order(xabs)
  si   <- si[ord]
  indj <- seq_len(n)
  j    <- n - indj + 1
  Ij   <- c(1L, as.integer(diff(si) != 0))
  
  df <- data.frame(indj = indj, j = j, si = si, Ij = Ij)
  
  make_row <- function(stat, stat.value, p.value) {
    data.frame(stat       = stat,
               stat.value = round(stat.value, 1),
               p.value    = round(p.value, 4),
               stringsAsFactors = FALSE)
  }
  
  stats <- list()
  
  # --- Bk (múltiples cortes) ---
  if ("Bk" %in% statis) {
    for (k in Bk) {
      mask <- df$j <= k
      svBk <- sum(df$Ij[mask])
      stats[[paste0("Bk", k)]] <- make_row(paste0("Bk", k),
                                           svBk + 1,
                                           sum(stats::dbinom(0:svBk, k, 0.5)))
    }
  }
  
  # --- Bkc ---
  if ("Bkc" %in% statis) {
    sub      <- df[df$j <= Bkc, ]
    svBkc    <- sum(sub$Ij)
    n1       <- sum(sub$si);  n0 <- nrow(sub) - n1
    stat.val <- svBkc + 1
    stats[["Bkc"]] <- make_row("Bkc",
                               stat.val,
                               sum(druns(seq_len(stat.val) - 1, n1, n0)))
  }
  
  # --- Jk (múltiples cortes) ---
  if ("Jk" %in% statis) {
    for (k in Jk) {
      mask <- df$j <= k
      svJk <- sum(df$indj[mask] * df$Ij[mask])
      EJk  <- k * (2*n - k + 1) / 4
      VJk  <- k * (6*n^2 + 6*n + 2*k^2 - 3*k - 6*n*k + 1) / 24
      stats[[paste0("Jk", k)]] <- make_row(paste0("Jk", k),
                                           svJk + 1,
                                           stats::pnorm((svJk - EJk) / sqrt(VJk)))
    }
  }
  
  # --- R ---
  if ("R" %in% statis) {
    svR <- sum(df$Ij[df$indj >= 2])
    stats[["R"]] <- make_row("R",
                             svR + 1,
                             sum(stats::dbinom(0:svR, n - 1, 0.5)))
  }
  
  # --- Rs ---
  if ("Rs" %in% statis) {
    svRs <- sum(df$Ij)
    n1   <- sum(df$si);  n0 <- n - n1
    ERs  <- 1 + (2*n0*n1) / n
    VRs  <- (2*n0*n1) * (2*n0*n1 - n) / (n^2 * (n - 1))
    stats[["Rs"]] <- make_row("Rs",
                              svRs,
                              stats::pnorm((svRs - ERs) / sqrt(VRs)))
  }
  
  # --- Mp ---
  if ("Mp" %in% statis) {
    for (p in Mp) {
      np   <- trunc(n * p / 100)
      mask <- df$indj >= np + 2
      svMp <- sum((df$indj[mask] - np) * df$Ij[mask])
      q    <- 1 - p/100
      EMp  <- (n*q - 1) * (n*q + 2) / 4
      VMp  <- (n*q - 1) * (2*n^2*q^2 + 5*n*q + 6) / 24
      stats[[paste0("M", p)]] <- make_row(paste0("M", p),
                                          svMp,
                                          stats::pnorm((svMp - EMp) / sqrt(VMp)))
    }
  }
  
  dplyr::bind_rows(stats)
}



