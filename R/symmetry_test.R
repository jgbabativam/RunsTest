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
#' (MyTest1 <- symmetry_test(x, statis = c("Bk", "Jk"), Bk = 5, Jk = 6))
#' (MyTest2 <- symmetry_test(x, statis = c("Mp", "Jk"), Jk = 6, Mp = c(10, 20, 25)))

symmetry_test <- function(x, statis = c("Bk", "Jk", "R", "Rs", "Mp", "Bkc"), type = "k",
                          Bk = 5, Jk = 6, Bkc = 11, Mp = c(10, 20, 25), median = 0){

  if(type == "u") median <- median(x)

  df <- data.frame(x) %>%
        mutate(xabs = abs(x-median),
               si = (sign(x-median)+1)/2) %>%
        arrange(xabs) %>%
        mutate(indj = row_number(),
               Ij = ifelse(indj == 1, 1,
                           ifelse(lag(si)==si,0, 1)),
               j = length(x)-indj+1)

  n <- length(x)
  stats = list()

  if("Bk" %in% statis) {
    stat <- "Bk"
    svBk  <- df %>%
             dplyr::filter(j <= Bk ) %>% dplyr::select(Ij) %>%
             sum()
    stat.value <- svBk + 1
    p.value <- sum(stats::dbinom(0:svBk, Bk, 0.5))
    out <- data.frame(stat = stat, stat.value=round(stat.value,1), p.value = round(p.value,4))
    out$stat <- as.character(out$stat)
    stats[[1]] <- out
  }

  if("Bkc" %in% statis) {
    stat <- "Bkc"
    temp  <- df %>%
             dplyr::filter(j <= Bkc) %>%
             dplyr::select(si, Ij) %>%
             summarise(svBkc = sum(Ij), n1 = sum(si), n0 = n() - n1)

    stat.value <- temp$svBkc + 1
    n1 <- temp$n1
    n0 <- temp$n0

    p.value <- sum(druns(0:stat.value - 1, n1, n0))
    out <- data.frame(stat = stat, stat.value=round(stat.value,1), p.value = round(p.value,4))
    out$stat <- as.character(out$stat)
    stats[[1]] <- out
  }

  if("Jk" %in% statis) {
    stat <- "Jk"
    svJk  <- df %>%
             dplyr::filter(j <= Jk ) %>% mutate(Ji = indj*Ij) %>% dplyr::select(Ji) %>%
             sum()
    EJk <- 1/4 * Jk * (2 * n - Jk + 1)
    VJk <- 1/24 * Jk * (6 * n^2 + 6 * n + 2 * Jk^2 - 3 * Jk - 6 * n * Jk + 1)

    stat.value <- svJk + 1
    p.value <- stats::pnorm((svJk-EJk)/sqrt(VJk))
    out <- data.frame(stat = stat, stat.value=round(stat.value,1), p.value = round(p.value,4))
    out$stat <- as.character(out$stat)
    stats[[2]] <- out
  }

  if("R" %in% statis) {
    stat <- "R"
    svR  <- df %>%
      dplyr::filter(indj >=  2) %>% dplyr::select(Ij) %>%
      sum()

    stat.value <- svR + 1
    p.value <- sum(stats::dbinom(0:svR, n-1, 0.5))
    out <- data.frame(stat = stat, stat.value=round(stat.value,1), p.value = round(p.value,4))
    out$stat <- as.character(out$stat)
    stats[[3]] <- out
  }

  if("Rs" %in% statis) {
    stat <- "Rs"
    res  <-  df %>%
             dplyr::select(si, Ij) %>%
             summarise(svRs = sum(Ij), n1 = sum(si))

    stat.value <- res$svRs
    n1 <- res$n1
    n0 <- n - n1

    ERs <- 1 + (2 * n0 * n1)/ n
    VRs <-  (2 * n0 * n1)*(2 * n0 * n1 - n)/(n^2 * (n - 1))

    p.value <- stats::pnorm((stat.value-ERs)/sqrt(VRs))

    out <- data.frame(stat = stat, stat.value=round(stat.value,1), p.value = round(p.value,4))
    out$stat <- as.character(out$stat)
    stats[[4]] <- out
  }

  if("Mp" %in% statis) {
    ps <- length(Mp)
    it <- 4
    for (p in Mp) {
     stat <- paste0("M", p)

     np <- trunc(n * p/100)
     svMp  <- df %>%
              dplyr::filter(indj >=  np + 2) %>%
              mutate(phiIk = (indj - np) * Ij) %>%
              dplyr::select(phiIk) %>%
              sum()

     EMp <- 1/4 * (n * (1 - p/100) - 1)*(n * (1- p/100) + 2)
     VMp <- 1/24 * (n * (1 - p/100) - 1)*(2 * n^2 * (1 - p/100)^2 + 5 * n * (1 - p/100) + 6)

     p.value <- stats::pnorm((svMp-EMp)/sqrt(VMp))
     out <- data.frame(stat = stat, stat.value=round(svMp,1), p.value = round(p.value,4))
     out$stat <- as.character(out$stat)
     it <- it + 1
     stats[[it]] <- out
    }
  }
  out <- dplyr::bind_rows(stats)
  return(out)

}





