#' @import dplyr
#' @import tidyr
#' @import ggpubr
#' @export
#'
#' @title
#' Compute empirical power using Monte Carlo methods
#' @description
#' This function compute the empirical power using Monte Carlo methods for a large number of tests for symmetry with known median. To see the available tests, see
#' details.
#' @return
#' A data frame with empirical powers for the test
#' @details
#' Compute empirical power on tribble or data frame that contains the distributions of the DLG.
#' Jk: Corzo, J. and Babativa, G. (2013),
#' R: McWiliams, P. (1990),
#' Rs: Baklizi, A. (2003)
#' Mp: Modarres and Gastwirth (1996).
#' @author Giovany Babativa <gbabativam@@gmail.com>
#' @param data data frame or tribble with name and lambda parameters of the DLG.
#' @param statis Test statistic to be used. By default \code{stat = c("Bk", "Jk")}
#' @param type Type of test. When the test is with know median, select \code{type = "k"} else \code{type = "u"} for unknown median.
#' @param alpha Size of test. By default \eqn{\alpha = 0.05}
#' @param nsize Sample Size.
#' @param rep Number of replicates for the Monte Carlo process.
#' @param plot draw empirical power function with results.
#' @references
#' Corzo, J., & Babativa, G. (2013). A modified runs test for symmetry. Journal of Statistical Computation and Simulation, 83(5), 984-991.
#' @examples
#' x <- rnorm(20)
#' #--- All test
#' (test <- symmetry_test(x))
#' #--- Choose any test
#' (Jk_test <- symmetry_test(x, stat = "Jk", Jk = 6))
#' #--- Choose severals tests
#' (MyTest1 <- symmetry_test(x, stat = c("Bk", "Jk"), Bk = 5, Jk = 6))
#' (MyTest2 <- symmetry_test(x, stat = c("Mp", "Jk"), Jk = 6, Mp = c(10, 20, 25)))

power_symmetry_test <- function(data, statis = c("Bk", "Jk", "R", "Rs", "Mp"), Bk = 5, Jk =6,
                  type = "k", alpha = 0.05, nsize, rep, plot = TRUE){
  ni <<- nsize
  B  <<- rep
  lsamples = data[,2:5] %>% purrr::pmap(samples)

  outPower <-  lapply(lsamples, function(x){
    result <- apply(x, 2, function(y){
      Stats = symmetry_test(y, statis = statis, Bk = Bk, Jk = Jk, type = type)
    })
    Stats <- dplyr::bind_rows(result) %>%
      mutate(reject =  ifelse(stat=="Bk" & stat.value - 1 < stats::qbinom(alpha, Bk, 0.5), 1,
                              ifelse(stat=="Bk" & stat.value - 1 == stats::qbinom(alpha, Bk, 0.5), (alpha- sum(stats::dbinom(0:stats::qbinom(alpha, Bk, 0.5)-1, Bk, 0.5)))/stats::dbinom(stats::qbinom(alpha, Bk, 0.5), Bk, 0.5),
                                     ifelse(stat == "R" & stat.value - 1 < stats::qbinom(alpha, ni-1, 0.5), 1,
                                            ifelse(stat == "R" & stat.value - 1 == stats::qbinom(alpha, ni-1, 0.5), (alpha - sum(stats::dbinom(0:stats::qbinom(alpha, ni-1, 0.5)-1, ni-1, 0.5)))/stats::dbinom(stats::qbinom(alpha, ni-1, 0.5), ni-1, 0.5),
                                                   ifelse(p.value < alpha, 1, 0))))))%>%
      group_by(stat) %>%
      summarise(power = mean(reject)) %>%
      pivot_wider(names_from = stat, values_from = power)
  })

  Power <- dplyr::bind_rows(outPower) %>%
    cbind(distributions[, "case"])

  if(plot){
    print(Power.plot(Power, ni))
  }
  return(Power)
}

