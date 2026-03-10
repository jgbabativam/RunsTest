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

power_symmetry_test <- function(data, statis = c("Bk", "Jk", "R", "Rs", "Mp"),
                                Bk = 5, Jk = 6, type = "k",
                                alpha = 0.05, nsize, rep, plot = TRUE, ncores = NULL) {
  
  ni <- nsize
  B  <- rep
  
  # --- Número de cores automático si no se especifica ---
  if (is.null(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)
  
  # --- Precalcular valores críticos UNA sola vez (fuera de las réplicas) ---
  crit <- list()
  if ("Bk" %in% statis) {
    for (k in Bk) {
      qb  <- stats::qbinom(alpha, k, 0.5)
      acc <- (alpha - sum(stats::dbinom(0:(qb - 1), k, 0.5))) /
        stats::dbinom(qb, k, 0.5)
      crit[[paste0("Bk", k)]] <- list(q = qb, acc = acc)
    }
  }
  if ("R" %in% statis) {
    qr  <- stats::qbinom(alpha, ni - 1, 0.5)
    acc <- (alpha - sum(stats::dbinom(0:(qr - 1), ni - 1, 0.5))) /
      stats::dbinom(qr, ni - 1, 0.5)
    crit[["R"]] <- list(q = qr, acc = acc)
  }
  
  # --- Función de rechazo vectorizada (sin ifelse anidados) ---
  compute_reject <- function(stat_name, stat_value, p_value) {
    sv <- stat_value - 1   # valor original sin el +1
    
    if (grepl("^Bk", stat_name)) {
      info <- crit[[stat_name]]
      return(ifelse(sv < info$q, 1,
                    ifelse(sv == info$q, info$acc, 0)))
    }
    if (stat_name == "R") {
      info <- crit[["R"]]
      return(ifelse(sv < info$q, 1,
                    ifelse(sv == info$q, info$acc, 0)))
    }
    # Para Jk, Rs, Mp: usar directamente p.value
    return(as.integer(p_value < alpha))
  }
  
  # --- Generar muestras ---
  lsamples <- data[, 2:5] %>% purrr::pmap(samples)
  
  # --- Clúster paralelo ---
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterExport(cl,
                          varlist = c("symmetry_testOP", "statis", "Bk", "Jk", "type",
                                      "alpha", "ni", "crit", "compute_reject"),
                          envir = environment())
  parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(randtests)
  })
  
  # --- Paralelizado por distribución ---
  outPower <- parallel::parLapply(cl, lsamples, function(x) {
    
    # Todas las réplicas de una distribución juntas
    result <- apply(x, 2, function(y) {
      symmetry_testOP(y, statis = statis, Bk = Bk, Jk = Jk, type = type)
    })
    
    dplyr::bind_rows(result) %>%
      group_by(stat) %>%
      mutate(reject = compute_reject(stat[1], stat.value, p.value)) %>%
      summarise(power = mean(reject), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = stat, values_from = power)
  })
  
  Power <- dplyr::bind_rows(outPower) %>%
    cbind(data[, "case"])
  
  if (plot) print(Power.plot(Power, ni))
  
  return(Power)
}


