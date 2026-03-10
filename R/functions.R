#' @import ggplot2
#' @import tidyr
#' @export
#'
#' @title
#' Generate random samples from the Generalised Lambda Distribution (GLD)
#' @description
#' Generates a matrix of random samples from the Generalised Lambda Distribution (GLD)
#' using the inverse CDF method. Each column corresponds to one Monte Carlo replicate,
#' centred at the theoretical median.
#' @return
#' A numeric matrix of dimension \code{ni x B}, where \code{ni} is the sample size and
#' \code{B} is the number of replicates. Each column is one centred random sample.
#' @param lambda1 Location parameter \eqn{\lambda_1} of the GLD.
#' @param lambda2 Scale parameter \eqn{\lambda_2} of the GLD.
#' @param lambda3 First shape parameter \eqn{\lambda_3} of the GLD.
#' @param lambda4 Second shape parameter \eqn{\lambda_4} of the GLD.
#' @param ni Sample size per replicate.
#' @param B Number of Monte Carlo replicates.
#' @details
#' The samples are centred by subtracting the theoretical median of the GLD.
#' @examples
#' \dontrun{
#' mat <- samples(lambda1 = 0, lambda2 = 0.1975, lambda3 = 0.1349, lambda4 = 0.1349,
#'                ni = 20, B = 100)
#' dim(mat)  # 20 x 100
#' }
samples <- function(lambda1, lambda2, lambda3, lambda4, ni, B){
  u      <- replicate(B, stats::runif(ni))
  xstar  <- lambda1 + ((u^lambda3 - (1-u)^lambda4) / lambda2)
  med    <- lambda1 + (0.5^lambda3 - (1-0.5)^lambda4) / lambda2
  x      <- xstar - med
  return(x)
}

#' @keywords internal
Power.plot <- function(df, n){
  tmp <- df %>%
    pivot_longer(-case, names_to = "statistic", values_to = "power") %>%
    ggplot2::ggplot(ggplot2::aes(x = case, y = power, group = statistic)) +
    ggplot2::geom_line(ggplot2::aes(linetype = statistic, color = statistic)) +
    ggplot2::geom_point(ggplot2::aes(color = statistic)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    ggplot2::labs(title = paste0("Empirical power with n=", n),
                  x = "Distribution", y = "Empirical power") +
    ggplot2::theme_bw()
  tmp +
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 90, hjust = 1),
                   plot.title   = ggplot2::element_text(hjust = 0.5, size = 14),
                   axis.title   = ggplot2::element_text(face = "bold", size = 11))
}

#' @keywords internal
samplesq <- function(lambda1, lambda2, lambda3, lambda4, ni, B, pti){
  u      <- replicate(B, stats::runif(ni))
  xstar  <- lambda1 + ((u^lambda3 - (1-u)^lambda4) / lambda2)
  q      <- lambda1 + (pti^lambda3 - (1-pti)^lambda4) / lambda2
  x      <- xstar - q
  return(x)
}
