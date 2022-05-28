#' @import ggplot2
#' @import tidyr
#' @export
samples <- function(lambda1, lambda2, lambda3, lambda4){
  u <-  replicate(B, stats::runif(ni))
  xstar <-  lambda1 + ((u^lambda3 - (1-u)^lambda4)/lambda2)
  med <- lambda1 + (0.5^lambda3 - (1-0.5)^lambda4)/lambda2
  x <- xstar - med
  return(x)
}

Power.plot <- function(df, n){
  tmp <- df %>%
    pivot_longer(-case, names_to = "statistic", values_to = "power")%>%
    ggplot(aes(x=case, y=power, group=statistic)) +
    geom_line(aes(linetype=statistic, color=statistic)) +
    geom_point(aes(color=statistic))+
    scale_y_continuous(breaks = seq(0,1,0.1)) +
    labs(title = paste0("Empirical power with n=", n),
         xlab = "Distribution", ylab= "Empirical power")+
    theme_bw()
  tmp +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.title=element_text(face="bold",size="11"))

}

samplesq <- function(lambda1, lambda2, lambda3, lambda4){
  u <-  replicate(B, stats::runif(ni))
  xstar <-  lambda1 + ((u^lambda3 - (1-u)^lambda4)/lambda2)
  q <- lambda1 + (pti^lambda3 - (1-pti)^lambda4)/lambda2
  x <- xstar - q
  return(x)
}
