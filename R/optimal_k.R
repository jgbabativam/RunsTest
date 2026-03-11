# =============================================================================
# Cálculo del recorte óptimo k* para el estadístico Ck
#
# El estadístico Ck(k) aplica el test condicional de corridas de Baklizi (2003)
# solo sobre las k observaciones con mayor valor absoluto, excluyendo las n-k
# más cercanas a la mediana especificada.
#
# CRITERIO: k* = max( ceil(c / (1-q)),  ceil(gamma * n) )
#
#   Componente 1 — ESTABILIDAD: ceil(c / (1-q))
#     Bajo el peor escenario de mediana mal especificada (m0 = Q_q),
#     P(X > m0) = 1-q y n1k ~ Bin(k, 1-q).
#     Se requiere E[n1k] = k*(1-q) >= c para que la distribución condicional
#     de corridas tenga soporte suficiente (con c=5, q=0.70 → cota = 17).
#
#   Componente 2 — POTENCIA: ceil(gamma * n)
#     Se retiene al menos la fracción gamma de la muestra.
#     La simulación con n in {20,30,50,100} muestra que gamma=0.60 reproduce
#     el k* empíricamente óptimo para n >= 30:
#       n=30  → k*=18 (60%n), n=50 → k*=30 (60%n), n=100 → k*=60 (60%n)
#     Para n < 30 domina la cota de estabilidad (k*=17).
#
# Interpretación del 60%:
#   Se descartan las n-k observaciones más CERCANAS a la mediana especificada,
#   que son las más afectadas por un error de especificación.
#   Las k más EXTREMAS (en valor absoluto) son las más informativas sobre
#   la forma de la distribución y, por tanto, sobre la simetría.
# =============================================================================


#' Recorte óptimo k* para el estadístico Ck
#'
#' @param n      Tamaño de muestra
#' @param q      Percentil de desplazamiento más exigente (default 0.70)
#' @param c      Mínimo de E[n1k] requerido para estabilidad (default 5)
#' @param gamma  Fracción mínima de la muestra a retener (default 0.60)
#' @param verbose Si TRUE imprime diagnóstico detallado
#'
#' @return k* (entero)
#' @references
#' Baklizi, A. (2003). A conditional distribution free runs test for symmetry. Journal of Nonparametric Statistics, 15(6), 713-718. 
#' @examples
#' cat("k* para n=30:", optimal_k(30), "\n")
#' cat("k* para n=50:", optimal_k(50), "\n")
#' 
#' # Diagnóstico detallado para los tamaños del estudio
#' optimal_k(30, verbose = TRUE)
#' optimal_k(50, verbose = TRUE)
#' optimal_k(100, verbose = TRUE)
#' 
#' # Tabla para un rango amplio de tamaños
#' cat("=== Tabla de recortes óptimos (q=0.70, c=5, gamma=0.60) ===\n")
#' print(optimal_k_table(c(20, 25, 30, 40, 50, 75, 100, 150, 200)), row.names = FALSE)

optimal_k <- function(n, q = 0.70, c = 5, gamma = 0.60, verbose = FALSE) {
  
  p_above <- 1 - q                                 # P(X > m0) bajo H0 incorrecta
  k_stab  <- ceiling(c / p_above)                  # cota de estabilidad
  k_power <- ceiling(gamma * n)                    # cota de potencia
  k_star  <- min(max(k_stab, k_power), n)          # criterio combinado
  
  if (verbose) {
    cat(sprintf("\n=== Recorte óptimo k* para n = %d ===\n", n))
    cat(sprintf("  q     = %.2f  (percentil de desplazamiento más exigente)\n", q))
    cat(sprintf("  p     = %.2f  (P(X > m0) bajo H0 incorrecta)\n", p_above))
    cat(sprintf("  c     = %.0f   (mínimo E[n1k] requerido)\n", c))
    cat(sprintf("  gamma = %.2f  (fracción mínima de la muestra)\n", gamma))
    cat(sprintf("  k_estabilidad = ceil(%.0f / %.2f) = %d\n", c, p_above, k_stab))
    cat(sprintf("  k_potencia    = ceil(%.2f * %d) = %d\n", gamma, n, k_power))
    cat(sprintf("  k* = max(%d, %d) = %d  (%.0f%% de n)\n\n",
                k_stab, k_power, k_star, 100 * k_star / n))
    
    cat(sprintf("%5s  %6s  %8s  %12s  %12s\n",
                "k", "k/n", "E[n1k]", "P(n1k<=2)", "P(n1k<=3)"))
    cat(strrep("-", 52), "\n")
    k_lo <- max(1, k_star - 5)
    k_hi <- min(n, k_star + 5)
    for (k in k_lo:k_hi) {
      marker <- if (k == k_star) " <-- k*" else ""
      cat(sprintf("%5d  %5.0f%%  %8.2f  %12.4f  %12.4f%s\n",
                  k, 100*k/n, k*p_above,
                  pbinom(2, k, p_above),
                  pbinom(3, k, p_above),
                  marker))
    }
    cat("\n")
  }
  
  return(k_star)
}


#' Tabla de k* para un vector de tamaños de muestra
#'
#' @param n_vec  Vector de tamaños de muestra
#' @param q      Percentil de desplazamiento (default 0.70)
#' @param c      Mínimo E[n1k] (default 5)
#' @param gamma  Fracción mínima (default 0.60)
#'
#' @return data.frame con n, k*, k*/n, E[n1k*], P(n1k*<=2)
optimal_k_table <- function(n_vec, q = 0.70, c = 5, gamma = 0.60) {
  p_above <- 1 - q
  results <- lapply(n_vec, function(n) {
    k <- optimal_k(n, q=q, c=c, gamma=gamma, verbose=FALSE)
    data.frame(
      n          = n,
      k_star     = k,
      prop       = paste0(round(100 * k / n), "%"),
      E_n1k      = round(k * p_above, 2),
      P_n1k_le2  = round(pbinom(2, k, p_above), 4)
    )
  })
  do.call(rbind, results)
}



