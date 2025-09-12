# 1. La función de densidad 
dULL <- function(x, mu, sigma, log = FALSE) {
  lambda <- mu
  if (any(sigma <= 0)) stop("El parámetro sigma debe ser positivo.")
  if (any(lambda < 0)) stop("El parámetro mu (lambda) no puede ser negativo.")
  if (any(x <= 0 | x >= 1)) stop("x debe estar en el intervalo (0, 1).")
  
  log_pdf <- 2 * log(sigma) - log(1 + lambda * sigma) + 
    log(lambda - log(x)) + (sigma - 1) * log(x)
  
  if (log) {
    return(log_pdf)
  } else {
    result <- exp(log_pdf)
    result[is.nan(result) | is.infinite(result)] <- 0
    return(result)
  }
}

# 2. La lista de parámetros
# Primero es SIGMA luego el segundo es MU
parametros <- list(
  p1 = c(sigma = 0.5, mu = 0),      # Etiqueta (0.5, 0)
  p2 = c(sigma = 0.5, mu = 10),     # Etiqueta (0.5, 10)
  p3 = c(sigma = 2,   mu = 0.1),    # Etiqueta (2, 0.1)
  p4 = c(sigma = 5,   mu = 0.01),   # Etiqueta (5, 0.01)
  p5 = c(sigma = 20,  mu = 0.001),  # Etiqueta (20, 0.001)
  p6 = c(sigma = 20,  mu = 5)       # Etiqueta (20, 5)
)

# 3. Grafica
colores <- c("black", "gray50", "green4", "blue", "orange", "red")
leyendas <- c("(0.5, 0)", "(0.5, 10)", "(2, 0.1)", "(5, 0.01)", "(20, 0.001)", "(20, 5)")
x_vals <- seq(0.001, 0.999, length.out = 1000)

y_densities <- lapply(parametros, function(p) {
  dULL(x_vals, mu = p["mu"], sigma = p["sigma"])
})

max_y <- 11 # Fijamos el eje Y como en la gráfica para una mejor comparación
plot(NA, xlim = c(0, 1), ylim = c(0, max_y),
     xlab = "x", ylab = "f(x)",
     main = "Réplica Final de la Gráfica Log-Lindley")
grid()

for (i in 1:length(parametros)) {
  lines(x_vals, y_densities[[i]], col = colores[i], lwd = 1.5)
}

legend("topleft", legend = leyendas, col = colores, lwd = 2, bty = "n", 
       title = expression(paste("Valores de (", sigma, ", ", mu, ")")))

