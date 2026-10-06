# =============================================================================
# Modelo de Solow-Swan: dinámica de la ecuación fundamental
# Basado en Sala-i-Martin, "Apuntes de crecimiento económico", capítulo 1
#
#   k_dot = s*A*k^alpha - (delta + n + x)*k          [1.15] y [1.24]
#   gamma_k = s*A*k^(alpha-1) - (delta + n + x)       [1.20']
#   k*  = (s*A / (delta + n + x))^(1/(1-alpha))       [1.16]
#   beta* = (1 - alpha)*(delta + n + x)               [1.27]-[1.28]
# =============================================================================

library(ggplot2)
library(tidyr)
library(dplyr)

# ---- 1. Funciones del modelo ------------------------------------------------
f      <- function(k, p) p$A * k^p$alpha
k_dot  <- function(k, p) p$s * f(k, p) - (p$delta + p$n + p$x) * k
k_star <- function(p) (p$s * p$A / (p$delta + p$n + p$x))^(1 / (1 - p$alpha))
k_oro  <- function(p) (p$alpha * p$A / (p$delta + p$n + p$x))^(1 / (1 - p$alpha))

# ---- 2. Simulación con Runge-Kutta de 4o orden ------------------------------
# choque: lista(var = "s", valor = 0.30, T = 10)  o  NULL para no aplicar choque
simular_solow <- function(p, k0_rel = 1, H = 80, dt = 0.05, choque = NULL) {
  pasos <- round(H / dt)
  k <- k0_rel * k_star(p)
  res <- vector("list", pasos + 1)
  for (i in 0:pasos) {
    t <- i * dt
    q <- p
    if (!is.null(choque) && t >= choque$T) q[[choque$var]] <- choque$valor
    y <- f(k, q)
    res[[i + 1]] <- data.frame(
      t = t, k = k, y = y, c = (1 - q$s) * y,
      inversion = q$s * y,            # s*A*k^alpha
      depreciacion = q$delta * k,     # delta*k
      dilucion = q$n * k,             # n*k
      tecnologia = q$x * k,           # x*k
      k_dot = k_dot(k, q),
      gamma_k = q$s * q$A * k^(q$alpha - 1) - (q$delta + q$n + q$x),
      k_estrella = k_star(q)
    )
    k1 <- k_dot(k, q); k2 <- k_dot(k + dt / 2 * k1, q)
    k3 <- k_dot(k + dt / 2 * k2, q); k4 <- k_dot(k + dt * k3, q)
    k <- k + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
  }
  out <- bind_rows(res)
  out$gamma_y <- p$alpha * out$gamma_k   # gamma_y = alpha * gamma_k
  out
}

# ---- 3. Parámetros del libro (sección 1.6) ----------------------------------
p <- list(s = 0.20, A = 1, alpha = 0.30, delta = 0.10, n = 0.01, x = 0)

cat("k*          =", round(k_star(p), 4), "\n")
cat("k oro       =", round(k_oro(p), 4), "  (s oro = alpha =", p$alpha, ")\n")
beta <- (1 - p$alpha) * (p$delta + p$n + p$x)
cat("beta*       =", round(beta, 4), "\n")
cat("Vida media  =", round(log(2) / beta, 1), "años\n")

# ---- 4. Escenario: aumento de s de 0.20 a 0.30 en t = 10 (gráfico 1.11) -----
sim <- simular_solow(p, k0_rel = 1, H = 80, choque = list(var = "s", valor = 0.30, T = 10))

# 4a. Descomposición de la ecuación fundamental en el tiempo
sim |>
  select(t, inversion, depreciacion, dilucion, k_dot) |>
  pivot_longer(-t, names_to = "termino", values_to = "valor") |>
  ggplot(aes(t, valor, colour = termino)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 10, linetype = "dotted") +
  labs(title = "Términos de la ecuación fundamental de Solow-Swan",
       subtitle = "Aumento permanente de s en t = 10",
       x = "Tiempo (años)", y = NULL, colour = NULL) +
  theme_minimal()

# 4b. Trayectorias de k, y, c
sim |>
  select(t, k, y, c) |>
  pivot_longer(-t) |>
  ggplot(aes(t, value, colour = name)) +
  geom_line(linewidth = 0.9) +
  labs(title = "Trayectorias de transición", x = "Tiempo (años)", y = NULL, colour = NULL) +
  theme_minimal()

# 4c. Tasas de crecimiento
sim |>
  select(t, gamma_k, gamma_y) |>
  pivot_longer(-t) |>
  ggplot(aes(t, value, colour = name)) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Tasas de crecimiento: gamma_y = alpha * gamma_k", x = "Tiempo (años)",
       y = NULL, colour = NULL) +
  theme_minimal()

# ---- 5. Diagrama de Solow (gráfico 1.1) -------------------------------------
kk <- seq(0.01, 1.5 * k_star(p), length.out = 300)
data.frame(k = kk, `f(k)` = f(kk, p), `s f(k)` = p$s * f(kk, p),
           `(delta+n)k` = (p$delta + p$n + p$x) * kk, check.names = FALSE) |>
  pivot_longer(-k) |>
  ggplot(aes(k, value, colour = name)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = k_star(p), linetype = "dashed") +
  labs(title = "Estado estacionario en el modelo de Solow-Swan", y = NULL, colour = NULL) +
  theme_minimal()
