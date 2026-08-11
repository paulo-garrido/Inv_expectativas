# ==============================================================================
# OBJETIVO 2
# Modelo de parámetros variantes en el tiempo (TVP) mediante filtro de Kalman
#
# Ecuación de observación:
#   Delta exp_inf_24m_t = alpha_t + beta_t * shock_inf_t + epsilon_t
#
# Ecuaciones de transición:
#   alpha_t = alpha_{t-1} + eta_alpha,t
#   beta_t  = beta_{t-1}  + eta_beta,t
#
# Interpretación:
#   beta_t cercano a cero  -> mayor anclaje
#   beta_t positivo y alto -> menor anclaje
#
# Archivo de entrada:
#   bd_sample_modelo(1).csv
# ==============================================================================


# 0. PAQUETES ------------------------------------------------------------------

paquetes <- c(
  "dplyr",
  "readr",
  "lubridate",
  "ggplot2",
  "tidyr",
  "tibble",
  "broom",
  "KFAS"
)

paquetes_faltantes <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(paquetes_faltantes) > 0) {
  install.packages(paquetes_faltantes, dependencies = TRUE)
}

invisible(lapply(paquetes, library, character.only = TRUE))


# 1. CONFIGURACIÓN -------------------------------------------------------------

archivo_entrada <- "bd_sample_modelo(1).csv"
carpeta_salida  <- "resultados_objetivo2_tvp"

# La investigación se ha delimitado al período 2011-2026.
fecha_inicio <- as.Date("2011-01-01")

# Nivel utilizado para las bandas de confianza.
nivel_confianza <- 0.95
z_critico <- qnorm(1 - (1 - nivel_confianza) / 2)

dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)


# 2. IMPORTACIÓN Y VALIDACIÓN --------------------------------------------------

if (!file.exists(archivo_entrada)) {
  stop(
    paste0(
      "No se encontró el archivo '", archivo_entrada, "'.\n",
      "Coloque el script y la base en la misma carpeta de trabajo,\n",
      "o modifique el objeto 'archivo_entrada'."
    )
  )
}

bd_original <- readr::read_csv(
  archivo_entrada,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", ".", "null")
)

variables_requeridas <- c(
  "fecha",
  "exp_inf_12m",
  "exp_inf_24m",
  "infl_gt"
)

variables_faltantes <- setdiff(variables_requeridas, names(bd_original))

if (length(variables_faltantes) > 0) {
  stop(
    paste(
      "Faltan las siguientes variables requeridas:",
      paste(variables_faltantes, collapse = ", ")
    )
  )
}


# 3. CONSTRUCCIÓN DE VARIABLES -------------------------------------------------

bd_modelo <- bd_original %>%
  mutate(
    fecha = lubridate::ymd(fecha)
  ) %>%
  arrange(fecha) %>%
  mutate(
    # Cambio mensual de las expectativas de inflación a 24 meses.
    d_exp24 = exp_inf_24m - lag(exp_inf_24m),

    # Expectativa a 12 meses disponible en t-1.
    exp12_lag1 = lag(exp_inf_12m),

    # Proxy propuesta para la sorpresa o brecha inflacionaria.
    shock_inf = infl_gt - exp12_lag1
  ) %>%
  filter(fecha >= fecha_inicio) %>%
  filter(
    is.finite(d_exp24),
    is.finite(shock_inf)
  )

if (nrow(bd_modelo) < 36) {
  stop(
    "La muestra efectiva tiene menos de 36 observaciones. ",
    "No es suficiente para estimar el modelo propuesto."
  )
}

if (sd(bd_modelo$shock_inf, na.rm = TRUE) == 0) {
  stop("La variable shock_inf no presenta variación en la muestra.")
}


# 4. MODELO DE REFERENCIA CON PARÁMETROS CONSTANTES ---------------------------

modelo_constante <- lm(
  d_exp24 ~ shock_inf,
  data = bd_modelo
)

coef_constante <- broom::tidy(modelo_constante, conf.int = TRUE) %>%
  mutate(modelo = "Parámetros constantes")

readr::write_csv(
  coef_constante,
  file.path(carpeta_salida, "01_modelo_constante.csv")
)


# 5. FUNCIÓN PARA CONSTRUIR EL MODELO TVP -------------------------------------

construir_modelo_tvp <- function(y, x, valores_iniciales_estado) {

  n <- length(y)

  # Matriz Z_t = [1, shock_inf_t].
  # Sus elementos cambian en cada período debido a x_t.
  Z_t <- array(
    NA_real_,
    dim = c(1, 2, n)
  )

  Z_t[1, 1, ] <- 1
  Z_t[1, 2, ] <- x

  # Los estados son:
  # estado 1 = alpha_t
  # estado 2 = beta_t
  modelo <- KFAS::SSModel(
    y ~ -1 +
      KFAS::SSMcustom(
        Z = Z_t,
        T = diag(2),
        R = diag(2),
        Q = diag(NA_real_, 2),
        n = n,
        a1 = valores_iniciales_estado,
        P1 = matrix(0, nrow = 2, ncol = 2),
        P1inf = diag(2)
      ),
    H = matrix(NA_real_, nrow = 1, ncol = 1)
  )

  return(modelo)
}


# 6. FUNCIÓN DE ACTUALIZACIÓN DE VARIANZAS ------------------------------------

# par[1] = log(H)
# par[2] = log(Q_alpha)
# par[3] = log(Q_beta)
#
# La transformación exponencial garantiza que las varianzas sean positivas.

actualizar_varianzas <- function(par, modelo) {

  modelo$H[1, 1, 1] <- exp(par[1])
  modelo$Q[1, 1, 1] <- exp(par[2])
  modelo$Q[2, 2, 1] <- exp(par[3])

  return(modelo)
}


# 7. ESTIMACIÓN CON MÚLTIPLES VALORES INICIALES -------------------------------

y <- bd_modelo$d_exp24
x <- bd_modelo$shock_inf

coef_ols <- coef(modelo_constante)

estado_inicial <- c(
  alpha = unname(coef_ols["(Intercept)"]),
  beta  = unname(coef_ols["shock_inf"])
)

modelo_tvp_sin_estimar <- construir_modelo_tvp(
  y = y,
  x = x,
  valores_iniciales_estado = estado_inicial
)

var_residuo_ols <- max(
  var(residuals(modelo_constante), na.rm = TRUE),
  1e-8
)

var_y <- max(
  var(y, na.rm = TRUE),
  1e-8
)

# Se utilizan varios puntos de partida para reducir el riesgo de que el
# optimizador converja a una solución local dependiente de un único inicio.
lista_inicios <- list(
  log(c(var_residuo_ols, var_y * 1e-2, var_y * 1e-3)),
  log(c(var_residuo_ols, var_y * 1e-3, var_y * 1e-4)),
  log(c(var_residuo_ols, var_y * 1e-4, var_y * 1e-5)),
  log(c(var_residuo_ols * 0.5, var_y * 1e-3, var_y * 1e-3)),
  log(c(var_residuo_ols * 2.0, var_y * 1e-2, var_y * 1e-4))
)

ajustes <- lapply(
  seq_along(lista_inicios),
  function(i) {

    intento <- tryCatch(
      KFAS::fitSSM(
        model = modelo_tvp_sin_estimar,
        inits = lista_inicios[[i]],
        updatefn = actualizar_varianzas,
        method = "L-BFGS-B",
        lower = rep(log(1e-10), 3),
        upper = rep(log(100), 3),
        control = list(
          maxit = 5000,
          factr = 1e7
        )
      ),
      error = function(e) NULL
    )

    if (is.null(intento)) {
      return(NULL)
    }

    log_verosimilitud <- tryCatch(
      as.numeric(logLik(intento$model)),
      error = function(e) NA_real_
    )

    list(
      ajuste = intento,
      inicio = i,
      convergencia = intento$optim.out$convergence,
      mensaje = intento$optim.out$message,
      logLik = log_verosimilitud
    )
  }
)

ajustes <- Filter(Negate(is.null), ajustes)

if (length(ajustes) == 0) {
  stop("Ninguno de los intentos de estimación pudo ejecutarse.")
}

tabla_intentos <- tibble::tibble(
  intento = vapply(ajustes, function(z) z$inicio, numeric(1)),
  convergencia = vapply(ajustes, function(z) z$convergencia, numeric(1)),
  logLik = vapply(ajustes, function(z) z$logLik, numeric(1)),
  mensaje = vapply(
    ajustes,
    function(z) {
      if (is.null(z$mensaje)) "" else as.character(z$mensaje)
    },
    character(1)
  )
) %>%
  arrange(desc(logLik))

readr::write_csv(
  tabla_intentos,
  file.path(carpeta_salida, "02_intentos_estimacion.csv")
)

# Se selecciona la estimación convergente con mayor log-verosimilitud.
ajustes_convergentes <- ajustes[
  vapply(
    ajustes,
    function(z) z$convergencia == 0 && is.finite(z$logLik),
    logical(1)
  )
]

if (length(ajustes_convergentes) == 0) {
  warning(
    "Ningún intento reportó convergencia igual a cero. ",
    "Se utilizará el ajuste con mayor log-verosimilitud, ",
    "pero debe revisarse la tabla de intentos."
  )
  ajustes_elegibles <- ajustes
} else {
  ajustes_elegibles <- ajustes_convergentes
}

indice_mejor <- which.max(
  vapply(ajustes_elegibles, function(z) z$logLik, numeric(1))
)

mejor_resultado <- ajustes_elegibles[[indice_mejor]]
ajuste_tvp <- mejor_resultado$ajuste
modelo_tvp <- ajuste_tvp$model


# 8. FILTRO Y SUAVIZADOR DE KALMAN --------------------------------------------

resultado_kalman <- KFAS::KFS(
  modelo_tvp,
  filtering = c("state", "mean"),
  smoothing = c("state", "mean", "disturbance")
)

# Estados suavizados.
alpha_t <- as.numeric(resultado_kalman$alphahat[, 1])
beta_t  <- as.numeric(resultado_kalman$alphahat[, 2])

# Errores estándar de los estados suavizados.
se_alpha_t <- sqrt(pmax(as.numeric(resultado_kalman$V[1, 1, ]), 0))
se_beta_t  <- sqrt(pmax(as.numeric(resultado_kalman$V[2, 2, ]), 0))

# Valor ajustado por el modelo.
ajustado_tvp <- alpha_t + beta_t * x

resultados_tvp <- bd_modelo %>%
  transmute(
    fecha,
    exp_inf_12m,
    exp_inf_24m,
    infl_gt,
    d_exp24,
    exp12_lag1,
    shock_inf,

    alpha_t = alpha_t,
    se_alpha_t = se_alpha_t,
    alpha_li_95 = alpha_t - z_critico * se_alpha_t,
    alpha_ls_95 = alpha_t + z_critico * se_alpha_t,

    beta_t = beta_t,
    se_beta_t = se_beta_t,
    beta_li_95 = beta_t - z_critico * se_beta_t,
    beta_ls_95 = beta_t + z_critico * se_beta_t,

    # Un valor mayor implica más anclaje únicamente bajo la transformación
    # propuesta A_t = -beta_t.
    indice_anclaje = -beta_t,

    ajustado_tvp = ajustado_tvp,
    residuo_tvp = d_exp24 - ajustado_tvp
  )

readr::write_csv(
  resultados_tvp,
  file.path(carpeta_salida, "03_resultados_tvp_kalman.csv")
)


# 9. PARÁMETROS DE VARIANZA ----------------------------------------------------

H_estimado       <- as.numeric(modelo_tvp$H[1, 1, 1])
Q_alpha_estimado <- as.numeric(modelo_tvp$Q[1, 1, 1])
Q_beta_estimado  <- as.numeric(modelo_tvp$Q[2, 2, 1])

tabla_varianzas <- tibble::tibble(
  parametro = c(
    "H",
    "Q_alpha",
    "Q_beta",
    "Q_alpha / H",
    "Q_beta / H"
  ),
  estimacion = c(
    H_estimado,
    Q_alpha_estimado,
    Q_beta_estimado,
    Q_alpha_estimado / H_estimado,
    Q_beta_estimado / H_estimado
  ),
  interpretacion = c(
    "Varianza del error de observación",
    "Varianza del cambio mensual del intercepto",
    "Varianza del cambio mensual de la sensibilidad beta_t",
    "Cambio del intercepto relativo al ruido de observación",
    "Cambio de beta_t relativo al ruido de observación"
  )
)

readr::write_csv(
  tabla_varianzas,
  file.path(carpeta_salida, "04_varianzas_estimadas.csv")
)


# 10. DIAGNÓSTICOS -------------------------------------------------------------

# Residuos recursivos estandarizados. Si la versión instalada de KFAS no
# permite generarlos mediante rstandard(), se conserva NA y el resto del
# script continúa ejecutándose.
residuos_estandarizados <- tryCatch(
  as.numeric(stats::rstandard(resultado_kalman, type = "recursive")),
  error = function(e) rep(NA_real_, nrow(resultados_tvp))
)

resultados_tvp$residuo_estandarizado <- residuos_estandarizados

readr::write_csv(
  resultados_tvp,
  file.path(carpeta_salida, "03_resultados_tvp_kalman.csv")
)

residuos_validos <- residuos_estandarizados[
  is.finite(residuos_estandarizados)
]

if (length(residuos_validos) >= 20) {

  rezago_lb <- min(12, floor(length(residuos_validos) / 5))

  prueba_ljung_box <- Box.test(
    residuos_validos,
    lag = rezago_lb,
    type = "Ljung-Box"
  )

  media_r <- mean(residuos_validos)
  desv_r  <- sd(residuos_validos)

  asimetria <- mean(
    ((residuos_validos - media_r) / desv_r)^3
  )

  curtosis <- mean(
    ((residuos_validos - media_r) / desv_r)^4
  )

  estadistico_jb <- length(residuos_validos) / 6 *
    (
      asimetria^2 +
        ((curtosis - 3)^2) / 4
    )

  p_valor_jb <- pchisq(
    estadistico_jb,
    df = 2,
    lower.tail = FALSE
  )

  tabla_diagnosticos <- tibble::tibble(
    diagnostico = c(
      "Ljung-Box",
      "Jarque-Bera aproximado"
    ),
    estadistico = c(
      unname(prueba_ljung_box$statistic),
      estadistico_jb
    ),
    grados_libertad = c(
      unname(prueba_ljung_box$parameter),
      2
    ),
    p_valor = c(
      prueba_ljung_box$p.value,
      p_valor_jb
    )
  )

} else {

  tabla_diagnosticos <- tibble::tibble(
    diagnostico = "No disponible",
    estadistico = NA_real_,
    grados_libertad = NA_real_,
    p_valor = NA_real_
  )
}

readr::write_csv(
  tabla_diagnosticos,
  file.path(carpeta_salida, "05_diagnosticos_residuales.csv")
)


# 11. RESUMEN DE LA TRAYECTORIA DE BETA ---------------------------------------

resumen_beta <- resultados_tvp %>%
  summarise(
    fecha_inicial = min(fecha),
    fecha_final = max(fecha),
    observaciones = n(),
    beta_inicial = first(beta_t),
    beta_final = last(beta_t),
    cambio_beta = beta_final - beta_inicial,
    beta_promedio = mean(beta_t),
    beta_minimo = min(beta_t),
    fecha_beta_minimo = fecha[which.min(beta_t)],
    beta_maximo = max(beta_t),
    fecha_beta_maximo = fecha[which.max(beta_t)],
    proporcion_beta_significativo_positivo = mean(beta_li_95 > 0),
    proporcion_beta_no_distinto_cero = mean(
      beta_li_95 <= 0 & beta_ls_95 >= 0
    )
  )

readr::write_csv(
  resumen_beta,
  file.path(carpeta_salida, "06_resumen_trayectoria_beta.csv")
)


# 12. GRÁFICAS -----------------------------------------------------------------

grafica_beta <- ggplot(
  resultados_tvp,
  aes(x = fecha, y = beta_t)
) +
  geom_ribbon(
    aes(ymin = beta_li_95, ymax = beta_ls_95),
    alpha = 0.20
  ) +
  geom_line(linewidth = 0.8) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Sensibilidad variable de las expectativas a 24 meses",
    subtitle = paste0(
      "Modelo TVP–Kalman; bandas de confianza del ",
      nivel_confianza * 100,
      "%"
    ),
    x = NULL,
    y = expression(beta[t]),
    caption = paste0(
      "Ecuación: Delta e[t]^24 == alpha[t] + beta[t] * shock[t] + epsilon[t]"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(carpeta_salida, "07_trayectoria_beta_t.png"),
  plot = grafica_beta,
  width = 10,
  height = 6,
  dpi = 300
)

grafica_anclaje <- ggplot(
  resultados_tvp,
  aes(x = fecha, y = indice_anclaje)
) +
  geom_line(linewidth = 0.8) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Índice complementario de anclaje",
    subtitle = expression(A[t] == -beta[t]),
    x = NULL,
    y = expression(A[t])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(carpeta_salida, "08_indice_anclaje.png"),
  plot = grafica_anclaje,
  width = 10,
  height = 6,
  dpi = 300
)

grafica_ajuste <- resultados_tvp %>%
  select(fecha, observado = d_exp24, ajustado = ajustado_tvp) %>%
  pivot_longer(
    cols = c(observado, ajustado),
    names_to = "serie",
    values_to = "valor"
  ) %>%
  ggplot(
    aes(x = fecha, y = valor, linetype = serie)
  ) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Cambio observado y ajustado de las expectativas a 24 meses",
    x = NULL,
    y = expression(Delta * e[t]^24),
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(carpeta_salida, "09_ajuste_modelo_tvp.png"),
  plot = grafica_ajuste,
  width = 10,
  height = 6,
  dpi = 300
)

if (length(residuos_validos) >= 20) {

  datos_residuos <- tibble::tibble(
    fecha = resultados_tvp$fecha[
      is.finite(resultados_tvp$residuo_estandarizado)
    ],
    residuo_estandarizado = residuos_validos
  )

  grafica_residuos <- ggplot(
    datos_residuos,
    aes(x = fecha, y = residuo_estandarizado)
  ) +
    geom_line(linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = c(-2, 2), linetype = "dotted") +
    labs(
      title = "Residuos recursivos estandarizados",
      x = NULL,
      y = "Residuo estandarizado"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title.position = "plot",
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(carpeta_salida, "10_residuos_estandarizados.png"),
    plot = grafica_residuos,
    width = 10,
    height = 5.5,
    dpi = 300
  )
}


# 13. INFORMACIÓN FINAL EN CONSOLA ---------------------------------------------

cat("\n")
cat("============================================================\n")
cat("MODELO TVP–KALMAN FINALIZADO\n")
cat("============================================================\n")
cat("Archivo de entrada: ", archivo_entrada, "\n", sep = "")
cat("Muestra efectiva: ", format(min(bd_modelo$fecha), "%Y-%m"),
    " a ", format(max(bd_modelo$fecha), "%Y-%m"), "\n", sep = "")
cat("Número de observaciones: ", nrow(bd_modelo), "\n", sep = "")
cat("Intento seleccionado: ", mejor_resultado$inicio, "\n", sep = "")
cat("Convergencia del optimizador: ",
    mejor_resultado$convergencia, "\n", sep = "")
cat("Log-verosimilitud: ",
    round(mejor_resultado$logLik, 4), "\n", sep = "")
cat("\nVarianzas estimadas:\n")
cat("H       = ", format(H_estimado, scientific = TRUE), "\n", sep = "")
cat("Q_alpha = ", format(Q_alpha_estimado, scientific = TRUE), "\n", sep = "")
cat("Q_beta  = ", format(Q_beta_estimado, scientific = TRUE), "\n", sep = "")
cat("Q_beta/H= ",
    format(Q_beta_estimado / H_estimado, scientific = TRUE), "\n", sep = "")
cat("\nTrayectoria de beta_t:\n")
cat("beta inicial = ", round(first(resultados_tvp$beta_t), 5), "\n", sep = "")
cat("beta final   = ", round(last(resultados_tvp$beta_t), 5), "\n", sep = "")
cat("cambio total = ",
    round(last(resultados_tvp$beta_t) - first(resultados_tvp$beta_t), 5),
    "\n", sep = "")
cat("\nResultados guardados en: ", carpeta_salida, "\n", sep = "")
cat("============================================================\n")
