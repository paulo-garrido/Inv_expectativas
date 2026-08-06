# ==============================================================================
# OBJETIVO 2 - SEGUNDA ESPECIFICACIÓN
# Modelo TVP-Kalman con persistencia AR(1)
#
# Ecuación de observación:
#   Delta exp_inf_24m_t = alpha + rho * Delta exp_inf_24m_(t-1)
#                         + beta_t * shock_inf_t + epsilon_t
#
# Ecuación de transición:
#   beta_t = beta_(t-1) + eta_beta,t
#
# Parámetros constantes:
#   alpha y rho
#
# Parámetro variante en el tiempo:
#   beta_t
#
# Interpretación:
#   beta_t cercano a cero  -> mayor anclaje
#   beta_t positivo y alto -> menor anclaje
#
# ==============================================================================
# ------------------------------------------------------------
# 0) Limpieza de entorno
# ------------------------------------------------------------
rm(list = ls())
cat("\014")
if (!is.null(dev.list())) dev.off()


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

# El script busca primero el archivo sin sufijo y después la versión con “(1)”.
archivos_candidatos <- c(
  "bd_sample_modelo.csv",
  "bd_sample_modelo(1).csv"
)

archivo_entrada <- archivos_candidatos[
  file.exists(archivos_candidatos)
][1]

if (is.na(archivo_entrada)) {
  stop(
    paste0(
      "No se encontró ninguno de estos archivos en el directorio de trabajo:\n",
      paste0("- ", archivos_candidatos, collapse = "\n"), "\n\n",
      "Directorio actual: ", getwd()
    )
  )
}

# Se usa una carpeta diferente para conservar intactos los resultados de la
# primera especificación.
carpeta_salida <- file.path(
  getwd(),
  "resultados_objetivo2_tvp_ar1"
)

fecha_inicio <- as.Date("2011-01-01")

nivel_confianza <- 0.95
z_critico <- qnorm(1 - (1 - nivel_confianza) / 2)

dir.create(
  carpeta_salida,
  showWarnings = FALSE,
  recursive = TRUE
)

if (!dir.exists(carpeta_salida)) {
  stop("No fue posible crear la carpeta de resultados: ", carpeta_salida)
}

if (file.access(carpeta_salida, 2) != 0) {
  stop("R no tiene permiso de escritura en: ", carpeta_salida)
}


# 2. IMPORTACIÓN Y VALIDACIÓN --------------------------------------------------

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

variables_faltantes <- setdiff(
  variables_requeridas,
  names(bd_original)
)

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

    # Persistencia de la revisión mensual de expectativas.
    d_exp24_lag1 = lag(d_exp24),

    # Expectativa de inflación a 12 meses observada en t-1.
    exp12_lag1 = lag(exp_inf_12m),

    # Brecha entre inflación observada y expectativa rezagada.
    shock_inf = infl_gt - exp12_lag1
  ) %>%
  filter(fecha >= fecha_inicio) %>%
  filter(
    is.finite(d_exp24),
    is.finite(d_exp24_lag1),
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

if (sd(bd_modelo$d_exp24_lag1, na.rm = TRUE) == 0) {
  stop("La variable d_exp24_lag1 no presenta variación en la muestra.")
}


# 4. MODELO DE REFERENCIA CON PARÁMETROS CONSTANTES ---------------------------

# Benchmark con la misma estructura dinámica, pero beta constante.
modelo_constante_ar1 <- lm(
  d_exp24 ~ d_exp24_lag1 + shock_inf,
  data = bd_modelo
)

coef_constante_ar1 <- broom::tidy(
  modelo_constante_ar1,
  conf.int = TRUE,
  conf.level = nivel_confianza
) %>%
  mutate(modelo = "AR(1) con parámetros constantes")

readr::write_csv(
  coef_constante_ar1,
  file.path(carpeta_salida, "01_modelo_constante_ar1.csv")
)


# 5. CONSTRUCCIÓN DEL MODELO TVP-AR(1) ----------------------------------------

construir_modelo_tvp_ar1 <- function(
    y,
    y_lag1,
    x,
    valores_iniciales_estado
) {

  n <- length(y)

  # Vector de estados:
  #   estado 1 = alpha, constante
  #   estado 2 = rho, constante
  #   estado 3 = beta_t, variante en el tiempo
  #
  # Matriz de medición:
  #   Z_t = [1, Delta exp24_(t-1), shock_inf_t]
  Z_t <- array(
    NA_real_,
    dim = c(1, 3, n)
  )

  Z_t[1, 1, ] <- 1
  Z_t[1, 2, ] <- y_lag1
  Z_t[1, 3, ] <- x

  # Solo beta_t recibe una innovación de estado.
  # alpha y rho permanecen constantes porque las dos primeras filas de R
  # son iguales a cero.
  R_beta <- matrix(
    c(0, 0, 1),
    nrow = 3,
    ncol = 1
  )

  modelo <- KFAS::SSModel(
    y ~ -1 +
      SSMcustom(
        Z = Z_t,
        T = diag(3),
        R = R_beta,
        Q = matrix(NA_real_, nrow = 1, ncol = 1),
        n = n,
        a1 = valores_iniciales_estado,
        P1 = matrix(0, nrow = 3, ncol = 3),
        P1inf = diag(3)
      ),
    H = matrix(NA_real_, nrow = 1, ncol = 1)
  )

  return(modelo)
}


# 6. ACTUALIZACIÓN DE HIPERPARÁMETROS -----------------------------------------

# par[1] = log(H)
# par[2] = log(Q_beta)
#
# La transformación exponencial garantiza varianzas positivas.

actualizar_varianzas_ar1 <- function(par, modelo) {

  modelo$H[1, 1, 1] <- exp(par[1])
  modelo$Q[1, 1, 1] <- exp(par[2])

  return(modelo)
}


# 7. ESTIMACIÓN CON MÚLTIPLES VALORES INICIALES -------------------------------

y      <- bd_modelo$d_exp24
y_lag1 <- bd_modelo$d_exp24_lag1
x      <- bd_modelo$shock_inf

coef_ols <- coef(modelo_constante_ar1)

if (any(!is.finite(coef_ols))) {
  stop(
    "El modelo OLS inicial produjo coeficientes no finitos. ",
    "Revise colinealidad o datos faltantes."
  )
}

estado_inicial <- c(
  alpha = unname(coef_ols["(Intercept)"]),
  rho   = unname(coef_ols["d_exp24_lag1"]),
  beta  = unname(coef_ols["shock_inf"])
)

modelo_tvp_sin_estimar <- construir_modelo_tvp_ar1(
  y = y,
  y_lag1 = y_lag1,
  x = x,
  valores_iniciales_estado = estado_inicial
)

var_residuo_ols <- max(
  var(residuals(modelo_constante_ar1), na.rm = TRUE),
  1e-8
)

var_y <- max(
  var(y, na.rm = TRUE),
  1e-8
)

# Múltiples puntos de partida para H y Q_beta.
lista_inicios <- list(
  log(c(var_residuo_ols, var_y * 1e-2)),
  log(c(var_residuo_ols, var_y * 1e-3)),
  log(c(var_residuo_ols, var_y * 1e-4)),
  log(c(var_residuo_ols * 0.5, var_y * 1e-3)),
  log(c(var_residuo_ols * 2.0, var_y * 1e-3)),
  log(c(var_residuo_ols, 1e-6))
)

ajustes <- lapply(
  seq_along(lista_inicios),
  function(i) {

    intento <- tryCatch(
      KFAS::fitSSM(
        model = modelo_tvp_sin_estimar,
        inits = lista_inicios[[i]],
        updatefn = actualizar_varianzas_ar1,
        method = "L-BFGS-B",
        lower = rep(log(1e-10), 2),
        upper = rep(log(100), 2),
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
  file.path(carpeta_salida, "02_intentos_estimacion_ar1.csv")
)

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
ajuste_tvp_ar1  <- mejor_resultado$ajuste
modelo_tvp_ar1  <- ajuste_tvp_ar1$model


# 8. FILTRO Y SUAVIZADOR DE KALMAN --------------------------------------------

resultado_kalman <- KFAS::KFS(
  modelo_tvp_ar1,
  filtering = c("state", "mean"),
  smoothing = c("state", "mean", "disturbance")
)

# Estados suavizados.
alpha_t <- as.numeric(resultado_kalman$alphahat[, 1])
rho_t   <- as.numeric(resultado_kalman$alphahat[, 2])
beta_t  <- as.numeric(resultado_kalman$alphahat[, 3])

# Errores estándar de los estados suavizados.
se_alpha_t <- sqrt(
  pmax(as.numeric(resultado_kalman$V[1, 1, ]), 0)
)

se_rho_t <- sqrt(
  pmax(as.numeric(resultado_kalman$V[2, 2, ]), 0)
)

se_beta_t <- sqrt(
  pmax(as.numeric(resultado_kalman$V[3, 3, ]), 0)
)

# Como alpha y rho no tienen innovación de estado, sus estimaciones suavizadas
# deben ser constantes salvo diferencias numéricas mínimas.
alpha_estimado <- mean(alpha_t, na.rm = TRUE)
rho_estimado   <- mean(rho_t, na.rm = TRUE)

se_alpha <- mean(se_alpha_t, na.rm = TRUE)
se_rho   <- mean(se_rho_t, na.rm = TRUE)

z_alpha <- alpha_estimado / se_alpha
z_rho   <- rho_estimado / se_rho

p_alpha <- 2 * pnorm(abs(z_alpha), lower.tail = FALSE)
p_rho   <- 2 * pnorm(abs(z_rho), lower.tail = FALSE)

ajustado_tvp <- alpha_t + rho_t * y_lag1 + beta_t * x

resultados_tvp <- bd_modelo %>%
  transmute(
    fecha,
    exp_inf_12m,
    exp_inf_24m,
    infl_gt,
    d_exp24,
    d_exp24_lag1,
    exp12_lag1,
    shock_inf,

    alpha = alpha_t,
    se_alpha = se_alpha_t,
    alpha_li_95 = alpha_t - z_critico * se_alpha_t,
    alpha_ls_95 = alpha_t + z_critico * se_alpha_t,

    rho = rho_t,
    se_rho = se_rho_t,
    rho_li_95 = rho_t - z_critico * se_rho_t,
    rho_ls_95 = rho_t + z_critico * se_rho_t,

    beta_t = beta_t,
    se_beta_t = se_beta_t,
    beta_li_95 = beta_t - z_critico * se_beta_t,
    beta_ls_95 = beta_t + z_critico * se_beta_t,

    indice_anclaje = -beta_t,

    ajustado_tvp = ajustado_tvp,
    residuo_tvp = d_exp24 - ajustado_tvp
  )

readr::write_csv(
  resultados_tvp,
  file.path(carpeta_salida, "03_resultados_tvp_kalman_ar1.csv")
)


# 9. COEFICIENTES CONSTANTES Y VARIANZAS --------------------------------------

tabla_coeficientes_fijos <- tibble::tibble(
  parametro = c("alpha", "rho"),
  estimacion = c(alpha_estimado, rho_estimado),
  error_estandar = c(se_alpha, se_rho),
  limite_inferior_95 = c(
    alpha_estimado - z_critico * se_alpha,
    rho_estimado - z_critico * se_rho
  ),
  limite_superior_95 = c(
    alpha_estimado + z_critico * se_alpha,
    rho_estimado + z_critico * se_rho
  ),
  estadistico_z = c(z_alpha, z_rho),
  p_valor = c(p_alpha, p_rho),
  interpretacion = c(
    "Intercepto constante de la ecuación de observación",
    "Persistencia AR(1) de las revisiones de expectativas"
  )
)

readr::write_csv(
  tabla_coeficientes_fijos,
  file.path(carpeta_salida, "04_coeficientes_fijos_ar1.csv")
)

H_estimado      <- as.numeric(modelo_tvp_ar1$H[1, 1, 1])
Q_beta_estimado <- as.numeric(modelo_tvp_ar1$Q[1, 1, 1])

tabla_varianzas <- tibble::tibble(
  parametro = c(
    "H",
    "Q_beta",
    "Q_beta / H",
    "sqrt(H)",
    "sqrt(Q_beta)"
  ),
  estimacion = c(
    H_estimado,
    Q_beta_estimado,
    Q_beta_estimado / H_estimado,
    sqrt(H_estimado),
    sqrt(Q_beta_estimado)
  ),
  interpretacion = c(
    "Varianza del error de observación",
    "Varianza del cambio mensual de beta_t",
    "Cambio de beta_t relativo al ruido de observación",
    "Desviación estándar del error de observación",
    "Desviación estándar de la innovación mensual de beta_t"
  )
)

readr::write_csv(
  tabla_varianzas,
  file.path(carpeta_salida, "05_varianzas_estimadas_ar1.csv")
)


# 10. DIAGNÓSTICOS RESIDUALES --------------------------------------------------

residuos_estandarizados <- tryCatch(
  as.numeric(stats::rstandard(resultado_kalman, type = "recursive")),
  error = function(e) rep(NA_real_, nrow(resultados_tvp))
)

resultados_tvp$residuo_estandarizado <- residuos_estandarizados

readr::write_csv(
  resultados_tvp,
  file.path(carpeta_salida, "03_resultados_tvp_kalman_ar1.csv")
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

  if (!is.finite(desv_r) || desv_r <= 0) {
    asimetria <- NA_real_
    curtosis <- NA_real_
    estadistico_jb <- NA_real_
    p_valor_jb <- NA_real_
  } else {
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
  }

  tabla_diagnosticos <- tibble::tibble(
    diagnostico = c(
      paste0("Ljung-Box, rezago ", rezago_lb),
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
    ),
    conclusion_5pct = c(
      ifelse(
        prueba_ljung_box$p.value < 0.05,
        "Se rechaza ausencia de autocorrelación",
        "No se rechaza ausencia de autocorrelación"
      ),
      ifelse(
        p_valor_jb < 0.05,
        "Se rechaza normalidad",
        "No se rechaza normalidad"
      )
    )
  )

} else {

  tabla_diagnosticos <- tibble::tibble(
    diagnostico = "No disponible",
    estadistico = NA_real_,
    grados_libertad = NA_real_,
    p_valor = NA_real_,
    conclusion_5pct = "No disponible"
  )
}

readr::write_csv(
  tabla_diagnosticos,
  file.path(carpeta_salida, "06_diagnosticos_residuales_ar1.csv")
)


# 11. RESUMEN DE LA TRAYECTORIA DE BETA ---------------------------------------

resumen_beta <- resultados_tvp %>%
  summarise(
    fecha_inicial = min(fecha),
    fecha_final = max(fecha),
    observaciones = n(),
    alpha_estimado = alpha_estimado,
    rho_estimado = rho_estimado,
    rho_estable = abs(rho_estimado) < 1,
    beta_inicial = first(beta_t),
    beta_final = last(beta_t),
    cambio_beta = beta_final - beta_inicial,
    cambio_porcentual_beta = 100 *
      (beta_final - beta_inicial) / abs(beta_inicial),
    beta_promedio = mean(beta_t),
    beta_minimo = min(beta_t),
    fecha_beta_minimo = fecha[which.min(beta_t)],
    beta_maximo = max(beta_t),
    fecha_beta_maximo = fecha[which.max(beta_t)],
    proporcion_beta_significativo_positivo = mean(beta_li_95 > 0),
    proporcion_beta_significativo_negativo = mean(beta_ls_95 < 0),
    proporcion_beta_no_distinto_cero = mean(
      beta_li_95 <= 0 & beta_ls_95 >= 0
    )
  )

readr::write_csv(
  resumen_beta,
  file.path(carpeta_salida, "07_resumen_trayectoria_beta_ar1.csv")
)


# 12. COMPARACIÓN INFORMATIVA CON EL BENCHMARK --------------------------------

loglik_ols <- as.numeric(logLik(modelo_constante_ar1))
loglik_tvp <- as.numeric(logLik(modelo_tvp_ar1))

# Conteo aproximado de parámetros para criterios de información:
# OLS: alpha, rho, beta y varianza residual = 4.
# TVP: alpha, rho, beta inicial, H y Q_beta = 5.
k_ols <- 4
k_tvp <- 5
n_obs <- nrow(bd_modelo)

tabla_comparacion <- tibble::tibble(
  modelo = c(
    "AR(1) con beta constante",
    "AR(1) con beta_t variante"
  ),
  log_verosimilitud = c(loglik_ols, loglik_tvp),
  parametros_aproximados = c(k_ols, k_tvp),
  AIC_aproximado = c(
    -2 * loglik_ols + 2 * k_ols,
    -2 * loglik_tvp + 2 * k_tvp
  ),
  BIC_aproximado = c(
    -2 * loglik_ols + log(n_obs) * k_ols,
    -2 * loglik_tvp + log(n_obs) * k_tvp
  ),
  nota = c(
    "Modelo lineal de referencia",
    "Criterios calculados con conteo aproximado de hiperparámetros y estados iniciales"
  )
)

readr::write_csv(
  tabla_comparacion,
  file.path(carpeta_salida, "08_comparacion_modelos_ar1.csv")
)


# 13. GRÁFICAS -----------------------------------------------------------------

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
      "Modelo TVP-Kalman con componente AR(1); bandas del ",
      nivel_confianza * 100,
      "%"
    ),
    x = NULL,
    y = expression(beta[t]),
    caption = paste0(
      "Ecuación: Delta e[t]^24 = alpha + rho Delta e[t-1]^24 + ",
      "beta[t] shock[t] + epsilon[t]"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(carpeta_salida, "09_trayectoria_beta_t_ar1.png"),
  plot = grafica_beta,
  width = 10,
  height = 6,
  dpi = 300
)

grafica_beta
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
  filename = file.path(carpeta_salida, "10_indice_anclaje_ar1.png"),
  plot = grafica_anclaje,
  width = 10,
  height = 6,
  dpi = 300
)


grafica_ajuste <- resultados_tvp %>%
  dplyr::select(
    fecha,
    observado = d_exp24,
    ajustado = ajustado_tvp
  ) %>%
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
    subtitle = "Modelo TVP-Kalman con persistencia AR(1)",
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
  filename = file.path(carpeta_salida, "11_ajuste_modelo_tvp_ar1.png"),
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
      subtitle = "Modelo TVP-Kalman con persistencia AR(1)",
      x = NULL,
      y = "Residuo estandarizado"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title.position = "plot",
      panel.grid.minor = element_blank()
    )
  grafica_residuos
  ggsave(
    filename = file.path(
      carpeta_salida,
      "12_residuos_estandarizados_ar1.png"
    ),
    plot = grafica_residuos,
    width = 10,
    height = 5.5,
    dpi = 300
  )

  png(
    filename = file.path(carpeta_salida, "13_acf_residuos_ar1.png"),
    width = 1800,
    height = 1100,
    res = 180
  )
  acf(
    residuos_validos,
    lag.max = 24,
    main = "ACF de residuos recursivos estandarizados"
  )
  dev.off()
}


# 14. INFORMACIÓN FINAL EN CONSOLA --------------------------------------------

cat("\n")
cat("============================================================\n")
cat("MODELO TVP-KALMAN AR(1) FINALIZADO\n")
cat("============================================================\n")
cat("Archivo de entrada: ", archivo_entrada, "\n", sep = "")
cat(
  "Muestra efectiva: ",
  format(min(bd_modelo$fecha), "%Y-%m"),
  " a ",
  format(max(bd_modelo$fecha), "%Y-%m"),
  "\n",
  sep = ""
)
cat("Número de observaciones: ", nrow(bd_modelo), "\n", sep = "")
cat("Intento seleccionado: ", mejor_resultado$inicio, "\n", sep = "")
cat(
  "Convergencia del optimizador: ",
  mejor_resultado$convergencia,
  "\n",
  sep = ""
)
cat(
  "Log-verosimilitud: ",
  round(mejor_resultado$logLik, 4),
  "\n",
  sep = ""
)

cat("\nCoeficientes constantes:\n")
cat("alpha = ", round(alpha_estimado, 6), "\n", sep = "")
cat("rho   = ", round(rho_estimado, 6), "\n", sep = "")
cat("p-valor de rho = ", format(p_rho, scientific = TRUE), "\n", sep = "")
cat("|rho| < 1: ", abs(rho_estimado) < 1, "\n", sep = "")

cat("\nVarianzas estimadas:\n")
cat("H       = ", format(H_estimado, scientific = TRUE), "\n", sep = "")
cat("Q_beta  = ", format(Q_beta_estimado, scientific = TRUE), "\n", sep = "")
cat(
  "Q_beta/H= ",
  format(Q_beta_estimado / H_estimado, scientific = TRUE),
  "\n",
  sep = ""
)

cat("\nTrayectoria de beta_t:\n")
cat(
  "beta inicial = ",
  round(first(resultados_tvp$beta_t), 5),
  "\n",
  sep = ""
)
cat(
  "beta final   = ",
  round(last(resultados_tvp$beta_t), 5),
  "\n",
  sep = ""
)
cat(
  "cambio total = ",
  round(
    last(resultados_tvp$beta_t) - first(resultados_tvp$beta_t),
    5
  ),
  "\n",
  sep = ""
)

if (exists("tabla_diagnosticos")) {
  cat("\nDiagnósticos:\n")
  print(tabla_diagnosticos)
}

cat("\nResultados guardados en:\n", carpeta_salida, "\n", sep = "")
cat("============================================================\n")

extremos_residuales <- resultados_tvp %>%
  dplyr::filter(is.finite(residuo_estandarizado)) %>%
  dplyr::arrange(desc(abs(residuo_estandarizado))) %>%
  dplyr::select(
    fecha,
    exp_inf_12m,
    exp_inf_24m,
    infl_gt,
    d_exp24,
    d_exp24_lag1,
    shock_inf,
    ajustado_tvp,
    residuo_tvp,
    residuo_estandarizado
  ) %>%
  dplyr::slice_head(n = 10)

print(extremos_residuales)

readr::write_csv(
  extremos_residuales,
  file.path(
    carpeta_salida,
    "14_observaciones_residuales_extremas.csv"
  )
)

resultados_tvp %>%
  dplyr::filter(abs(residuo_estandarizado) > 2) %>%
  dplyr::select(
    fecha,
    d_exp24,
    shock_inf,
    ajustado_tvp,
    residuo_tvp,
    residuo_estandarizado
  )

fecha_extrema <- extremos_residuales$fecha[1]

ventana_extrema <- resultados_tvp %>%
  dplyr::filter(
    fecha >= fecha_extrema %m-% months(3),
    fecha <= fecha_extrema %m+% months(3)
  ) %>%
  dplyr::select(
    fecha,
    exp_inf_12m,
    exp_inf_24m,
    infl_gt,
    d_exp24,
    d_exp24_lag1,
    shock_inf,
    ajustado_tvp,
    residuo_tvp,
    residuo_estandarizado
  )

print(ventana_extrema)
