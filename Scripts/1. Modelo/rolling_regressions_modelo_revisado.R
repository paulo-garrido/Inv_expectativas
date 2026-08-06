# ============================================================
# ROLLING REGRESSIONS: MODELO REVISADO
# Anclaje de las expectativas de inflación en Guatemala
# ============================================================
#
# Especificación principal:
#
#   (e24_t - meta_t) = alpha_w
#                    + rho_w  (e24_{t-1} - meta_{t-1})
#                    + beta_w (pi_{t-1}  - meta_{t-1})
#                    + gamma_w q_{t-1}
#                    + u_t
#
# El coeficiente de interés es beta_w. Una reducción de su
# magnitud se interpreta como evidencia compatible con un
# fortalecimiento del anclaje de las expectativas.
#
# Ventana principal: 36 meses.
# Robustez: 30 meses.
# Inferencia: errores estándar Newey-West HAC.
# ============================================================

# ------------------------------------------------------------
# 0) Paquetes y limpieza del entorno
# ------------------------------------------------------------
paquetes <- c(
  "dplyr", "tidyr", "readr", "lubridate", "ggplot2",
  "purrr", "sandwich", "lmtest"
)

faltantes <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltantes) > 0) {
  stop(
    "Faltan los siguientes paquetes: ", paste(faltantes, collapse = ", "),
    ". Instálelos con install.packages(c(",
    paste(sprintf('"%s"', faltantes), collapse = ", "), "))."
  )
}

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ggplot2)
library(purrr)
library(sandwich)
library(lmtest)

rm(list = setdiff(ls(), c("paquetes", "faltantes")))
cat("\014")
if (!is.null(dev.list())) dev.off()

# ------------------------------------------------------------
# 1) Configuración
# ------------------------------------------------------------
archivo_datos <- "bd_sample_modelo.csv"
carpeta_salida <- "output"

# Nombres esperados de las variables principales.
col_fecha <- "fecha"
col_exp24 <- "exp_inf_24m"
col_inflacion <- "infl_gt"

# El script buscará automáticamente la primera variable de tipo
# de cambio disponible entre estas opciones. Agregue aquí el
# nombre de su columna si es distinto.
candidatas_tc <- c(
  "tcdep_var", "depreciacion_tc", "tc_var", "q_t",
  "tipo_cambio_var", "dep_tc"
)

# Ventanas móviles
ventana_principal <- 36L
ventana_robustez <- 30L

# Rezago para la matriz HAC Newey-West.
# Tres rezagos constituyen una elección parsimoniosa para datos
# mensuales y ventanas relativamente cortas.
hac_lag <- 3L

# ------------------------------------------------------------
# 2) Funciones auxiliares
# ------------------------------------------------------------
primera_columna_disponible <- function(datos, candidatas) {
  disponibles <- candidatas[candidatas %in% names(datos)]
  if (length(disponibles) == 0) return(NA_character_)
  disponibles[[1]]
}

verificar_columnas <- function(datos, columnas) {
  faltan <- setdiff(columnas, names(datos))
  if (length(faltan) > 0) {
    stop(
      "No se encontraron las siguientes columnas en la base: ",
      paste(faltan, collapse = ", "),
      ". Revise la sección 1) Configuración."
    )
  }
}

# ------------------------------------------------------------
# 3) Cargar y ordenar la base
# ------------------------------------------------------------
if (!file.exists(archivo_datos)) {
  stop("No se encontró el archivo: ", archivo_datos)
}

df_model <- read_csv(archivo_datos, show_col_types = FALSE)

verificar_columnas(df_model, c(col_fecha, col_exp24, col_inflacion))

col_tc <- primera_columna_disponible(df_model, candidatas_tc)
if (is.na(col_tc)) {
  stop(
    "No se encontró una variable de depreciación cambiaria. ",
    "Agregue su nombre al vector 'candidatas_tc'. Columnas disponibles: ",
    paste(names(df_model), collapse = ", ")
  )
}

message("Variable de tipo de cambio seleccionada: ", col_tc)

df_model <- df_model |>
  mutate(fecha = ymd(.data[[col_fecha]])) |>
  arrange(fecha)

if (anyNA(df_model$fecha)) {
  stop("Existen fechas que no pudieron convertirse al formato Date.")
}

# ------------------------------------------------------------
# 4) Construcción de la meta y de las brechas
# ------------------------------------------------------------
# Si la base ya contiene 'meta_inflacion', se utiliza esa serie.
# En caso contrario, se construye según la trayectoria vigente
# durante la muestra:
#   2010-2011: 5.0%
#   2012:      4.5%
#   2013-...:  4.0%

if (!"meta_inflacion" %in% names(df_model)) {
  df_model <- df_model |>
    mutate(
      meta_inflacion = case_when(
        fecha < ymd("2012-01-01") ~ 5.0,
        fecha < ymd("2013-01-01") ~ 4.5,
        TRUE ~ 4.0
      )
    )
} else {
  message("Se utilizará la columna existente 'meta_inflacion'.")
}

df_est <- df_model |>
  transmute(
    fecha = fecha,
    exp_inf_24m = .data[[col_exp24]],
    infl_gt = .data[[col_inflacion]],
    meta_inflacion = meta_inflacion,
    tcdep = .data[[col_tc]],

    # Variables expresadas como desviaciones respecto de la meta
    exp_gap_24 = exp_inf_24m - meta_inflacion,
    infl_gap = infl_gt - meta_inflacion,

    # Rezagos incluidos en la especificación
    exp_gap_24_lag1 = lag(exp_gap_24, 1),
    infl_gap_lag1 = lag(infl_gap, 1),
    tcdep_lag1 = lag(tcdep, 1)
  ) |>
  drop_na(
    fecha, exp_gap_24, exp_gap_24_lag1,
    infl_gap_lag1, tcdep_lag1
  )

if (nrow(df_est) < ventana_principal) {
  stop(
    "La muestra efectiva tiene ", nrow(df_est),
    " observaciones, menos que la ventana principal de ",
    ventana_principal, " meses."
  )
}

# Verificación informativa de continuidad mensual
saltos_meses <- diff(year(df_est$fecha) * 12 + month(df_est$fecha))
if (any(saltos_meses != 1)) {
  warning(
    "La muestra efectiva presenta saltos entre meses. ",
    "Revise observaciones faltantes antes de interpretar las ventanas móviles."
  )
}

# ------------------------------------------------------------
# 5) Especificación econométrica
# ------------------------------------------------------------
formula_principal <- exp_gap_24 ~
  exp_gap_24_lag1 + infl_gap_lag1 + tcdep_lag1

coeficiente_interes <- "infl_gap_lag1"

# ------------------------------------------------------------
# 6) Función de estimación rolling con errores HAC
# ------------------------------------------------------------
estimar_ventana <- function(data_window, window_size, hac_lag = 3L) {
  modelo <- lm(formula_principal, data = data_window)

  vcov_hac <- sandwich::NeweyWest(
    modelo,
    lag = hac_lag,
    prewhite = FALSE,
    adjust = TRUE
  )

  tabla_hac <- lmtest::coeftest(modelo, vcov. = vcov_hac)

  if (!coeficiente_interes %in% rownames(tabla_hac)) {
    stop("No fue posible estimar el coeficiente ", coeficiente_interes, ".")
  }

  beta <- tabla_hac[coeficiente_interes, "Estimate"]
  se_hac <- tabla_hac[coeficiente_interes, "Std. Error"]
  estadistico <- tabla_hac[coeficiente_interes, "t value"]
  p_valor <- tabla_hac[coeficiente_interes, "Pr(>|t|)"]

  tibble(
    fecha_inicio = min(data_window$fecha),
    fecha_fin = max(data_window$fecha),
    ventana = window_size,
    n_obs = nobs(modelo),
    estimate = unname(beta),
    std.error = unname(se_hac),
    statistic = unname(estadistico),
    p.value = unname(p_valor),
    ic_inf = unname(beta - qnorm(0.975) * se_hac),
    ic_sup = unname(beta + qnorm(0.975) * se_hac),
    r_squared = summary(modelo)$r.squared,
    adj_r_squared = summary(modelo)$adj.r.squared
  )
}

correr_rolling <- function(datos, window_size, hac_lag = 3L) {
  if (nrow(datos) < window_size) {
    stop("La ventana solicitada excede la muestra disponible.")
  }

  purrr::map_dfr(seq.int(window_size, nrow(datos)), function(i) {
    data_window <- datos[(i - window_size + 1):i, , drop = FALSE]
    estimar_ventana(data_window, window_size, hac_lag)
  })
}

# ------------------------------------------------------------
# 7) Estimaciones: ventana principal y robustez
# ------------------------------------------------------------
rolling_36m <- correr_rolling(
  datos = df_est,
  window_size = ventana_principal,
  hac_lag = hac_lag
)

rolling_30m <- correr_rolling(
  datos = df_est,
  window_size = ventana_robustez,
  hac_lag = hac_lag
)

rolling_comparacion <- bind_rows(
  rolling_36m |> mutate(especificacion = "Ventana principal: 36 meses"),
  rolling_30m |> mutate(especificacion = "Robustez: 30 meses")
)

# ------------------------------------------------------------
# 8) Guardar resultados
# ------------------------------------------------------------
if (!dir.exists(carpeta_salida)) dir.create(carpeta_salida, recursive = TRUE)

write_csv(
  rolling_36m,
  file.path(carpeta_salida, "rolling_modelo_principal_36m.csv")
)

write_csv(
  rolling_30m,
  file.path(carpeta_salida, "rolling_robustez_30m.csv")
)

write_csv(
  rolling_comparacion,
  file.path(carpeta_salida, "rolling_comparacion_36m_30m.csv")
)

write_csv(
  df_est,
  file.path(carpeta_salida, "datos_utilizados_rolling.csv")
)

# ------------------------------------------------------------
# 9) Gráfico principal: ventana de 36 meses
# ------------------------------------------------------------
p_roll_36 <- ggplot(rolling_36m, aes(x = fecha_fin, y = estimate)) +
  geom_ribbon(aes(ymin = ic_inf, ymax = ic_sup), alpha = 0.20) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(
    xintercept = as.numeric(ymd(c("2012-01-01", "2013-01-01"))),
    linetype = "dotted",
    linewidth = 0.5
  ) +
  labs(
    title = "Sensibilidad de las expectativas a la brecha de inflación",
    subtitle = paste0(
      "Regresiones móviles de ", ventana_principal,
      " meses; bandas de confianza HAC al 95%"
    ),
    x = "Fecha final de la ventana",
    y = expression(beta[w])
  ) +
  theme_minimal()

ggsave(
  file.path(carpeta_salida, "rolling_modelo_principal_36m.png"),
  plot = p_roll_36,
  width = 10,
  height = 5.5,
  dpi = 300
)

# ------------------------------------------------------------
# 10) Gráfico de robustez: comparación 36 vs. 30 meses
# ------------------------------------------------------------
p_comparacion <- ggplot(
  rolling_comparacion,
  aes(x = fecha_fin, y = estimate, linetype = especificacion)
) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Robustez a la longitud de la ventana móvil",
    subtitle = "Coeficiente de la brecha de inflación rezagada",
    x = "Fecha final de la ventana",
    y = expression(beta[w]),
    linetype = NULL
  ) +
  theme_minimal()

ggsave(
  file.path(carpeta_salida, "rolling_comparacion_36m_30m.png"),
  plot = p_comparacion,
  width = 10,
  height = 5.5,
  dpi = 300
)

print(p_roll_36)
print(p_comparacion)

# ------------------------------------------------------------
# 11) Resumen en consola
# ------------------------------------------------------------
cat("\n============================================================\n")
cat("MODELO ROLLING FINALIZADO\n")
cat("============================================================\n")
cat("Variable dependiente: brecha de expectativa a 24 meses\n")
cat("Coeficiente de interés: brecha de inflación rezagada\n")
cat("Control dinámico: brecha de expectativa rezagada\n")
cat("Control externo: depreciación cambiaria rezagada\n")
cat("Ventana principal:", ventana_principal, "meses\n")
cat("Ventana de robustez:", ventana_robustez, "meses\n")
cat("Errores estándar: Newey-West HAC, lag =", hac_lag, "\n")
cat("Resultados guardados en:", carpeta_salida, "\n")
cat("============================================================\n")


