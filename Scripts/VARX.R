# ============================================================
# TESINA PES (Guatemala)
# ETAPA EXPLORATORIA: VARX de expectativas, inflación y política monetaria
# Autor: Paulo Garrido
# Fecha: 2026-03-11
#
# OBJETIVO:
#   Estimar un VARX pequeño para analizar la dinámica conjunta entre:
#   - expectativas de inflación a 24 meses
#   - inflación observada
#   - tasa de política monetaria
#   - depreciación / tipo de cambio
#
#   incluyendo variables externas como exógenas:
#   - d_effr
#   - d_infl_us
#   - oil_var
#
# INSUMO:
#   - bd_sample_modelo.csv
#
# SALIDAS:
#   - output/varx/seleccion_rezagos.txt
#   - output/varx/resumen_varx.txt
#   - output/varx/irf_exp24_shock_tasa.png
#   - output/varx/irf_inflacion_shock_tasa.png
#   - output/varx/fevd_exp24.png
#   - output/varx/fevd_inflacion.png
# ============================================================

# ------------------------------------------------------------
# 0) Limpieza de entorno
# ------------------------------------------------------------
rm(list = ls())
cat("\014")
if (!is.null(dev.list())) dev.off()

# ------------------------------------------------------------
# 1) Paquetes
# ------------------------------------------------------------
library(readr)
library(dplyr)
library(lubridate)
library(vars)
library(ggplot2)
library(tseries)

select <- dplyr::select

# ------------------------------------------------------------
# 2) Carpetas de trabajo
# ------------------------------------------------------------
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/varx")) dir.create("output/varx")

# ------------------------------------------------------------
# 3) Cargar base
# ------------------------------------------------------------
df_model <- read_csv("bd_sample_modelo.csv", show_col_types = FALSE) |>
    mutate(fecha = ymd(fecha)) |>
    arrange(fecha)

# ------------------------------------------------------------
# 4) Definir muestra de trabajo
#    Ajusta estas fechas según lo que quieras analizar
# ------------------------------------------------------------
usar_submuestra <- TRUE
fecha_inicio <- ymd("2011-01-01")
fecha_fin    <- ymd("2025-12-01")

if (usar_submuestra) {
    df_model <- df_model |>
        filter(fecha >= fecha_inicio & fecha <= fecha_fin)
}

# ------------------------------------------------------------
# 5) Selección de variables
# ------------------------------------------------------------
# Endógenas:
#   exp_inf_12m  -> expectativas 24 meses
#   infl_gt      -> inflación observada
#   tasa_pol     -> tasa de política
#   tcdep_var    -> variación interanual del tipo de cambio
#
# Exógenas:
#   d_effr       -> cambio tasa Fed
#   d_infl_us    -> cambio inflación USA
#   oil_var      -> variación interanual petróleo

df_model <- df_model |>
    arrange(fecha) |>
    mutate(d_tasa_pol_lag3 = lag(d_tasa_pol, 3))

vars_needed <- c(
    "fecha",
    "exp_inf_12m", "infl_gt", "d_tasa_pol", "tcdep_var",
    "d_effr", "d_infl_us", "oil_var"
)


df_var <- df_model |>
    dplyr::select(dplyr::all_of(vars_needed)) |>
    dplyr::filter(dplyr::if_all(-fecha, ~ !is.na(.x)))




cat("Observaciones efectivas para VARX:", nrow(df_var), "\n\n")

# ------------------------------------------------------------
# 6) Preparar matrices
# ------------------------------------------------------------

Y <- df_var |>
    select(exp_inf_12m, infl_gt, d_tasa_pol, tcdep_var)

X_exog <- df_var |>
    select(d_effr, d_infl_us, oil_var)

# Convertir a matrices / ts
Y_mat <- as.matrix(Y)
X_mat <- as.matrix(X_exog)




# ------------------------------------------------------------
# 7) Revisión rápida de estacionariedad (exploratoria)
# ------------------------------------------------------------
# Nota:
# Esto es solo una guía. No se usa automáticamente para transformar.
# Lo hacemos para documentar el comportamiento de las series.

sink("output/varx/tests_adf.txt")
cat("========================================\n")
cat("TESTS ADF EXPLORATORIOS\n")
cat("========================================\n\n")

for (v in colnames(Y_mat)) {
    cat("Variable:", v, "\n")
    print(adf.test(Y_mat[, v]))
    cat("\n-----------------------------\n\n")
}

for (v in colnames(X_mat)) {
    cat("Variable exógena:", v, "\n")
    print(adf.test(X_mat[, v]))
    cat("\n-----------------------------\n\n")
}

sink()

# ------------------------------------------------------------
# 8) Selección de rezagos
# ------------------------------------------------------------
# Se usa hasta 12 rezagos por ser datos mensuales.
# Puedes bajar a 6 si la muestra es corta.

lag_select <- VARselect(
    y = Y_mat,
    lag.max = 12,
    type = "const",
    exogen = X_mat
)

lag_select
sink("output/varx/seleccion_rezagos.txt")
cat("========================================\n")
cat("SELECCIÓN DE REZAGOS VARX\n")
cat("========================================\n\n")
print(lag_select)
sink()

print(lag_select)

# Elegir rezago según AIC, HQ o SC.
# Aquí usaré AIC por defecto; puedes cambiarlo si prefieres parsimonia.
p_opt <- lag_select$selection["AIC(n)"]
p_opt <- as.numeric(p_opt)
p_opt
cat("Rezago seleccionado por AIC:", p_opt, "\n\n")

# ------------------------------------------------------------
# 9) Estimar VARX
# ------------------------------------------------------------
modelo_varx <- VAR(
    y = Y_mat,
    p = p_opt,
    type = "const",
    exogen = X_mat
)

sink("output/varx/resumen_varx.txt")
cat("========================================\n")
cat("RESUMEN VARX\n")
cat("========================================\n\n")
print(summary(modelo_varx))
sink()

print(summary(modelo_varx))

# ------------------------------------------------------------
# 10) Diagnósticos del VARX
# ------------------------------------------------------------

# 10.1 Autocorrelación serial
serial_test <- serial.test(modelo_varx, lags.pt = 12, type = "PT.asymptotic")

# 10.2 Heterocedasticidad
arch_test <- arch.test(modelo_varx, lags.multi = 5)

# 10.3 Normalidad
normality_test <- normality.test(modelo_varx)

sink("output/varx/diagnosticos_varx.txt")
cat("========================================\n")
cat("DIAGNÓSTICOS VARX\n")
cat("========================================\n\n")

cat("---- Serial correlation test ----\n")
print(serial_test)
cat("\n\n")

cat("---- ARCH test ----\n")
print(arch_test)
cat("\n\n")

cat("---- Normality test ----\n")
print(normality_test)
cat("\n")
sink()

# ------------------------------------------------------------
# 11) Estabilidad del sistema
# ------------------------------------------------------------
roots_varx <- roots(modelo_varx, modulus = TRUE)

roots_varx
sink("output/varx/estabilidad_varx.txt")
cat("========================================\n")
cat("RAÍCES DEL VARX\n")
cat("========================================\n\n")
print(roots_varx)
cat("\n\n")
cat("Regla: el sistema es estable si todas las raíces en módulo son < 1.\n")
sink()

print(roots_varx)

# ------------------------------------------------------------
# 12) Impulse Response Functions (IRF)
# ------------------------------------------------------------
# Ordenamiento implícito de Cholesky:
# 1 exp_inf_12m
# 2 infl_gt
# 3 tasa_pol
# 4 tcdep_var
#
# Si luego quieres, lo cambiamos.
#
# Aquí vamos a mirar:
# - respuesta de expectativas a shock de tasa
# - respuesta de inflación a shock de tasa

irf_exp24_tasa <- irf(
    modelo_varx,
    impulse = "d_tasa_pol",
    response = "exp_inf_12m",
    n.ahead = 12,
    boot = TRUE,
    ci = 0.95,
    runs = 500
)

irf_infl_tasa <- irf(
    modelo_varx,
    impulse = "d_tasa_pol",
    response = "infl_gt",
    n.ahead = 12,
    boot = TRUE,
    ci = 0.95,
    runs = 500
)

# ------------------------------------------------------------
# 13) Función auxiliar para graficar IRF
# ------------------------------------------------------------
plot_irf_to_df <- function(irf_obj, impulse_name, response_name) {
    tibble(
        horizonte = 0:(length(irf_obj$irf[[1]]) - 1),
        irf = as.numeric(irf_obj$irf[[1]]),
        lower = as.numeric(irf_obj$Lower[[1]]),
        upper = as.numeric(irf_obj$Upper[[1]])
    ) |>
        mutate(
            impulse = impulse_name,
            response = response_name
        )
}

df_irf_exp24_tasa <- plot_irf_to_df(irf_exp24_tasa, "d_tasa_pol", "exp_inf_12m")
df_irf_infl_tasa  <- plot_irf_to_df(irf_infl_tasa, "d_tasa_pol", "infl_gt")

# ------------------------------------------------------------
# 14) Gráficos IRF
# ------------------------------------------------------------
p_irf_exp24_tasa <- ggplot(df_irf_exp24_tasa, aes(x = horizonte, y = irf)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "IRF: respuesta de expectativas a shock de tasa de política",
        subtitle = "Variable respuesta: exp_inf_12m",
        x = "Horizonte (meses)",
        y = "Respuesta"
    ) +
    theme_minimal()

ggsave(
    "output/varx/irf_exp24_shock_tasa.png",
    plot = p_irf_exp24_tasa,
    width = 10,
    height = 5,
    dpi = 300
)

p_irf_infl_tasa <- ggplot(df_irf_infl_tasa, aes(x = horizonte, y = irf)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "IRF: respuesta de inflación a shock de tasa de política",
        subtitle = "Variable respuesta: infl_gt",
        x = "Horizonte (meses)",
        y = "Respuesta"
    ) +
    theme_minimal()

ggsave(
    "output/varx/irf_inflacion_shock_tasa.png",
    plot = p_irf_infl_tasa,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_irf_exp24_tasa)
print(p_irf_infl_tasa)

# ------------------------------------------------------------
# 15) Forecast Error Variance Decomposition (FEVD)
# ------------------------------------------------------------
fevd_varx <- fevd(modelo_varx, n.ahead = 12)

# Extraer FEVD de expectativas
df_fevd_exp24 <- as.data.frame(fevd_varx$exp_inf_12m)
df_fevd_exp24$horizonte <- 1:nrow(df_fevd_exp24)

df_fevd_exp24_long <- reshape(
    df_fevd_exp24,
    varying = colnames(df_fevd_exp24)[1:(ncol(df_fevd_exp24)-1)],
    v.names = "proporcion",
    timevar = "shock",
    times = colnames(df_fevd_exp24)[1:(ncol(df_fevd_exp24)-1)],
    direction = "long"
)

p_fevd_exp24 <- ggplot(df_fevd_exp24_long,
                       aes(x = horizonte, y = proporcion, fill = shock)) +
    geom_area(position = "stack") +
    labs(
        title = "FEVD de expectativas de inflación a 24 meses",
        x = "Horizonte (meses)",
        y = "Proporción de varianza"
    ) +
    theme_minimal()

ggsave(
    "output/varx/fevd_exp24.png",
    plot = p_fevd_exp24,
    width = 10,
    height = 5,
    dpi = 300
)

# Extraer FEVD de inflación
df_fevd_infl <- as.data.frame(fevd_varx$infl_gt)
df_fevd_infl$horizonte <- 1:nrow(df_fevd_infl)

df_fevd_infl_long <- reshape(
    df_fevd_infl,
    varying = colnames(df_fevd_infl)[1:(ncol(df_fevd_infl)-1)],
    v.names = "proporcion",
    timevar = "shock",
    times = colnames(df_fevd_infl)[1:(ncol(df_fevd_infl)-1)],
    direction = "long"
)

p_fevd_infl <- ggplot(df_fevd_infl_long,
                      aes(x = horizonte, y = proporcion, fill = shock)) +
    geom_area(position = "stack") +
    labs(
        title = "FEVD de inflación observada",
        x = "Horizonte (meses)",
        y = "Proporción de varianza"
    ) +
    theme_minimal()

ggsave(
    "output/varx/fevd_inflacion.png",
    plot = p_fevd_infl,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_fevd_exp24)
print(p_fevd_infl)

# ------------------------------------------------------------
# 16) Guardar insumo final usado
# ------------------------------------------------------------
write_csv(df_var, "output/varx/base_varx_usada.csv")

cat("\nScript VARX completado con éxito.\n")

