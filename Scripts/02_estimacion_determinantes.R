# ============================================================
# TESINA PES (Guatemala)
# ETAPA 1: Estimación del modelo baseline (referencia FMI)
# Autor: Paulo Garrido
# Fecha: 2026-03-09
#
# OBJETIVO:
#   Estimar un modelo lineal baseline para identificar los
#   determinantes promedio de las expectativas de inflación
#   a 24 meses en Guatemala.
#
# INSUMO:
#   - bd_sample_modelo.csv
#
# SALIDAS:
#   - output/resumen_modelo_baseline.txt
#   - output/coeficientes_baseline.csv
#   - output/df_model_used.csv
#   - output/grafico_residuos_tiempo.png
#   - output/grafico_actual_vs_ajustado.png
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
library(ggplot2)
library(broom)
library(car)
library(lmtest)
library(sandwich)

# ------------------------------------------------------------
# 2) Carpetas de trabajo
# ------------------------------------------------------------
if (!dir.exists("output")) dir.create("output")

# ------------------------------------------------------------
# 3) Cargar base de datos
# ------------------------------------------------------------
df_model <- read_csv("bd_sample_modelo.csv", show_col_types = FALSE)

# ------------------------------------------------------------
# 3.1) Definir período de estimación
# ------------------------------------------------------------
fecha_inicio <- ymd("2015-01-01")
fecha_fin    <- ymd("2025-11-01")

df_model <- df_model |>
    mutate(fecha = ymd(fecha)) |>
    filter(fecha >= fecha_inicio & fecha <= fecha_fin) |>
    arrange(fecha)


# Guardar la muestra efectivamente usada
write_csv(df_model, "output/df_model_used.csv")

# ------------------------------------------------------------
# 5) Definición del modelo baseline
# ------------------------------------------------------------
# Variable dependiente principal:
#   exp_inf_24m = expectativas de inflación a 24 meses
#
# Especificación baseline inspirada en el working paper del FMI:
#   exp_inf_24m_t = f(
#      infl_gt_lag1,
#      d_tasa_pol,
#      exp24_lag1,
#      imae_lag2,
#      tcdep_var,
#      deficit_fiscal,
#      d_effr,
#      d_infl_us,
#      oil_var
#   )

modelo_baseline <- exp_inf_24m ~
    infl_gt_lag1 +
    d_tasa_pol +
    exp24_lag1 +
    imae_lag2 +
    tcdep_var +
    deficit_fiscal +
    d_effr +
    d_infl_us +
    oil_var

modelo_baseline <- lm(modelo_baseline, data = df_model)
summary(modelo_baseline)

# ------------------------------------------------------------
# 6) Estimación OLS del baseline
# ------------------------------------------------------------
modelo_baseline <- lm(modelo_baseline, data = df_model)
summary(modelo_baseline)


# ------------------------------------------------------------
# 7) Errores estándar robustos
# ------------------------------------------------------------
# Para una serie macro mensual, conviene al menos reportar
# errores robustos a heterocedasticidad.

vcov_robust <- vcovHC(modelo_baseline, type = "HC1")
resultados_robustos <- coeftest(modelo_baseline, vcov = vcov_robust)
resultados_robustos
cat("========================================\n")
cat("RESULTADOS CON ERRORES ESTÁNDAR ROBUSTOS (HC1)\n")
cat("========================================\n\n")
print(resultados_robustos)
cat("\n")

library(sandwich)
library(lmtest)

vcov_nw <- NeweyWest(modelo_baseline, lag = 4, prewhite = FALSE, adjust = TRUE)
resultados_nw <- coeftest(modelo_baseline, vcov = vcov_nw)

print(resultados_nw)

# ------------------------------------------------------------
# 8) Guardar tablas de resultados
# ------------------------------------------------------------
# Coeficientes OLS estándar
tabla_ols <- tidy(modelo_baseline)

# Coeficientes con errores robustos
tabla_robusta <- tidy(coeftest(modelo_baseline, vcov = vcov_robust)) |>
    rename(
        estimate = estimate,
        std.error = std.error,
        statistic = statistic,
        p.value = p.value
    )

write_csv(tabla_ols, "output/coeficientes_baseline_ols.csv")
write_csv(tabla_robusta, "output/coeficientes_baseline_robustos.csv")

# ------------------------------------------------------------
# 9) Diagnósticos básicos del modelo
# ------------------------------------------------------------

# 9.1 Multicolinealidad: VIF
vif_baseline <- vif(modelo_baseline)
vif_baseline
#•	VIF < 5 → sin mayor preocupación;
#•	VIF entre 5 y 10 → atención;
#•	VIF > 10 → problema serio.

# 9.2 Heterocedasticidad: Breusch-Pagan
bp_test <- bptest(modelo_baseline)
bp_test
cat("========================================\n")
cat("TEST DE BREUSCH-PAGAN\n")
cat("========================================\n\n")
print(bp_test)
cat("\n")

# 9.3 Autocorrelación: Durbin-Watson
dw_test <- dwtest(modelo_baseline)

cat("========================================\n")
cat("TEST DE DURBIN-WATSON\n")
cat("========================================\n\n")
print(dw_test)
cat("\n")

# 9.4 Normalidad de residuos: Jarque-Bera
jb_test <- jarque.bera.test(residuals(modelo_baseline))

cat("========================================\n")
cat("TEST DE JARQUE-BERA\n")
cat("========================================\n\n")
print(jb_test)
cat("\n")

# Guardar diagnóstico en txt adicional
sink("output/tests_diagnostico_baseline.txt")
cat("========================================\n")
cat("DIAGNÓSTICOS DEL MODELO BASELINE\n")
cat("========================================\n\n")

cat("Breusch-Pagan:\n")
print(bp_test)
cat("\n\n")

cat("Durbin-Watson:\n")
print(dw_test)
cat("\n\n")

cat("Jarque-Bera:\n")
print(jb_test)
cat("\n\n")

cat("VIF:\n")
print(vif_baseline)
cat("\n")
sink()

# ------------------------------------------------------------
# 10) Valores ajustados y residuos
# ------------------------------------------------------------
df_model <- df_model |>
    mutate(
        y_hat = fitted(modelo_baseline),
        resid = residuals(modelo_baseline)
    )

# ------------------------------------------------------------
# 11) Gráficos de diagnóstico
# ------------------------------------------------------------

# 11.1 Residuos en el tiempo
p1 <- ggplot(df_model, aes(x = fecha, y = resid)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "Residuos del modelo baseline en el tiempo",
        x = "Fecha",
        y = "Residuo"
    ) +
    theme_minimal()

ggsave(
    filename = "output/grafico_residuos_tiempo.png",
    plot = p1,
    width = 10,
    height = 5,
    dpi = 300
)

# 11.2 Observado vs ajustado
p2 <- ggplot(df_model, aes(x = fecha)) +
    geom_line(aes(y = exp_inf_24m, linetype = "Observado")) +
    geom_line(aes(y = y_hat, linetype = "Ajustado")) +
    labs(
        title = "Expectativas de inflación a 24 meses: observado vs ajustado",
        x = "Fecha",
        y = "Expectativa de inflación (%)",
        linetype = ""
    ) +
    theme_minimal()

ggsave(
    filename = "output/grafico_actual_vs_ajustado.png",
    plot = p2,
    width = 10,
    height = 5,
    dpi = 300
)

# ------------------------------------------------------------
# 12) Guardar base con fitted values y residuos
# ------------------------------------------------------------
write_csv(df_model, "output/df_model_with_fitted_resid.csv")

# ------------------------------------------------------------
# 13) Interpretación básica automática
# ------------------------------------------------------------
cat("========================================\n")
cat("INTERPRETACIÓN PRELIMINAR\n")
cat("========================================\n\n")

cat("1. Este modelo estima los determinantes promedio de las expectativas\n")
cat("   de inflación a 24 meses en Guatemala.\n\n")

cat("2. El coeficiente de exp24_lag1 mide la persistencia/inercia de las\n")
cat("   expectativas de mediano plazo.\n\n")

cat("3. El coeficiente de infl_gt_lag1 permite evaluar qué tanto las\n")
cat("   expectativas siguen a la inflación observada pasada, lo cual es\n")
cat("   central para discutir el anclaje.\n\n")

cat("4. El coeficiente de d_tasa_pol muestra si cambios en la tasa de\n")
cat("   política monetaria se asocian con movimientos en las expectativas.\n\n")

cat("5. Las variables externas (d_effr, d_infl_us, oil_var) capturan el\n")
cat("   entorno internacional y su posible transmisión a expectativas.\n\n")

cat("6. Esta etapa NO prueba aún asimetrías; solo establece el baseline\n")
cat("   lineal de referencia.\n\n")

cat("Script completado con éxito.\n")


meta_central <- 4
meta_inf     <- 3
meta_sup     <- 5

meta_central

df_model <- df_model |>
    mutate(
        gap_meta     = infl_gt - meta_central,
        gap_meta_pos = pmax(gap_meta, 0),
        gap_meta_neg = pmax(-gap_meta, 0),
        dentro_banda = if_else(infl_gt >= meta_inf & infl_gt <= meta_sup, 1, 0),
        sobre_banda  = if_else(infl_gt > meta_sup, 1, 0),
        bajo_banda   = if_else(infl_gt < meta_inf, 1, 0)
    )

modelo_banda <- lm(
    exp_inf_24m ~ exp24_lag1 + sobre_banda + bajo_banda +
        lag(d_tasa_pol, 3) + imae_lag2 + tcdep_var +
        deficit_fiscal + d_effr + d_infl_us + oil_var,
    data = df_model
)

summary(modelo_banda)

linearHypothesis(modelo_banda, "sobre_banda + bajo_banda = 0")

modelo_simetria <- lm(
    exp_inf_24m ~ exp24_lag1 + gap_meta_pos + gap_meta_neg +
        lag(d_tasa_pol, 3) + imae_lag2 + tcdep_var +
        deficit_fiscal + d_effr + d_infl_us + oil_var,
    data = df_model
)

summary(modelo_simetria)

library(car)

linearHypothesis(modelo_simetria, "gap_meta_pos + gap_meta_neg = 0")


#magnitud
meta_inf <- 3
meta_sup <- 5

df_model <- df_model |>
    mutate(
        gap_sobre_banda = pmax(infl_gt - meta_sup, 0),
        gap_bajo_banda  = pmax(meta_inf - infl_gt, 0)
    )

modelo_banda_gap <- lm(
    exp_inf_24m ~ exp24_lag1 + gap_sobre_banda + gap_bajo_banda +
        lag(d_tasa_pol, 3) + imae_lag2 + tcdep_var +
        deficit_fiscal + d_effr + d_infl_us + oil_var,
    data = df_model
)

summary(modelo_banda_gap)

library(car)

linearHypothesis(modelo_banda_gap, "gap_sobre_banda + gap_bajo_banda = 0")


df_model <- df_model |>
    mutate(
        infl_alta = if_else(infl_gt > 7, 1, 0),
        infl_lag1_x_alta = infl_gt_lag1 * infl_alta
    )

modelo_desanclaje <- lm(
    exp_inf_24m ~ exp24_lag1 + infl_gt_lag1 + infl_alta + infl_lag1_x_alta +
        lag(d_tasa_pol, 3) + imae_lag2 + tcdep_var +
        deficit_fiscal + d_effr + d_infl_us + oil_var,
    data = df_model
)

summary(modelo_desanclaje)
