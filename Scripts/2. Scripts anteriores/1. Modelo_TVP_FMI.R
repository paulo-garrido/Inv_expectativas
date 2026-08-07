# ============================================================
# TESINA PES - MODELO BASE FMI / WORLD BANK
# Expectativas de inflación a 24 meses - Guatemala
# ============================================================

# ------------------------------------------------------------
# 0) Limpieza de entorno
# ------------------------------------------------------------
rm(list = ls())
gc()
cat("\014")

# ------------------------------------------------------------
# 1) Paquetes
# ------------------------------------------------------------
paquetes <- c(
    "readr",
    "dplyr",
    "lubridate",
    "ggplot2",
    "broom",
    "lmtest",
    "sandwich",
    "car"
)

instalados <- rownames(installed.packages())

for (p in paquetes) {
    if (!(p %in% instalados)) install.packages(p)
}

invisible(lapply(paquetes, library, character.only = TRUE))

# ------------------------------------------------------------
# 2) Cargar base de datos
#    Cambia la ruta según donde guardes el archivo
# ------------------------------------------------------------
ruta_archivo <- "bd_sample_modelo.csv"

df <- read_csv(ruta_archivo, show_col_types = FALSE)

# ------------------------------------------------------------
# 3) Inspección inicial
# ------------------------------------------------------------
glimpse(df)
summary(df)

# ------------------------------------------------------------
# 4) Preparación de datos
# ------------------------------------------------------------
df <- df %>%
    mutate(
        fecha = as.Date(fecha)
    ) %>%
    arrange(fecha)

# Crear rezago de 3 meses para el cambio en tasa de política
# Se usa dplyr::lag explícitamente para evitar conflictos
df <- df %>%
    mutate(
        d_tasa_pol_lag3 = dplyr::lag(d_tasa_pol, 3)
    )

    # ------------------------------------------------------------
# 5) Selección de variables del modelo
#    Modelo base tipo FMI / World Bank:
#    exp_inf_24m_t = f(exp24_lag1, infl_gt_lag1, política monetaria,
#                      actividad, tipo de cambio, fiscal, factores externos)
# ------------------------------------------------------------
vars_modelo <- c(
    "fecha",
    "exp_inf_24m",
    "exp24_lag1",
    "infl_gt_lag1",
    "d_tasa_pol_lag3",
    "imae_lag2",
    "tcdep_var",
    "deficit_fiscal",
    "d_effr",
    "d_infl_us",
    "oil_var"
)

df_modelo <- df %>%
    dplyr::select(all_of(vars_modelo)) %>%
    filter(complete.cases(.))

# Verificar tamaño de muestra final
cat("Observaciones utilizadas en el modelo:", nrow(df_modelo), "\n")

# ------------------------------------------------------------
# 6) Estimación del modelo base
# ------------------------------------------------------------
modelo_fmi <- lm(
    exp_inf_24m ~ exp24_lag1 +
        infl_gt_lag1 +
        d_tasa_pol_lag3 +
        imae_lag2 +
        tcdep_var +
        deficit_fiscal +
        d_effr +
        d_infl_us +
        oil_var,
    data = df_modelo
)

# ------------------------------------------------------------
# 7) Resultados básicos
# ------------------------------------------------------------
cat("\n================ RESUMEN DEL MODELO ================\n")
summary(modelo_fmi)

cat("\n================ COEFICIENTES ORDENADOS ================\n")
broom::tidy(modelo_fmi)

# ------------------------------------------------------------
# 8) Errores estándar robustos
# ------------------------------------------------------------
cat("\n================ ERRORES ROBUSTOS HC1 ================\n")
coeftest(modelo_fmi, vcov. = vcovHC(modelo_fmi, type = "HC1"))

# ------------------------------------------------------------
# 9) Diagnósticos econométricos
# ------------------------------------------------------------

# 9.1 Heterocedasticidad: Breusch-Pagan
cat("\n================ BREUSCH-PAGAN ================\n")
bptest(modelo_fmi)

# 9.2 Autocorrelación: Breusch-Godfrey
cat("\n================ BREUSCH-GODFREY ================\n")
bgtest(modelo_fmi, order = 4)

# 9.3 Multicolinealidad: VIF
cat("\n================ VIF ================\n")
car::vif(modelo_fmi)

# ------------------------------------------------------------
# 10) Residuos y valores ajustados
# ------------------------------------------------------------
df_modelo <- df_modelo %>%
    mutate(
        fitted   = fitted(modelo_fmi),
        residuo  = resid(modelo_fmi)
    )

# ------------------------------------------------------------
# 11) Gráficas de ajuste
# ------------------------------------------------------------

# Serie observada vs ajustada
graf_ajuste <- ggplot(df_modelo, aes(x = fecha)) +
    geom_line(aes(y = exp_inf_24m, color = "Observada"), linewidth = 1) +
    geom_line(aes(y = fitted, color = "Ajustada"), linewidth = 1, linetype = "dashed") +
    labs(
        title = "Expectativa de inflación a 24 meses: observada vs ajustada",
        x = NULL,
        y = "%",
        color = NULL
    ) +
    theme_minimal()

print(graf_ajuste)

# Residuos en el tiempo
graf_residuos <- ggplot(df_modelo, aes(x = fecha, y = residuo)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "Residuos del modelo",
        x = NULL,
        y = "Residuo"
    ) +
    theme_minimal()

print(graf_residuos)

# Histograma de residuos
graf_hist <- ggplot(df_modelo, aes(x = residuo)) +
    geom_histogram(bins = 20) +
    labs(
        title = "Histograma de residuos",
        x = "Residuo",
        y = "Frecuencia"
    ) +
    theme_minimal()

print(graf_hist)

# ------------------------------------------------------------
# 12) Tabla final de resultados
# ------------------------------------------------------------
tabla_resultados <- broom::tidy(modelo_fmi) %>%
    mutate(
        signif = case_when(
            p.value < 0.01 ~ "***",
            p.value < 0.05 ~ "**",
            p.value < 0.10 ~ "*",
            TRUE ~ ""
        )
    )

print(tabla_resultados)

# ------------------------------------------------------------
# 13) Guardar resultados opcionalmente
# ------------------------------------------------------------
write.csv(tabla_resultados, "resultados_modelo_fmi.csv", row.names = FALSE)
write.csv(df_modelo, "base_con_ajuste_modelo_fmi.csv", row.names = FALSE)

cat("\nScript ejecutado correctamente.\n")