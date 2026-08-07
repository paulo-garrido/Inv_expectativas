# ============================================================
# TESINA PES (Guatemala)
# FASE 1: Medición del anclaje de expectativas de inflación
# Base real: bd.xlsx | hoja: BD
# Método: Rolling regressions sobre exp_24
# ============================================================

# ------------------------------------------------------------
# 0) Limpieza de entorno
# ------------------------------------------------------------
rm(list = ls())
gc()
cat("\014")
if (!is.null(dev.list())) dev.off()
# ------------------------------------------------------------
# 1) Carga de paquetes
# ------------------------------------------------------------
packages <- c(
    "readxl",
    "dplyr",
    "tidyr",
    "lubridate",
    "ggplot2",
    "purrr",
    "broom",
    "zoo",
    "stringr",
    "scales",
    "writexl"
    )

lag <- dplyr::lag

installed <- rownames(installed.packages())
for (p in packages) {
    if (!(p %in% installed)) install.packages(p)
}
invisible(lapply(packages, library, character.only = TRUE))

# ------------------------------------------------------------
# 2) Parámetros de trabajo
# ------------------------------------------------------------
ruta_archivo <- "bd.xlsx"   

hoja_datos  <- "BD"
window_size <- 60   # podemos probar también 48

# Crear carpetas de salida
#dir.create("output", showWarnings = FALSE)
#dir.create("output/graficos", showWarnings = FALSE)

# ------------------------------------------------------------
# 3) Carga de base
# ------------------------------------------------------------
df_raw <- read_excel(ruta_archivo, sheet = hoja_datos)


# ------------------------------------------------------------
# 4) Renombrar columnas
# ------------------------------------------------------------
df <- df_raw %>%
    rename(
        fecha            = fecha,
        exp_12           = exp_12,
        exp_24           = exp_24,
        infl             =  inf,
        tasa_politica    = Tasa_pol,
        imae             = imae,
        deficit_pct_pib  = deficit_pct_pib,
        tc               = tc,
        ff_rates         = ff_rates,
        inf_usa          = inf_usa,
        wti_usd_bbl      = wti_usd_bbl,
        def_fisc         = def_fisc
    )

# ------------------------------------------------------------
# 5) Formato de fecha, orden y filtro
# ------------------------------------------------------------
df <- df %>%
    mutate(
        fecha = as.Date(fecha)
    ) %>%
    arrange(fecha) %>%
    filter(
        fecha >= ymd("2011-01-01"),
        fecha <= ymd("2025-02-28")
    )

# ------------------------------------------------------------
# 6) Conversión numérica
# ------------------------------------------------------------
df <- df %>%
    mutate(
        across(
            .cols = c(
                exp_12, exp_24, infl, tasa_politica, imae,
                deficit_pct_pib, tc, ff_rates, inf_usa, wti_usd_bbl, def_fisc
            ),
            .fns = ~ as.numeric(.)
        )
    )

# ------------------------------------------------------------
# 7) Meta de inflación
# ------------------------------------------------------------
# Aquí la dejaremos fija en 4.0 como aproximación inicial.
# Luego, si quisieramos podriamos hacerla variable por período.
df <- df %>%
    mutate(meta_inf = 4)

# ------------------------------------------------------------
# 8) Transformación de variables
# ------------------------------------------------------------

df <- df %>%
    mutate(
        # Desalineamiento respecto de la meta
        gap_exp24     = exp_24 - meta_inf,
        abs_gap_exp24 = abs(exp_24 - meta_inf),
        
        # Tipo de cambio
        l_tc         = log(tc),
        d_l_tc       = l_tc - lag(l_tc, 1),
        tc_yoy       = 100 * (l_tc - lag(l_tc, 12)),
        
        # Petróleo
        l_wti        = log(wti_usd_bbl),
        d_l_wti      = l_wti - lag(l_wti, 1),
        wti_yoy      = 100 * (l_wti - lag(l_wti, 12)),
        
        # Actividad
        l_imae      = log(imae),
        imae_yoy     = 100 * (log(imae) - lag(log(imae), 12)),
        
        # Cambios / rezagos adicionales
        d_exp12      = exp_12 - lag(exp_12, 1),
        d_exp24      = exp_24 - lag(exp_24, 1),
        d_inflacion  = infl - lag(infl, 1),
        d_ff         = ff_rates - lag(ff_rates, 1),
        
        exp24_lag1   = lag(exp_24, 1),
        exp12_lag1   = lag(exp_12, 1),
        infl_lag1    = lag(infl, 1),
        tc_yoy_lag1  = lag(tc_yoy, 1),
        wti_yoy_lag1 = lag(wti_yoy, 1),
        inf_usa_lag1 = lag(inf_usa, 1),
        imae_yoy_lag1 = lag(imae_yoy, 1)
    )

# Filtro HP para el log del IMAE
#hp_imae <- mFilter::hpfilter(df$l_imae, freq = 14400, type = "lambda")

#df$imae_gap_hp <- 100 * as.numeric(hp_imae$cycle)

# ------------------------------------------------------------
# 9) Base de modelado
# ------------------------------------------------------------
# Especificación inicial del anclaje:
# exp_24 ~ exp_12 + infl + tc_yoy + wti_yoy 
df_model <- df %>%
    select(
        fecha,
        exp_24,
        exp_12,
        infl,
        tc_yoy,
        inf_usa,
        wti_yoy,
        imae_yoy,
        abs_gap_exp24
          ) %>%
    filter(
        !is.na(exp_24),
        !is.na(exp_12),
        !is.na(infl),
        !is.na(tc_yoy),
        !is.na(wti_yoy)
    )


write_xlsx(df_model, "output/df_model_anclaje.xlsx")

# ------------------------------------------------------------
# 10) Modelo base estático
# ------------------------------------------------------------
modelo_base <- lm(
    exp_24 ~ exp_12 + infl + tc_yoy  + wti_yoy,
    data = df_model
)

summary(modelo_base)

cat("\n================ MODELO BASE ================\n")
print(summary(modelo_base))

writeLines(
    capture.output(summary(modelo_base)),
    "output/resumen_modelo_base_anclaje.txt"
)

cor(df_model[, c("exp_12", "infl", "tc_yoy", "wti_yoy", "imae_yoy")],
    use = "complete.obs")

# ------------------------------------------------------------
# 11) Función rolling
# ------------------------------------------------------------
rolling_fun <- function(data_window) {
    
    mod <- lm(
        exp_24 ~ exp_12 + infl + tc_yoy + wti_yoy,
        data = data_window
    )
    
    coefs <- broom::tidy(mod)
    
    tibble(
        beta_exp12      = coefs %>% filter(term == "exp_12") %>% pull(estimate),
        beta_inflacion  = coefs %>% filter(term == "infl") %>% pull(estimate),
        beta_tc         = coefs %>% filter(term == "tc_yoy") %>% pull(estimate),
        beta_wti        = coefs %>% filter(term == "wti_yoy") %>% pull(estimate),
        r2              = summary(mod)$r.squared,
        adj_r2          = summary(mod)$adj.r.squared,
        sigma           = summary(mod)$sigma
    )
}

# ------------------------------------------------------------
# 12) Estimación rolling
# ------------------------------------------------------------
rolling_results <- map_dfr(window_size:nrow(df_model), function(i) {
    
    data_window <- df_model[(i - window_size + 1):i, ]
    res <- rolling_fun(data_window)
    
    tibble(
        fecha_fin = max(data_window$fecha)
    ) %>%
        bind_cols(res)
})

write_xlsx(rolling_results, "output/rolling_results_anclaje.xlsx")

# ------------------------------------------------------------
# 13) Índice preliminar de anclaje
# ------------------------------------------------------------
# Mayor sensibilidad = menor anclaje
rolling_results <- rolling_results %>%
    mutate(
        z_beta_exp12     = as.numeric(scale(abs(beta_exp12))),
        z_beta_inflacion = as.numeric(scale(abs(beta_inflacion))),
        z_beta_tc        = as.numeric(scale(abs(beta_tc))),
        z_wti            = as.numeric(scale(abs(beta_wti))),
        indice_anclaje = -(
            z_beta_exp12 +
                z_beta_inflacion +
                z_beta_tc 
                
        ) / 4
    )

write_xlsx(
    rolling_results,
    "output/rolling_results_con_indice_anclaje.xlsx"
)

rolling_results <- rolling_results %>%
    mutate(
        indice_anclaje_100 = 100 * (
            indice_anclaje - min(indice_anclaje, na.rm = TRUE)
        ) / (
            max(indice_anclaje, na.rm = TRUE) - min(indice_anclaje, na.rm = TRUE)
        )
    )

# ------------------------------------------------------------
# 14) Desalineamiento respecto de la meta
# ------------------------------------------------------------
rolling_gap <- df_model %>%
    mutate(
        abs_gap_roll_12 = zoo::rollmean(abs_gap_exp24, k = 12, fill = NA, align = "right")
    ) %>%
    select(fecha, abs_gap_exp24, abs_gap_roll_12)

write_xlsx(rolling_gap, "output/rolling_gap_exp24_meta.xlsx")

# ------------------------------------------------------------
# 15) Gráficos
# ------------------------------------------------------------

# Índice de anclaje
g1 <- ggplot(rolling_results, aes(x = fecha_fin, y = indice_anclaje)) +
    geom_line(linewidth = 1) +
    labs(
        title = "Índice preliminar de anclaje de expectativas",
        subtitle = paste0("Ventana móvil de ", window_size, " meses"),
        x = "Fecha",
        y = "Índice (mayor valor = mejor anclaje)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold")
    )

g1
ggsave("output/graficos/indice_anclaje.png", g1, width = 10, height = 6, dpi = 300)

# Sensibilidad a exp_12
g2 <- ggplot(rolling_results, aes(x = fecha_fin, y = beta_exp12)) +
    geom_line(linewidth = 1) +
    labs(
        title = "Sensibilidad rolling de exp_24 a exp_12",
        subtitle = paste0("Ventana móvil de ", window_size, " meses"),
        x = "Fecha",
        y = "Coeficiente"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold")
    )
g2
ggsave("output/graficos/beta_exp12.png", g2, width = 10, height = 6, dpi = 300)

# Sensibilidad a inflación observada
g3 <- ggplot(rolling_results, aes(x = fecha_fin, y = beta_inflacion)) +
    geom_line(linewidth = 1) +
    labs(
        title = "Sensibilidad rolling de exp_24 a inflación observada",
        subtitle = paste0("Ventana móvil de ", window_size, " meses"),
        x = "Fecha",
        y = "Coeficiente"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold")
    )

g3
ggsave("output/graficos/beta_inflacion.png", g3, width = 10, height = 6, dpi = 300)

# Sensibilidad a tipo de cambio
g4 <- ggplot(rolling_results, aes(x = fecha_fin, y = beta_tc)) +
    geom_line(linewidth = 1) +
    labs(
        title = "Sensibilidad rolling de exp_24 a depreciación cambiaria",
        subtitle = paste0("Ventana móvil de ", window_size, " meses"),
        x = "Fecha",
        y = "Coeficiente"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold")
    )
g4
ggsave("output/graficos/beta_tc.png", g4, width = 10, height = 6, dpi = 300)

# Desviación respecto de la meta
g5 <- ggplot(rolling_gap, aes(x = fecha, y = abs_gap_roll_12)) +
    geom_line(linewidth = 1) +
    labs(
        title = "Desviación absoluta de exp_24 respecto de la meta",
        subtitle = "Promedio móvil de 12 meses",
        x = "Fecha",
        y = "|exp_24 - meta|"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold")
    )
g5
ggsave("output/graficos/abs_gap_exp24_meta.png", g5, width = 10, height = 6, dpi = 300)

# ------------------------------------------------------------
# 16) Resumen interpretativo
# ------------------------------------------------------------
resumen_txt <- c(
    "FASE 1: MEDICION PRELIMINAR DEL ANCLAJE DE EXPECTATIVAS",
    "-------------------------------------------------------",
    paste("Observaciones utilizadas:", nrow(df_model)),
    paste("Ventana rolling:", window_size, "meses"),
    "",
    "Interpretacion economica:",
    "- Menor coeficiente de exp_12 sobre exp_24 sugiere menor pass-through de expectativas de corto a mediano plazo.",
    "- Menor coeficiente de inflacion observada sugiere menor sensibilidad adaptativa.",
    "- Menor coeficiente del tipo de cambio sugiere menor contaminacion del horizonte de 24 meses por shocks cambiarios.",
    "- El indice_anclaje aumenta cuando cae la sensibilidad de exp_24 a esas variables.",
    "",
    "Proximo paso sugerido:",
    "- Definir periodos de anclaje debil / fuerte.",
    "- Estimar un VAR pequeno para comparar persistencia e impulso-respuesta."
)

writeLines(resumen_txt, "output/resumen_fase1_anclaje.txt")

cat("\nProceso completado con éxito.\n")
cat("Revisa la carpeta 'output' para tablas, resúmenes y gráficos.\n")