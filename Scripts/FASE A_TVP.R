# ============================================================
# TESIS - ANCLAJE DE EXPECTATIVAS DE INFLACIÓN EN GUATEMALA
# FASE A: Medida base de anclaje
# FASE B: Descomposición en shocks externos y domésticos
# ============================================================

# ------------------------------------------------------------
# 0) Limpieza de entorno
# ------------------------------------------------------------
rm(list = ls())
gc()
cat("\014")
if (!is.null(dev.list())) dev.off()

# ------------------------------------------------------------
# 1) Paquetes
# ------------------------------------------------------------
paquetes <- c(
    "readr",
    "dplyr",
    "ggplot2",
    "lubridate",
    "tidyr",
    "purrr",
    "broom",
    "scales",
    "zoo",
    "KFAS"
)

instalados <- rownames(installed.packages())
for (p in paquetes) {
    if (!(p %in% instalados)) install.packages(p, dependencies = TRUE)
}
invisible(lapply(paquetes, library, character.only = TRUE))

# Para evitar conflictos
lag <- dplyr::lag

# ------------------------------------------------------------
# 2) Carga de datos
# ------------------------------------------------------------
# Ajusta esta ruta según donde tengas guardado el archivo
ruta_datos <- "bd_sample_modelo.csv"

df_raw <- read_csv(ruta_datos, show_col_types = FALSE)

# ------------------------------------------------------------
# 3) Preparación general de la base
# ------------------------------------------------------------
df <- df_raw %>%
    mutate(
        fecha = as.Date(fecha)
    ) %>%
    arrange(fecha)

# Revisión rápida
glimpse(df)
summary(df)

# ============================================================
# FASE A
# MEDIDA BASE DE ANCLAJE
# ============================================================

# ------------------------------------------------------------
# Intuición econométrica
# ------------------------------------------------------------
# Se quiere medir si los cambios en la expectativa a 24 meses
# reaccionan a un shock inflacionario.
#
# Ecuación base:
#   d_exp24_t = alpha_t + beta_t * shock_inf_t + error_t
#
# donde:
#   d_exp24_t   = exp_inf_24m_t - exp_inf_24m_(t-1)
#   shock_inf_t = infl_gt_t - exp_inf_12m_(t-1)
#
# Interpretación:
# - beta_t cercano a 0  -> expectativas más ancladas
# - beta_t alto/positivo -> expectativas más sensibles, menor anclaje
#
# alpha_t y beta_t evolucionan con el tiempo vía random walk.

# ------------------------------------------------------------
# 4) Construcción de variables FASE A
# ------------------------------------------------------------
df_faseA <- df %>%
    mutate(
        exp12_lag1 = lag(exp_inf_12m, 1),
        d_exp24    = exp_inf_24m - lag(exp_inf_24m, 1),
        shock_inf  = infl_gt - exp12_lag1
    ) %>%
    filter(
        !is.na(d_exp24),
        !is.na(shock_inf)
    )

# Estadísticos descriptivos básicos
cat("\n=============================\n")
cat("FASE A - Estadísticos básicos\n")
cat("=============================\n")
summary(select(df_faseA, d_exp24, shock_inf))

# ------------------------------------------------------------
# 5) Benchmark estático: OLS
# ------------------------------------------------------------
modelo_ols_A <- lm(d_exp24 ~ shock_inf, data = df_faseA)

cat("\n=============================\n")
cat("FASE A - OLS benchmark\n")
cat("=============================\n")
print(summary(modelo_ols_A))

# ------------------------------------------------------------
# 6) Modelo TVP con filtro de Kalman (KFAS)
# ------------------------------------------------------------
# Especificación:
#   y_t = mu_t + beta_t * shock_inf_t + eps_t
#   mu_t   = mu_(t-1) + xi_t
#   beta_t = beta_(t-1) + eta_t
#
# En KFAS:
# - SSMtrend(degree = 1) modela un nivel local (mu_t)
# - SSMregression(..., Q = NA) hace el coeficiente TVP

y_A <- df_faseA$d_exp24
x_A <- df_faseA$shock_inf

modelo_ssm_A <- SSModel(
    y_A ~
        SSMtrend(degree = 1, Q = list(NA)) + #nivel local que cambia en el tiempo
        SSMregression(~ x_A, Q = matrix(NA)), #shock como regresor, pero time-varying.
    H = NA
)

# Valores iniciales para optimización:
# usamos log-varianzas para garantizar positividad
inits_A <- log(c(var(y_A, na.rm = TRUE), 0.01, 0.01))

fit_A <- fitSSM( #calibra cuánto se permite variar en el tiempo el intercepto y cuánto se permite variar el coeficiente del shock.
    inits = inits_A,
    model = modelo_ssm_A,
    method = "BFGS" #<- máximo de la log-verosimilitud.
)

modelo_kfas_A <- fit_A$model

# Suavizado de estados
kfs_A <- KFS(modelo_kfas_A, smoothing = c("state", "mean")) #aplica el filtro de Kalman y, sobre todo aquí, el suavizador de Kalman al modelo ya estimado

# Extraer estados suavizados
# En este modelo:
# estado 1 = alpha_t (nivel local)
# estado 2 = beta_t  (coeficiente TVP del shock)
estados_A <- as.data.frame(kfs_A$alphahat)
colnames(estados_A) <- c("alpha_t", "beta_t")

df_result_A <- bind_cols(
    df_faseA %>% select(fecha, d_exp24, shock_inf),
    estados_A
)

# ------------------------------------------------------------
# 7) Visualización FASE A
# ------------------------------------------------------------
graf_beta_A <- ggplot(df_result_A, aes(x = fecha, y = beta_t)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "FASE A: Sensibilidad TVP de expectativas a shock inflacionario",
        subtitle = "beta_t cercano a cero sugiere mejor anclaje",
        x = NULL,
        y = expression(beta[t])
    ) +
    theme_minimal(base_size = 13)

graf_alpha_A <- ggplot(df_result_A, aes(x = fecha, y = alpha_t)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "FASE A: Componente autónomo de revisión de expectativas",
        x = NULL,
        y = expression(alpha[t])
    ) +
    theme_minimal(base_size = 13)

print(graf_beta_A)
print(graf_alpha_A)

# ------------------------------------------------------------
# 8) Indicador simple de anclaje
# ------------------------------------------------------------
# Para interpretación práctica:
# mayor beta_t  => menor anclaje
# una transformación simple es:
#   indice_anclaje = -beta_t
#
# Si quieres, luego lo podemos normalizar a [0,1] o [0,100].

df_result_A <- df_result_A %>%
    mutate(
        indice_anclaje_A = -beta_t
    )

graf_indice_A <- ggplot(df_result_A, aes(x = fecha, y = indice_anclaje_A)) +
    geom_line(linewidth = 0.9) +
    labs(
        title = "FASE A: Índice simple de anclaje",
        subtitle = "Valores mayores indican mejor anclaje (porque beta_t es menor)",
        x = NULL,
        y = "Índice de anclaje"
    ) +
    theme_minimal(base_size = 13)

print(graf_indice_A)

# ------------------------------------------------------------
# 9) Guardar resultados FASE A
# ------------------------------------------------------------
write_csv(df_result_A, "resultados_fase_A_kalman.csv")

ggsave("faseA_beta_t.png", graf_beta_A, width = 10, height = 5, dpi = 300)
ggsave("faseA_alpha_t.png", graf_alpha_A, width = 10, height = 5, dpi = 300)
ggsave("faseA_indice_anclaje.png", graf_indice_A, width = 10, height = 5, dpi = 300)

# ============================================================
# CHEQUEO A: VARIANZAS ESTIMADAS DEL MODELO TVP
# ============================================================

cat("\n====================================================\n")
cat("CHEQUEO A: VARIANZAS ESTIMADAS DEL MODELO TVP\n")
cat("====================================================\n")

# Extraer varianzas estimadas
H_est <- modelo_kfas_A$H[1, 1, 1]
Q_alpha_est <- modelo_kfas_A$Q[1, 1, 1]
Q_beta_est  <- modelo_kfas_A$Q[2, 2, 1]

cat("H       (error observacional) = ", format(H_est, scientific = TRUE), "\n", sep = "")
cat("Q_alpha (estado alpha_t)      = ", format(Q_alpha_est, scientific = TRUE), "\n", sep = "")
cat("Q_beta  (estado beta_t)       = ", format(Q_beta_est, scientific = TRUE), "\n", sep = "")

cat("\nRazones diagnósticas:\n")
cat("Q_alpha / H = ", format(Q_alpha_est / H_est, scientific = TRUE), "\n", sep = "")
cat("Q_beta  / H = ", format(Q_beta_est  / H_est, scientific = TRUE), "\n", sep = "")

if (Q_beta_est < 1e-6) {
  cat("\nALERTA: Q_beta es extremadamente pequeño.\n")
  cat("El modelo impone una trayectoria muy rígida para beta_t.\n")
} else if (Q_beta_est < 1e-5) {
  cat("\nADVERTENCIA: Q_beta es muy pequeño.\n")
  cat("La trayectoria de beta_t probablemente será muy suave y con cambios graduales.\n")
} else {
  cat("\nQ_beta no es particularmente pequeño.\n")
  cat("La suavidad de beta_t parecería venir más de los datos que de la restricción del modelo.\n")
}

# ============================================================
# CHEQUEO B: BANDAS DE CONFIANZA PARA beta_t
# ============================================================

# kfs_A$V es un arreglo 3D:
# dimensión 1 = estado
# dimensión 2 = estado
# dimensión 3 = tiempo

# Varianza suavizada de alpha_t en cada t
var_alpha_t <- sapply(1:dim(kfs_A$V)[3], function(i) kfs_A$V[1, 1, i])

# Varianza suavizada de beta_t en cada t
var_beta_t  <- sapply(1:dim(kfs_A$V)[3], function(i) kfs_A$V[2, 2, i])

# Errores estándar
se_alpha_t <- sqrt(var_alpha_t)
se_beta_t  <- sqrt(var_beta_t)

# Agregar al data frame de resultados
df_result_A <- df_result_A %>%
  mutate(
    se_alpha_t  = se_alpha_t,
    se_beta_t   = se_beta_t,
    beta_low68  = beta_t - 1.00 * se_beta_t,
    beta_up68   = beta_t + 1.00 * se_beta_t,
    alpha_low68 = alpha_t - 1.00 * se_alpha_t,
    alpha_up68  = alpha_t + 1.00 * se_alpha_t
  )

# Gráfico con banda de confianza para beta_t
graf_beta_ci_A <- ggplot(df_result_A, aes(x = fecha, y = beta_t)) +
  geom_ribbon(aes(ymin = beta_low68, ymax = beta_up68), alpha = 0.20) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "FASE A: beta_t con banda de confianza al 68%",
    subtitle = "Sensibilidad TVP de expectativas a shock inflacionario",
    x = NULL,
    y = expression(beta[t])
  ) +
  theme_minimal(base_size = 13)

print(graf_beta_ci_A)

# Diagnóstico automático simple:
df_result_A <- df_result_A %>%
    mutate(
        beta_significativo_95 = ifelse(beta_low95 > 0 | beta_up95 < 0, 1, 0)
    )

cat("\n====================================================\n")
cat("CHEQUEO B: SIGNIFICANCIA DE beta_t\n")
cat("====================================================\n")
cat("Porcentaje de periodos en que beta_t es significativamente distinto de 0 (95%): ",
    round(mean(df_result_A$beta_significativo_95, na.rm = TRUE) * 100, 2), "%\n", sep = "")

# Guardar gráfico
ggsave("faseA_beta_t_con_bandas.png", graf_beta_ci_A, width = 10, height = 5, dpi = 300)

# ============================================================
# CHEQUEO C: COMPARACIÓN CON ROLLING REGRESSION
# ============================================================

library(purrr)
library(broom)

# Tamaño de ventana rolling
window_size <- 60

# Función para estimar beta en cada ventana
rolling_fun <- function(data_window) {
    modelo <- lm(d_exp24 ~ shock_inf, data = data_window)
    beta_row <- broom::tidy(modelo) %>%
        dplyr::filter(term == "shock_inf")
    
    tibble(
        estimate  = beta_row$estimate,
        std.error = beta_row$std.error,
        statistic = beta_row$statistic,
        p.value   = beta_row$p.value
    )
}

# Aplicar rolling
rolling_results <- map_dfr(window_size:nrow(df_faseA), function(i) {
    
    data_window <- df_faseA[(i - window_size + 1):i, ]
    res <- rolling_fun(data_window)
    
    tibble(
        fecha_fin  = max(data_window$fecha),
        beta_roll  = res$estimate,
        se_roll    = res$std.error,
        beta_low95 = res$estimate - 1.96 * res$std.error,
        beta_up95  = res$estimate + 1.96 * res$std.error
    )
})

# Gráfico rolling beta
graf_beta_roll <- ggplot(rolling_results, aes(x = fecha_fin, y = beta_roll)) +
    geom_ribbon(aes(ymin = beta_low95, ymax = beta_up95), alpha = 0.20) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "Rolling regression: sensibilidad de expectativas al shock",
        subtitle = paste0("Ventana móvil de ", window_size, " meses"),
        x = NULL,
        y = expression(hat(beta)[roll])
    ) +
    theme_minimal(base_size = 13)

print(graf_beta_roll)

# ------------------------------------------------------------
# Comparación gráfica TVP vs rolling
# ------------------------------------------------------------

comparacion_beta <- df_result_A %>%
    select(fecha, beta_t) %>%
    rename(fecha_fin = fecha) %>%
    inner_join(rolling_results %>% select(fecha_fin, beta_roll), by = "fecha_fin") %>%
    pivot_longer(
        cols = c(beta_t, beta_roll),
        names_to = "metodo",
        values_to = "beta"
    ) %>%
    mutate(
        metodo = dplyr::recode(
            metodo,
            beta_t = "Kalman TVP",
            beta_roll = "Rolling OLS"
        )
    )

graf_comp_beta <- ggplot(comparacion_beta, aes(x = fecha_fin, y = beta, linetype = metodo)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "Comparación: beta_t Kalman vs rolling OLS",
        subtitle = "Si ambas trayectorias son parecidas, el hallazgo gana robustez",
        x = NULL,
        y = expression(beta[t]),
        linetype = NULL
    ) +
    theme_minimal(base_size = 13)

print(graf_comp_beta)

# ------------------------------------------------------------
# Correlación entre ambas series
# ------------------------------------------------------------
comp_wide <- df_result_A %>%
    select(fecha, beta_t) %>%
    rename(fecha_fin = fecha) %>%
    inner_join(rolling_results %>% select(fecha_fin, beta_roll), by = "fecha_fin")

cor_beta <- cor(comp_wide$beta_t, comp_wide$beta_roll, use = "complete.obs")

cat("\n====================================================\n")
cat("CHEQUEO C: COMPARACIÓN TVP vs ROLLING\n")
cat("====================================================\n")
cat("Correlación entre beta_t (Kalman) y beta_roll (rolling OLS): ",
    round(cor_beta, 4), "\n", sep = "")

# Guardar outputs
write_csv(rolling_results, "rolling_results_fase_A.csv")
ggsave("faseA_beta_rolling.png", graf_beta_roll, width = 10, height = 5, dpi = 300)
ggsave("faseA_comparacion_kalman_vs_rolling.png", graf_comp_beta, width = 10, height = 5, dpi = 300)



# ============================================================
# FASE B
# DESCOMPOSICIÓN EN SHOCKS EXTERNOS Y DOMÉSTICOS
# ============================================================

# ------------------------------------------------------------
# Intuición econométrica
# ------------------------------------------------------------
# En vez de usar un solo shock agregado, construimos dos bloques:
#
# 1) Shock externo:
#    - d_effr
#    - d_infl_us
#    - oil_var
#
# 2) Shock doméstico:
#    - infl_gt_lag1
#    - d_tasa_pol_lag3
#    - imae_lag2
#    - tcdep_var
#    - deficit_fiscal
#
# Luego sintetizamos cada bloque en un factor principal (PCA).
# Esto evita meter demasiados coeficientes TVP al mismo tiempo.
#
# Modelo:
#   d_exp24_t = alpha_t + b_ext_t * shock_ext_t + b_dom_t * shock_dom_t + e_t
#
# Interpretación:
# - b_ext_t alto  -> expectativas reaccionan más a shocks externos
# - b_dom_t alto  -> expectativas reaccionan más a shocks domésticos

# ------------------------------------------------------------
# 10) Construcción de variables FASE B
# ------------------------------------------------------------
df_faseB <- df %>%
    mutate(
        d_exp24         = exp_inf_24m - lag(exp_inf_24m, 1),
        d_tasa_pol_lag3 = lag(d_tasa_pol, 3)
    ) %>%
    filter(
        !is.na(d_exp24),
        !is.na(d_effr),
        !is.na(d_infl_us),
        !is.na(oil_var),
        !is.na(infl_gt_lag1),
        !is.na(d_tasa_pol_lag3),
        !is.na(imae_lag2),
        !is.na(tcdep_var),
        !is.na(deficit_fiscal)
    )

# ------------------------------------------------------------
# 11) Construcción de shocks externos y domésticos con PCA
# ------------------------------------------------------------
# Estandarizamos antes de PCA para evitar problemas de escala

# Bloque externo
mat_ext <- df_faseB %>%
    select(d_effr, d_infl_us, oil_var) %>%
    scale(center = TRUE, scale = TRUE)

pca_ext <- prcomp(mat_ext, center = FALSE, scale. = FALSE)

# Bloque doméstico
mat_dom <- df_faseB %>%
    select(infl_gt_lag1, d_tasa_pol_lag3, imae_lag2, tcdep_var, deficit_fiscal) %>%
    scale(center = TRUE, scale = TRUE)

pca_dom <- prcomp(mat_dom, center = FALSE, scale. = FALSE)

df_faseB <- df_faseB %>%
    mutate(
        shock_externo  = as.numeric(pca_ext$x[, 1]),
        shock_domestico = as.numeric(pca_dom$x[, 1])
    )

cat("\n=========================================\n")
cat("FASE B - Varianza explicada por PCA externo\n")
cat("=========================================\n")
print(summary(pca_ext))

cat("\n===========================================\n")
cat("FASE B - Varianza explicada por PCA doméstico\n")
cat("===========================================\n")
print(summary(pca_dom))

# ------------------------------------------------------------
# 12) Benchmark estático FASE B
# ------------------------------------------------------------
modelo_ols_B <- lm(
    d_exp24 ~ shock_externo + shock_domestico,
    data = df_faseB
)

cat("\n=============================\n")
cat("FASE B - OLS benchmark\n")
cat("=============================\n")
print(summary(modelo_ols_B))

# ------------------------------------------------------------
# 13) Modelo TVP con dos coeficientes variables
# ------------------------------------------------------------
y_B     <- df_faseB$d_exp24
x_ext_B <- df_faseB$shock_externo
x_dom_B <- df_faseB$shock_domestico

modelo_ssm_B <- SSModel(
    y_B ~
        SSMtrend(degree = 1, Q = list(NA)) +
        SSMregression(~ x_ext_B + x_dom_B, Q = diag(NA, 2)),
    H = NA
)

# Parámetros:
# H + Q_nivel + Q_beta_ext + Q_beta_dom = 4 parámetros
inits_B <- log(c(var(y_B, na.rm = TRUE), 0.01, 0.01, 0.01))

fit_B <- fitSSM(
    inits = inits_B,
    model = modelo_ssm_B,
    method = "BFGS"
)

modelo_kfas_B <- fit_B$model
kfs_B <- KFS(modelo_kfas_B, smoothing = c("state", "mean"))

estados_B <- as.data.frame(kfs_B$alphahat)
colnames(estados_B) <- c("alpha_t", "beta_ext_t", "beta_dom_t")

df_result_B <- bind_cols(
    df_faseB %>% select(fecha, d_exp24, shock_externo, shock_domestico),
    estados_B
)

# ------------------------------------------------------------
# 14) Visualización FASE B
# ------------------------------------------------------------
graf_beta_ext <- ggplot(df_result_B, aes(x = fecha, y = beta_ext_t)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "FASE B: Sensibilidad TVP a shocks externos",
        x = NULL,
        y = expression(beta[ext][t])
    ) +
    theme_minimal(base_size = 13)

graf_beta_dom <- ggplot(df_result_B, aes(x = fecha, y = beta_dom_t)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "FASE B: Sensibilidad TVP a shocks domésticos",
        x = NULL,
        y = expression(beta[dom][t])
    ) +
    theme_minimal(base_size = 13)

print(graf_beta_ext)
print(graf_beta_dom)

round(pca_ext$rotation, 3)
round(pca_dom$rotation, 3)

# Serie comparativa
df_long_B <- df_result_B %>%
    select(fecha, beta_ext_t, beta_dom_t) %>%
    pivot_longer(
        cols = c(beta_ext_t, beta_dom_t),
        names_to = "tipo_shock",
        values_to = "beta_t"
    ) %>%
    mutate(
        tipo_shock = recode(
            tipo_shock,
            beta_ext_t = "Shock externo",
            beta_dom_t = "Shock doméstico"
        )
    )

graf_comparativo_B <- ggplot(df_long_B, aes(x = fecha, y = beta_t, linetype = tipo_shock)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "FASE B: Comparación de sensibilidades TVP",
        x = NULL,
        y = expression(beta[t]),
        linetype = NULL
    ) +
    theme_minimal(base_size = 13)

print(graf_comparativo_B)

# ------------------------------------------------------------
# 15) Índices interpretativos FASE B
# ------------------------------------------------------------
df_result_B <- df_result_B %>%
    mutate(
        indice_anclaje_externo  = -beta_ext_t,
        indice_anclaje_domestico = -beta_dom_t
    )

# ------------------------------------------------------------
# 16) Guardar resultados FASE B
# ------------------------------------------------------------
write_csv(df_result_B, "resultados_fase_B_kalman.csv")

ggsave("faseB_beta_externo.png", graf_beta_ext, width = 10, height = 5, dpi = 300)
ggsave("faseB_beta_domestico.png", graf_beta_dom, width = 10, height = 5, dpi = 300)
ggsave("faseB_comparativo.png", graf_comparativo_B, width = 10, height = 5, dpi = 300)

# ============================================================
# 17) Resumen interpretativo automático
# ============================================================

cat("\n====================================================\n")
cat("RESUMEN INTERPRETATIVO AUTOMÁTICO\n")
cat("====================================================\n")

cat("\nFASE A:\n")
cat("Promedio beta_t =", round(mean(df_result_A$beta_t, na.rm = TRUE), 4), "\n")
cat("Desv. estándar beta_t =", round(sd(df_result_A$beta_t, na.rm = TRUE), 4), "\n")
cat("Último beta_t =", round(tail(df_result_A$beta_t, 1), 4), "\n")

cat("\nFASE B:\n")
cat("Promedio beta_ext_t =", round(mean(df_result_B$beta_ext_t, na.rm = TRUE), 4), "\n")
cat("Promedio beta_dom_t =", round(mean(df_result_B$beta_dom_t, na.rm = TRUE), 4), "\n")
cat("Último beta_ext_t =", round(tail(df_result_B$beta_ext_t, 1), 4), "\n")
cat("Último beta_dom_t =", round(tail(df_result_B$beta_dom_t, 1), 4), "\n")

# ============================================================
# FIN
# ============================================================