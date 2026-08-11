# ============================================================
# 01. ROLLING REGRESSIONS — MODELO DE STROHSAL
# Medición dinámica del anclaje de las expectativas de inflación
# ============================================================

# Ecuación de formación de expectativas:
#
# e24_t = (1 - theta1 - theta2) * meta_t
#         + theta1 * inflacion_(t-1)
#         + theta2 * e12_(t-1)
#
# Forma estimable:
#
# [Delta e24_t - (meta_t - e24_(t-1))]
#     = theta1 * [inflacion_(t-1) - meta_t]
#       + theta2 * [e12_(t-1) - meta_t]
#       + u_t
#
# Indicador de anclaje:
#
# A = 1 - theta1 - theta2
#
# Interpretación:
# - theta1: sensibilidad frente a la inflación observada rezagada.
# - theta2: sensibilidad frente a las expectativas a 12 meses.
# - A: peso de la meta de inflación.
# ============================================================

# ------------------------------------------------------------
# Paquetes
# ------------------------------------------------------------
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(broom)
library(purrr)
library(lmtest)
library(sandwich)
library(tidyr)

# ------------------------------------------------------------
# 0) Limpieza de entorno
# ------------------------------------------------------------
rm(list = ls())
cat("\014")
if (!is.null(dev.list())) dev.off()

# ------------------------------------------------------------
# 1) Parámetros generales
# ------------------------------------------------------------
window_size <- 60   # Ventana principal: 36 meses
nw_lag      <- 3    # Rezagos para errores Newey-West HAC

# Prueba de robustez:
# window_size <- 30

# ------------------------------------------------------------
# 2) Cargar base y construir variables
# ------------------------------------------------------------
df_model <- read_csv(
    "bd_sample_modelo.csv",
    show_col_types = FALSE
) |>
    mutate(
        fecha = ymd(fecha)
    ) |>
    arrange(fecha) |>
    mutate(
        # Meta central vigente
        meta_inflacion = case_when(
            fecha < ymd("2012-01-01") ~ 5.0,
            fecha < ymd("2013-01-01") ~ 4.5,
            TRUE                       ~ 4.0
        ),

        # Rezagos de las expectativas
        exp12_lag1 = lag(exp_inf_12m, 1),
        exp24_lag1_modelo = lag(exp_inf_24m, 1),

        # Cambio mensual de la expectativa a 24 meses
        d_exp24 = exp_inf_24m - exp24_lag1_modelo,

        # Variable dependiente de la forma estimable
        y_strohsal =
            d_exp24 -
            (meta_inflacion - exp24_lag1_modelo),

        # Fuente de desanclaje 1:
        # inflación observada rezagada respecto de la meta vigente
        x_inflacion =
            infl_gt_lag1 - meta_inflacion,

        # Fuente de desanclaje 2:
        # expectativa a 12 meses rezagada respecto de la meta vigente
        x_exp12 =
            exp12_lag1 - meta_inflacion
    ) |>
    filter(
        complete.cases(
            y_strohsal,
            x_inflacion,
            x_exp12
        )
    )

cat(
    "\nPrimera fecha efectiva: ",
    as.character(min(df_model$fecha)),
    "\n"
)

cat(
    "Última fecha efectiva: ",
    as.character(max(df_model$fecha)),
    "\n"
)

cat(
    "Observaciones efectivas: ",
    nrow(df_model),
    "\n"
)

if (nrow(df_model) < window_size) {
    stop(
        "La muestra efectiva tiene menos observaciones que la ventana móvil."
    )
}

# ------------------------------------------------------------
# 3) Función de estimación para cada ventana
# ------------------------------------------------------------
rolling_fun <- function(data_window) {

    # La regresión se estima sin constante.
    # La ausencia de intercepto incorpora la restricción teórica
    # de que los pesos de la meta, la inflación y las expectativas
    # de corto plazo suman uno.
    modelo <- lm(
        y_strohsal ~ 0 + x_inflacion + x_exp12,
        data = data_window
    )

    # Matriz de varianzas y covarianzas Newey-West HAC
    vcov_hac <- sandwich::NeweyWest(
        modelo,
        lag = nw_lag,
        prewhite = FALSE,
        adjust = TRUE
    )

    resultado_hac <- lmtest::coeftest(
        modelo,
        vcov. = vcov_hac
    )

    resultados <- broom::tidy(resultado_hac)
    ajuste     <- broom::glance(modelo)

    coef_theta1 <- resultados |>
        filter(term == "x_inflacion")

    coef_theta2 <- resultados |>
        filter(term == "x_exp12")

    if (nrow(coef_theta1) != 1 || nrow(coef_theta2) != 1) {
        stop(
            "No fue posible identificar de forma única los coeficientes ",
            "x_inflacion y x_exp12 en una de las ventanas."
        )
    }

    # Extraer los resultados como escalares antes de construir el tibble.
    theta1_est <- coef_theta1$estimate[[1]]
    theta1_se  <- coef_theta1$std.error[[1]]
    theta1_p   <- coef_theta1$p.value[[1]]

    theta2_est <- coef_theta2$estimate[[1]]
    theta2_se  <- coef_theta2$std.error[[1]]
    theta2_p   <- coef_theta2$p.value[[1]]

    # --------------------------------------------------------
    # Indicador de anclaje
    # A = 1 - theta1 - theta2
    # --------------------------------------------------------
    anclaje_est <- 1 - theta1_est - theta2_est

    # Error estándar de A mediante el método delta:
    # Var(A) = Var(theta1) + Var(theta2) + 2 Cov(theta1, theta2)
    var_anclaje <-
        vcov_hac["x_inflacion", "x_inflacion"] +
        vcov_hac["x_exp12", "x_exp12"] +
        2 * vcov_hac["x_inflacion", "x_exp12"]

    se_anclaje_est <- sqrt(max(var_anclaje, 0))

    # --------------------------------------------------------
    # Prueba conjunta de anclaje perfecto
    # H0: theta1 = 0 y theta2 = 0
    # --------------------------------------------------------
    theta_vector <- c(
        theta1_est,
        theta2_est
    )

    vcov_theta <- vcov_hac[
        c("x_inflacion", "x_exp12"),
        c("x_inflacion", "x_exp12")
    ]

    inversa_vcov <- tryCatch(
        solve(vcov_theta),
        error = function(e) NULL
    )

    if (is.null(inversa_vcov)) {
        wald_perfecto <- NA_real_
        p_anclaje_perfecto <- NA_real_
    } else {
        wald_perfecto <- as.numeric(
            t(theta_vector) %*%
                inversa_vcov %*%
                theta_vector
        )

        p_anclaje_perfecto <- pchisq(
            wald_perfecto,
            df = 2,
            lower.tail = FALSE
        )
    }

    tibble(
        # Sensibilidad frente a la inflación pasada
        theta1 = theta1_est,
        se_theta1 = theta1_se,
        p_theta1 = theta1_p,

        # Sensibilidad frente a las expectativas a 12 meses
        theta2 = theta2_est,
        se_theta2 = theta2_se,
        p_theta2 = theta2_p,

        # Indicador de anclaje
        anclaje = anclaje_est,
        se_anclaje = se_anclaje_est,

        # Prueba conjunta H0: theta1 = theta2 = 0
        wald_anclaje_perfecto = wald_perfecto,
        p_anclaje_perfecto = p_anclaje_perfecto,

        # Estadísticos del modelo
        r2_ajustado = ajuste$adj.r.squared,
        rmse = sqrt(mean(residuals(modelo)^2)),
        n = nobs(modelo)
    )
}

# ------------------------------------------------------------
# 4) Ciclo móvil
# ------------------------------------------------------------
rolling_results <- purrr::map_dfr(
    window_size:nrow(df_model),
    function(i) {

        data_window <- df_model[
            (i - window_size + 1):i,
        ]

        res <- rolling_fun(data_window)

        tibble(
            fecha_inicio = min(data_window$fecha),
            fecha_fin = max(data_window$fecha)
        ) |>
            bind_cols(res)
    }
)

# ------------------------------------------------------------
# 5) Intervalos de confianza y validación de restricciones
# ------------------------------------------------------------
rolling_results <- rolling_results |>
    mutate(
        # Intervalos de theta1
        theta1_ic_inf = theta1 - 1.96 * se_theta1,
        theta1_ic_sup = theta1 + 1.96 * se_theta1,

        # Intervalos de theta2
        theta2_ic_inf = theta2 - 1.96 * se_theta2,
        theta2_ic_sup = theta2 + 1.96 * se_theta2,

        # Intervalos del indicador de anclaje
        anclaje_ic_inf = anclaje - 1.96 * se_anclaje,
        anclaje_ic_sup = anclaje + 1.96 * se_anclaje,

        # Restricciones compatibles con la interpretación de ponderaciones
        restricciones_validas =
            theta1 >= 0 &
            theta2 >= 0 &
            anclaje >= 0 &
            anclaje <= 1
    )

# ------------------------------------------------------------
# 6) Guardar resultados
# ------------------------------------------------------------
if (!dir.exists("resultados_obj1")) {
    dir.create("resultados_obj1")
}

write_csv(
    rolling_results,
    paste0(
        "resultados_obj1/01_rolling_strohsal_",
        window_size,
        "m.csv"
    )
)

# ------------------------------------------------------------
# 7) Gráfico de theta1
# ------------------------------------------------------------
p_theta1 <- ggplot(
    rolling_results,
    aes(x = fecha_fin, y = theta1)
) +
    geom_ribbon(
        aes(
            ymin = theta1_ic_inf,
            ymax = theta1_ic_sup
        ),
        alpha = 0.2
    ) +
    geom_line(linewidth = 1) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    labs(
        title = paste0(
            "Sensibilidad frente a la inflación pasada ",
            "(ventana de ",
            window_size,
            " meses)"
        ),
        subtitle = expression(theta[1*w]),
        x = "Fecha final de la ventana",
        y = expression(theta[1*w])
    ) +
    theme_minimal()

ggsave(
    paste0(
        "resultados_obj1/01_theta1_strohsal_",
        window_size,
        "m.png"
    ),
    plot = p_theta1,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_theta1)

# ------------------------------------------------------------
# 8) Gráfico de theta2
# ------------------------------------------------------------
p_theta2 <- ggplot(
    rolling_results,
    aes(x = fecha_fin, y = theta2)
) +
    geom_ribbon(
        aes(
            ymin = theta2_ic_inf,
            ymax = theta2_ic_sup
        ),
        alpha = 0.2
    ) +
    geom_line(linewidth = 1) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    labs(
        title = paste0(
            "Sensibilidad frente a las expectativas a 12 meses ",
            "(ventana de ",
            window_size,
            " meses)"
        ),
        subtitle = expression(theta[2*w]),
        x = "Fecha final de la ventana",
        y = expression(theta[2*w])
    ) +
    theme_minimal()

ggsave(
    paste0(
        "resultados_obj1/01_theta2_strohsal_",
        window_size,
        "m.png"
    ),
    plot = p_theta2,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_theta2)

# ------------------------------------------------------------
# 9) Gráfico del indicador de anclaje
# ------------------------------------------------------------
p_anclaje <- ggplot(
    rolling_results,
    aes(x = fecha_fin, y = anclaje)
) +
    geom_ribbon(
        aes(
            ymin = anclaje_ic_inf,
            ymax = anclaje_ic_sup
        ),
        alpha = 0.2
    ) +
    geom_line(linewidth = 1) +
    geom_hline(
        yintercept = 1,
        linetype = "dotted"
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    labs(
        title = paste0(
            "Indicador de anclaje de las expectativas ",
            "(ventana de ",
            window_size,
            " meses)"
        ),
        subtitle = expression(
            A[w] == 1 - theta[1*w] - theta[2*w]
        ),
        x = "Fecha final de la ventana",
        y = expression(A[w])
    ) +
    theme_minimal()

ggsave(
    paste0(
        "resultados_obj1/01_anclaje_strohsal_",
        window_size,
        "m.png"
    ),
    plot = p_anclaje,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_anclaje)

# ------------------------------------------------------------
# 10) Gráfico conjunto de los tres pesos
# ------------------------------------------------------------
pesos_largos <- rolling_results |>
    select(
        fecha_fin,
        theta1,
        theta2,
        anclaje
    ) |>
    pivot_longer(
        cols = c(theta1, theta2, anclaje),
        names_to = "componente",
        values_to = "peso"
    ) |>
    mutate(
        componente = recode(
            componente,
            theta1 = "Inflación pasada",
            theta2 = "Expectativas a 12 meses",
            anclaje = "Meta de inflación"
        )
    )

p_pesos <- ggplot(
    pesos_largos,
    aes(
        x = fecha_fin,
        y = peso,
        linetype = componente
    )
) +
    geom_line(linewidth = 0.9) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    geom_hline(
        yintercept = 1,
        linetype = "dotted"
    ) +
    labs(
        title = paste0(
            "Composición de las expectativas a 24 meses ",
            "(ventana de ",
            window_size,
            " meses)"
        ),
        subtitle = paste0(
            "Pesos estimados de la meta, la inflación pasada ",
            "y las expectativas a 12 meses"
        ),
        x = "Fecha final de la ventana",
        y = "Peso estimado",
        linetype = "Componente"
    ) +
    theme_minimal()

ggsave(
    paste0(
        "resultados_obj1/01_pesos_strohsal_",
        window_size,
        "m.png"
    ),
    plot = p_pesos,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_pesos)

# ------------------------------------------------------------
# 11) Resumen en consola
# ------------------------------------------------------------
cat(
    "\nVentana móvil: ",
    window_size,
    " meses\n"
)

cat(
    "Número de estimaciones: ",
    nrow(rolling_results),
    "\n"
)

cat(
    "Primera estimación: ",
    as.character(min(rolling_results$fecha_fin)),
    "\n"
)

cat(
    "Última estimación: ",
    as.character(max(rolling_results$fecha_fin)),
    "\n"
)

cat(
    "Ventanas con restricciones válidas: ",
    sum(rolling_results$restricciones_validas, na.rm = TRUE),
    " de ",
    nrow(rolling_results),
    "\n"
)

cat(
    "\nPromedio del indicador de anclaje: ",
    round(mean(rolling_results$anclaje, na.rm = TRUE), 3),
    "\n"
)

cat(
    "\nArchivos generados en la carpeta resultados_obj1:\n",
    "- Resultados completos en CSV\n",
    "- Gráfico de theta1\n",
    "- Gráfico de theta2\n",
    "- Gráfico del indicador de anclaje\n",
    "- Gráfico conjunto de los tres pesos\n"
)

