# ============================================================
# 01. ROLLING REGRESSIONS
# Sensibilidad de las expectativas a 24 meses frente a la inflación
# Modelo ampliado: persistencia, brecha de inflación y tipo de cambio
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

# ------------------------------------------------------------
# 0) Limpieza de entorno
# ------------------------------------------------------------
rm(list = ls())
cat("\014")
if (!is.null(dev.list())) dev.off()

# ------------------------------------------------------------
# 1) Cargar base y construir variables del modelo
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
        # Meta central vigente:
        # 5.0% durante 2010-2011
        # 4.5% durante 2012
        # 4.0% desde enero de 2013
        meta_inflacion = case_when(
            fecha < ymd("2012-01-01") ~ 5.0,
            fecha < ymd("2013-01-01") ~ 4.5,
            TRUE                       ~ 4.0
        ),

        # Meta correspondiente al período t-1
        meta_inflacion_lag1 = lag(meta_inflacion, 1),

        # Desviación de la expectativa a 24 meses respecto de la meta
        brecha_exp24 =
            exp_inf_24m - meta_inflacion,

        # Persistencia de la expectativa respecto de la meta
        brecha_exp24_lag1 =
            exp24_lag1 - meta_inflacion_lag1,

        # Desviación de la inflación rezagada respecto de la meta
        brecha_infl_lag1 =
            infl_gt_lag1 - meta_inflacion_lag1,

        # Variación interanual rezagada del tipo de cambio
        tcdep_var_lag1 =
            lag(tcdep_var, 1)
    ) |>
    filter(
        complete.cases(
            brecha_exp24,
            brecha_exp24_lag1,
            brecha_infl_lag1,
            tcdep_var_lag1
        )
    )

cat(
    "\nPrimera fecha efectiva de la muestra: ",
    as.character(min(df_model$fecha)),
    "\n"
)

cat(
    "Última fecha efectiva de la muestra: ",
    as.character(max(df_model$fecha)),
    "\n"
)

cat(
    "Observaciones efectivas: ",
    nrow(df_model),
    "\n"
)

# ------------------------------------------------------------
# 2) Definir ventana móvil
# ------------------------------------------------------------
window_size <- 36   # 36 meses = 3 años

# Prueba de robustez:
# cambiar únicamente la línea anterior por:
# window_size <- 30

if (nrow(df_model) < window_size) {
    stop(
        "La muestra efectiva tiene menos observaciones que la ventana móvil."
    )
}

# ------------------------------------------------------------
# 3) Función para estimar la regresión en cada ventana
# ------------------------------------------------------------
rolling_fun <- function(data_window) {

    modelo <- lm(
        brecha_exp24 ~
            brecha_exp24_lag1 +
            brecha_infl_lag1 +
            tcdep_var_lag1,
        data = data_window
    )

    # Errores estándar Newey-West HAC
    resultado_hac <- lmtest::coeftest(
        modelo,
        vcov. = sandwich::NeweyWest(
            modelo,
            lag = 3,
            prewhite = FALSE,
            adjust = TRUE
        )
    )

    resultados <- broom::tidy(resultado_hac)

    coef_beta <- resultados |>
        filter(term == "brecha_infl_lag1")

    coef_rho <- resultados |>
        filter(term == "brecha_exp24_lag1")

    coef_gamma <- resultados |>
        filter(term == "tcdep_var_lag1")

    ajuste <- broom::glance(modelo)

    tibble(
        # Sensibilidad frente a la brecha de inflación
        beta = coef_beta$estimate,
        se_beta = coef_beta$std.error,
        p_beta = coef_beta$p.value,

        # Persistencia de la brecha de expectativas
        rho = coef_rho$estimate,
        se_rho = coef_rho$std.error,
        p_rho = coef_rho$p.value,

        # Efecto de la depreciación cambiaria
        gamma = coef_gamma$estimate,
        se_gamma = coef_gamma$std.error,
        p_gamma = coef_gamma$p.value,

        # Ajuste del modelo
        r2_ajustado = ajuste$adj.r.squared
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
            fecha_fin = max(data_window$fecha),

            beta = res$beta,
            se_beta = res$se_beta,
            p_beta = res$p_beta,

            rho = res$rho,
            se_rho = res$se_rho,
            p_rho = res$p_rho,

            gamma = res$gamma,
            se_gamma = res$se_gamma,
            p_gamma = res$p_gamma,

            r2_ajustado = res$r2_ajustado
        )
    }
)

# ------------------------------------------------------------
# 5) Intervalos de confianza y medidas derivadas
# ------------------------------------------------------------
rolling_results <- rolling_results |>
    mutate(
        # Intervalos de confianza de beta
        beta_ic_inf = beta - 1.96 * se_beta,
        beta_ic_sup = beta + 1.96 * se_beta,

        # Intervalos de confianza de rho
        rho_ic_inf = rho - 1.96 * se_rho,
        rho_ic_sup = rho + 1.96 * se_rho,

        # Intervalos de confianza de gamma
        gamma_ic_inf = gamma - 1.96 * se_gamma,
        gamma_ic_sup = gamma + 1.96 * se_gamma,

        # Velocidad mensual de ajuste hacia la meta
        velocidad_ajuste = 1 - rho,

        # Vida media de una desviación, en meses
        vida_media = if_else(
            rho > 0 & rho < 1,
            log(0.5) / log(rho),
            NA_real_
        ),

        # Multiplicador de largo plazo
        multiplicador_lp = if_else(
            abs(rho) < 1,
            beta / (1 - rho),
            NA_real_
        ),

        # Indicador de cautela para el multiplicador:
        # TRUE cuando el intervalo de rho incluye la unidad
        rho_ic_incluye_uno =
            rho_ic_inf <= 1 & rho_ic_sup >= 1
    )

# ------------------------------------------------------------
# 6) Guardar resultados
# ------------------------------------------------------------
if (!dir.exists("output")) {
    dir.create("output")
}

write_csv(
    rolling_results,
    paste0(
        "output/rolling_brecha_infl_tc_",
        window_size,
        "m.csv"
    )
)

# ------------------------------------------------------------
# 7) Gráfico de beta
# ------------------------------------------------------------
p_beta <- ggplot(
    rolling_results,
    aes(x = fecha_fin, y = beta)
) +
    geom_ribbon(
        aes(
            ymin = beta_ic_inf,
            ymax = beta_ic_sup
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
            "Rolling regression (ventana de ",
            window_size,
            " meses)"
        ),
        subtitle = paste0(
            "Coeficiente de la brecha de inflación rezagada ",
            "sobre la brecha de expectativas a 24 meses"
        ),
        x = "Fecha final de la ventana",
        y = expression(beta[w])
    ) +
    theme_minimal()

ggsave(
    paste0(
        "output/rolling_beta_tc_",
        window_size,
        "m.png"
    ),
    plot = p_beta,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_beta)

# ------------------------------------------------------------
# 8) Gráfico de rho
# ------------------------------------------------------------
p_rho <- ggplot(
    rolling_results,
    aes(x = fecha_fin, y = rho)
) +
    geom_ribbon(
        aes(
            ymin = rho_ic_inf,
            ymax = rho_ic_sup
        ),
        alpha = 0.2
    ) +
    geom_line(linewidth = 1) +
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
            "Persistencia de las expectativas (ventana de ",
            window_size,
            " meses)"
        ),
        subtitle = paste0(
            "Coeficiente de la brecha de expectativas ",
            "a 24 meses rezagada"
        ),
        x = "Fecha final de la ventana",
        y = expression(rho[w])
    ) +
    theme_minimal()

ggsave(
    paste0(
        "output/rolling_rho_tc_",
        window_size,
        "m.png"
    ),
    plot = p_rho,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_rho)

# ------------------------------------------------------------
# 9) Resumen en consola
# ------------------------------------------------------------
cat(
    "\nVentana móvil utilizada: ",
    window_size,
    " meses\n"
)

cat(
    "Número de estimaciones móviles: ",
    nrow(rolling_results),
    "\n"
)

cat(
    "Primer coeficiente estimado: ",
    as.character(min(rolling_results$fecha_fin)),
    "\n"
)

cat(
    "Último coeficiente estimado: ",
    as.character(max(rolling_results$fecha_fin)),
    "\n"
)

cat(
    "\nArchivos generados en la carpeta output:\n",
    "- Resultados completos en CSV\n",
    "- Gráfico de beta\n",
    "- Gráfico de rho\n"
)
