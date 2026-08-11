# 01. ROLLING REGRESSIONS
# Sensibilidad de las expectativas a 24 meses frente a la inflación
# ============================================================
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

    ) |>
    filter(
        complete.cases(
            brecha_exp24,
            brecha_exp24_lag1,
            brecha_infl_lag1
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
window_size <- 60  # 60 meses = 5 años

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
            brecha_infl_lag1,
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
# Estilo gráfico común del working paper
# ------------------------------------------------------------

library(ggplot2)
library(scales)

# Paleta coherente con las figuras anteriores
col_serie  <- "#2F4A67"
col_banda  <- "#AFC4D6"
col_ref    <- "#4D4D4D"

tema_wp <- theme_minimal(base_size = 10.5) +
    theme(
        # Fondo y cuadrícula
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(
            colour = "grey92",
            linewidth = 0.30
        ),
        panel.grid.major.y = element_line(
            colour = "grey88",
            linewidth = 0.35
        ),
        
        # Ejes
        axis.title.x = element_blank(),
        axis.title.y = element_text(
            size = 10.5,
            colour = "black",
            margin = margin(r = 8)
        ),
        axis.text = element_text(
            size = 9,
            colour = "grey25"
        ),
        axis.ticks = element_blank(),
        
        # Sin títulos internos:
        # el título y la nota se administran desde Quarto
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        
        # Márgenes compactos
        plot.margin = margin(
            t = 4,
            r = 8,
            b = 4,
            l = 4
        )
    )


# ------------------------------------------------------------
# Figura: sensibilidad de las expectativas
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
        fill = "#AFC4D6",
        alpha = 0.28,
        linewidth = 0
    ) +
    
    geom_line(
        colour = "#2F4A67",
        linewidth = 0.80,
        lineend = "round"
    ) +
    
    geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "grey30",
        linewidth = 0.45
    ) +
    
    scale_x_date(
        date_breaks = "2 years",
        date_labels = "%Y",
        date_minor_breaks = "1 year",
        expand = expansion(mult = c(0.01, 0.015))
    ) +
    
    scale_y_continuous(
        breaks = seq(0, 0.4, by = 0.1),
        expand = expansion(mult = c(0.03, 0.05))
    ) +
    
    coord_cartesian(
        ylim = c(-0.015, NA),
        clip = "off"
    ) +
    
    labs(
        x = NULL,
        y = expression(hat(beta)[t])
    ) +
    
    theme_minimal(base_size = 10.5) +
    
    theme(
        panel.grid.minor.x = element_line(
            colour = "grey94",
            linewidth = 0.25
        ),
        panel.grid.minor.y = element_blank(),
        
        panel.grid.major.x = element_line(
            colour = "grey90",
            linewidth = 0.30
        ),
        panel.grid.major.y = element_line(
            colour = "grey87",
            linewidth = 0.35
        ),
        
        axis.title.y = element_text(
            size = 10.5,
            colour = "black",
            margin = margin(r = 8)
        ),
        
        axis.text = element_text(
            size = 9,
            colour = "grey25"
        ),
        
        axis.ticks = element_blank(),
        
        plot.margin = margin(
            t = 4,
            r = 8,
            b = 4,
            l = 4
        )
    )


ggsave(
    paste0(
        "output/rolling_beta_",
        window_size,
        "m.png"
    ),
    plot = p_beta,
    width = 8.5,
    height = 4.3,
    dpi = 300,
    bg = "white"
)

print(p_beta)


# ------------------------------------------------------------
# Figura: persistencia de las expectativas
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
        fill = "#AFC4D6",
        alpha = 0.28,
        linewidth = 0
    ) +
    
    geom_line(
        colour = "#2F4A67",
        linewidth = 0.80,
        lineend = "round"
    ) +
    
    geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "grey30",
        linewidth = 0.45
    ) +
    
    geom_hline(
        yintercept = 1,
        linetype = "dotted",
        colour = "grey50",
        linewidth = 0.45
    ) +
    
    scale_x_date(
        date_breaks = "2 years",
        date_labels = "%Y",
        date_minor_breaks = "1 year",
        expand = expansion(mult = c(0.01, 0.015))
    ) +
    
    scale_y_continuous(
        breaks = seq(0, 1, by = 0.2),
        expand = expansion(mult = c(0.03, 0.05))
    ) +
    
    coord_cartesian(
        ylim = c(-0.02, 1.05),
        clip = "off"
    ) +
    
    labs(
        x = NULL,
        y = expression(hat(rho)[t])
    ) +
    
    theme_minimal(base_size = 10.5) +
    
    theme(
        panel.grid.minor.x = element_line(
            colour = "grey94",
            linewidth = 0.25
        ),
        panel.grid.minor.y = element_blank(),
        
        panel.grid.major.x = element_line(
            colour = "grey90",
            linewidth = 0.30
        ),
        panel.grid.major.y = element_line(
            colour = "grey87",
            linewidth = 0.35
        ),
        
        axis.title.y = element_text(
            size = 10.5,
            colour = "black",
            margin = margin(r = 8)
        ),
        
        axis.text = element_text(
            size = 9,
            colour = "grey25"
        ),
        
        axis.ticks = element_blank(),
        
        plot.margin = margin(
            t = 4,
            r = 8,
            b = 4,
            l = 4
        )
    )

ggsave(
    paste0(
        "output/rolling_rho_",
        window_size,
        "m.png"
    ),
    plot = p_rho,
    width = 8.5,
    height = 4.3,
    dpi = 300,
    bg = "white"
)

ggsave(
    paste0(
        "output/rolling_rho_",
        window_size,
        "m.png"
    ),
    plot = p_rho,
    width = 8.5,
    height = 4.3,
    dpi = 300,
    bg = "white"
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

