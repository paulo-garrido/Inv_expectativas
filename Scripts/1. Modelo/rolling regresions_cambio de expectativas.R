# ============================================================
# ROLLING REGRESSIONS
# Sensibilidad de expectativas a 24 meses a la inflación pasada
# ============================================================

library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(broom)
library(purrr)

# ------------------------------------------------------------
# 0) Limpieza de entorno
# ------------------------------------------------------------
rm(list = ls())
cat("\014")
if (!is.null(dev.list())) dev.off()

# ------------------------------------------------------------
# 1) Cargar base
# ------------------------------------------------------------
df_model <- read_csv("bd_sample_modelo.csv", show_col_types = FALSE) |>
    mutate(fecha = ymd(fecha)) |>
    arrange(fecha)

# ------------------------------------------------------------
# 2) Definir ventana móvil
# ------------------------------------------------------------
window_size <- 60   # 60 meses = 5 años

# ------------------------------------------------------------
# 3) Función para estimar regresión en cada ventana
# ------------------------------------------------------------
rolling_fun <- function(data_window) {
    modelo <- lm(exp_inf_24m ~ infl_gt_lag1, data = data_window)
    
    broom::tidy(modelo) |>
        dplyr::filter(term == "infl_gt_lag1") |>
        dplyr::select(term, estimate, std.error, statistic, p.value)
}

rolling_results <- purrr::map_dfr(window_size:nrow(df_model), function(i) {
    
    data_window <- df_model[(i - window_size + 1):i, ]
    
    res <- rolling_fun(data_window)
    
    tibble::tibble(
        fecha_fin = max(data_window$fecha),
        estimate  = res$estimate,
        std.error = res$std.error,
        statistic = res$statistic,
        p.value   = res$p.value
    )
})


# ------------------------------------------------------------
# 4) Loop rolling
# ------------------------------------------------------------
rolling_results <- map_dfr(window_size:nrow(df_model), function(i) {
    
    data_window <- df_model[(i - window_size + 1):i, ]
    
    res <- rolling_fun(data_window)
    
    tibble(
        fecha_fin = max(data_window$fecha),
        estimate  = res$estimate,
        std.error = res$std.error,
        statistic = res$statistic,
        p.value   = res$p.value
    )
})

# ------------------------------------------------------------
# 5) Intervalos de confianza al 95%
# ------------------------------------------------------------
rolling_results <- rolling_results |>
    mutate(
        ic_inf = estimate - 1.96 * std.error,
        ic_sup = estimate + 1.96 * std.error
    )

# ------------------------------------------------------------
# 6) Guardar resultados
# ------------------------------------------------------------
if (!dir.exists("output")) dir.create("output")

write_csv(rolling_results, "output/rolling_infl_gt_lag1_60m.csv")

# ------------------------------------------------------------
# 7) Gráfico
# ------------------------------------------------------------
p_roll <- ggplot(rolling_results, aes(x = fecha_fin, y = estimate)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = ic_inf, ymax = ic_sup), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "Rolling regression (ventana de 60 meses)",
        subtitle = "Coeficiente de infl_gt_lag1 sobre exp_inf_24m",
        x = "Fecha final de la ventana",
        y = "Coeficiente estimado"
    ) +
    theme_minimal()

ggsave(
    "output/rolling_infl_gt_lag1_60m.png",
    plot = p_roll,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_roll)

