# ============================================================
# TESINA PES (Guatemala)
# Tema: Expectativas de inflación (EEE) — determinantes y anclaje
# Autor: Paulo Garrido
# Objetivo del script:
#   Preparar una base mensual casi completa desde enero de 2010,
#   utilizando la información de 2009 únicamente para construir
#   la variación interanual del tipo de cambio.
# ============================================================

# ------------------------------------------------------------
# 0) Limpieza de entorno
# ------------------------------------------------------------
rm(list = ls())
cat("\014")
if (!is.null(dev.list())) dev.off()

# ------------------------------------------------------------
# 1) Paquetes y opciones
# ------------------------------------------------------------
library(readxl)
library(dplyr)
library(lubridate)
library(janitor)
library(stringr)
library(zoo)
library(readr)
library(ggplot2)

# ------------------------------------------------------------
# 2) Leer Excel y renombrar variables por orden de columnas
#
# Orden esperado:
#  1 fecha
#  2 exp_12
#  3 exp_24
#  4 inf
#  5 Tasa_pol
#  6 imae
#  7 def_fisc
#  8 tc
#  9 ff_rates
# 10 inf_usa
# 11 dco
# ------------------------------------------------------------

path_file <- "bd.xlsx"

df <- read_excel(path_file, sheet = 1) |>
    clean_names() |>
    dplyr::select(1:11) |>
    rlang::set_names(c(
        "fecha",
        "exp_inf_12m",
        "exp_inf_24m",
        "infl_gt",
        "tasa_pol",
        "imae",
        "deficit_fiscal",
        "tc",
        "effr",
        "infl_us",
        "oil"
    )) |>
    mutate(
        fecha = floor_date(as.Date(fecha), unit = "month")
    ) |>
    arrange(fecha)

# ------------------------------------------------------------
# 3) Validación básica: duplicados y huecos mensuales
# ------------------------------------------------------------

# Fechas duplicadas
dups <- df |>
    count(fecha) |>
    filter(n > 1)

if (nrow(dups) > 0) {
    warning(
        "Hay fechas duplicadas en la base. Revisa: ",
        paste(dups$fecha, collapse = ", ")
    )
}

# Huecos mensuales
full_index <- tibble(
    fecha = seq(
        min(df$fecha, na.rm = TRUE),
        max(df$fecha, na.rm = TRUE),
        by = "month"
    )
)

missing_months <- anti_join(
    full_index,
    df |> distinct(fecha),
    by = "fecha"
)

if (nrow(missing_months) > 0) {
    warning(
        "Faltan meses en la base. Ejemplos: ",
        paste(head(missing_months$fecha, 6), collapse = ", "),
        ifelse(nrow(missing_months) > 6, " ...", "")
    )
}

# ------------------------------------------------------------
# 4) Construcción de variables
# ------------------------------------------------------------
# Supuestos de unidades:
# - exp_inf_* e infl_* están en porcentaje (%).
# - tasa_pol y effr están en porcentaje anual.
# - tc y oil son niveles.
# - imae es un índice.
#
# Transformaciones:
# - d_tasa_pol: cambio mensual de la tasa de política (p.p.).
# - d_effr: cambio mensual de la tasa efectiva federal (p.p.).
# - tcdep_var: variación interanual del tipo de cambio,
#   calculada como diferencia logarítmica por 100.
# - oil_var: variación interanual del precio del petróleo,
#   calculada como diferencia logarítmica por 100.
# - d_infl_us: cambio mensual de la inflación de EE. UU. (p.p.).
# - exp24_lag1: expectativa a 24 meses rezagada un mes.
# - infl_gt_lag1: inflación de Guatemala rezagada un mes.
# - imae_lag2: IMAE rezagado dos meses.
#
# Importante:
# La información del tipo de cambio desde enero de 2009 permite
# calcular tcdep_var desde enero de 2010.
# ------------------------------------------------------------

df_2 <- df |>
    mutate(
        # Diferencias mensuales en tasas de interés
        d_tasa_pol = tasa_pol - lag(tasa_pol, 1),
        d_effr     = effr - lag(effr, 1),

        # Variaciones interanuales: log-diferencia por 100
        tcdep_var = 100 * (
            log(tc) - lag(log(tc), 12)
        ),

        oil_var = 100 * (
            log(oil) - lag(log(oil), 12)
        ),

        # Cambio mensual de la inflación de Estados Unidos
        d_infl_us = infl_us - lag(infl_us, 1),

        # Rezagos
        exp24_lag1   = lag(exp_inf_24m, 1),
        infl_gt_lag1 = lag(infl_gt, 1),
        imae_lag2    = lag(imae, 2)
    )

# ------------------------------------------------------------
# 5) Definir muestra general
# ------------------------------------------------------------
# Se conserva 2009 únicamente durante la construcción de
# tcdep_var. El archivo final comienza en enero de 2010.
#
# No se eliminan observaciones mediante un filtro global de
# casos completos. Cada modelo econométrico debe seleccionar
# únicamente las variables que utiliza y definir su propia
# muestra efectiva.
# ------------------------------------------------------------

df_model <- df_2 |>
    filter(fecha >= ymd("2010-01-01")) |>
    select(
        fecha,
        exp_inf_12m,
        exp_inf_24m,
        infl_gt,
        tasa_pol,
        imae,
        deficit_fiscal,
        tc,
        effr,
        infl_us,
        oil,

        # Variables construidas
        d_tasa_pol,
        imae_lag2,
        tcdep_var,
        d_effr,
        d_infl_us,
        oil_var,
        exp24_lag1,
        infl_gt_lag1
    )

# ------------------------------------------------------------
# 6) Comprobación de la muestra efectiva para rolling regressions
# ------------------------------------------------------------
# Modelo propuesto:
# brecha_exp24 ~ brecha_exp24_lag1 +
#                brecha_infl_lag1 +
#                tcdep_var_lag1
#
# Esta comprobación no modifica df_model; únicamente informa
# la primera y última fecha potencialmente utilizables.
# ------------------------------------------------------------

df_check_rolling <- df_model |>
    mutate(
        tcdep_var_lag1 = lag(tcdep_var, 1)
    ) |>
    filter(
        if_all(
            all_of(c(
                "exp_inf_24m",
                "exp24_lag1",
                "infl_gt_lag1",
                "tcdep_var_lag1"
            )),
            ~ !is.na(.x)
        )
    )

cat(
    "\nMuestra general guardada: ",
    format(min(df_model$fecha), "%Y-%m"),
    " a ",
    format(max(df_model$fecha), "%Y-%m"),
    "\n"
)

cat(
    "Observaciones en bd_sample_modelo: ",
    nrow(df_model),
    "\n"
)

cat(
    "Primera fecha potencial para rolling regressions: ",
    format(min(df_check_rolling$fecha), "%Y-%m"),
    "\n"
)

cat(
    "Última fecha potencial para rolling regressions: ",
    format(max(df_check_rolling$fecha), "%Y-%m"),
    "\n"
)

cat(
    "Observaciones potenciales para rolling regressions: ",
    nrow(df_check_rolling),
    "\n"
)

# ------------------------------------------------------------
# 7) Gráfica de control
# ------------------------------------------------------------

p_tc <- ggplot(df_model, aes(x = fecha, y = tc)) +
    geom_line(linewidth = 1) +
    labs(
        title = "Tipo de cambio promedio mensual",
        x = "Fecha",
        y = "Quetzales por dólar estadounidense"
    ) +
    scale_x_date(
        limits = c(
            min(df_model$fecha),
            max(df_model$fecha)
        ),
        date_breaks = "1 year",
        date_labels = "%Y",
        expand = c(0, 0)
    ) +
    theme_minimal(base_size = 13) +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        )
    )

print(p_tc)

# ------------------------------------------------------------
# 8) Guardar output
# ------------------------------------------------------------

write_csv(
    df_model,
    "bd_sample_modelo.csv",
    na = ""
)

cat(
    "\nArchivo generado correctamente: bd_sample_modelo.csv\n"
)
