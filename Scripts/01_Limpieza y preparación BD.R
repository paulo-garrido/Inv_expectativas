# ============================================================
# TESINA PES (Guatemala) 
# Tema: Expectativas de inflación (EEE) — determinantes y anclaje
# Autor: Paulo Garrido (plantilla)
# Fecha: 9/03/26
# ============================================================

#-1) #limpieza entorno
rm(list = ls())
cat("\014")
if (!is.null(dev.list())) dev.off()

# ============================================================
# 0) Paquetes y opciones
# ============================================================
library(readxl)
library(dplyr)
library(lubridate)
library(janitor)
library(stringr)
library(zoo)
library(readr)

# ============================================================
# 1) Leer Excel y renombrar variables por ORDEN de columnas
#    Orden esperado:
#    1 fecha
#    2 exp_12m
#    3 exp_24m  ("24 meses")
#    4 infl_gt
#    5 tasa_pol (tasa de política)
#    6 imae
#    7 deficit_fiscal
#    8 tc (tipo de cambio)
#    9 effr (efective federal found rates)
#    10 infl_us
#    11 oil (precio petróleo)
# ============================================================

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
    ))

df <- df |>
    mutate(
        fecha = floor_date(as.Date(fecha), unit = "month")  # deja todo en YYYY-MM-01
    ) |>
    arrange(fecha)
df


# ============================================================
# 3) Validación básica: duplicados y huecos mensuales
# ============================================================

# Duplicados
dups <- df |> count(fecha) |> filter(n > 1)
if (nrow(dups) > 0) {
    warning("Hay fechas duplicadas en la base. Revisa: ", paste(dups$fecha, collapse = ", "))
}

# Huecos mensuales
full_index <- tibble(fecha = seq(min(df$fecha), max(df$fecha), by = "month"))
missing_months <- anti_join(full_index, df |> distinct(fecha), by = "fecha")
if (nrow(missing_months) > 0) {
    warning("Faltan meses en la base (huecos). Ejemplos: ",
            paste(head(missing_months$fecha, 6), collapse = ", "),
            ifelse(nrow(missing_months) > 6, " ...", ""))
}

# ============================================================
# 4) Construcción de variables (consistentes con FMI)
# ============================================================
# Supuestos de unidades:
# - exp_inf_* e infl_* están en porcentaje (%)
# - tasa_pol y effr están en porcentaje (% anual)
# - tc y oil son niveles (Q/USD y USD/barril, respectivamente)
# - imae es índice/nivel
# - deficit_fiscal: se transformó como % PIB
## Construcciones:
# - d_tasa_pol: cambio mensual (p.p.)
# - d_effr: cambio mensual (p.p.)
# - tcdep_var: variación interanual del tipo de cambio (en %, con log-dif)
# - oil_var: variación interanual del petróleo (en %, con log-dif)
# - d_infl_us: cambio mensual de inflación USA (en p.p.), si infl_us ya es tasa (var)
# - Rezagos: exp24_lag1, infl_gt_lag1, imae_lag2
# ============================================================

df_2 <- df |>
    mutate(
        # Diferencias en tasas de interés (puntos porcentuales)
        d_tasa_pol  = tasa_pol  - lag(tasa_pol, 1),
        d_effr = effr - lag(effr, 1),
        
        # Variaciones interanuales (log-dif * 100)
        tcdep_var  = 100 * (log(tc)  - lag(log(tc), 12)),
        oil_var    = 100 * (log(oil) - lag(log(oil), 12)),
        
        # Si infl_us ya es inflación, su cambio mensual:
        d_infl_us  = infl_us - lag(infl_us, 1),
        
        # Rezagos:
        exp24_lag1    = lag(exp_inf_24m, 1),
        infl_gt_lag1  = lag(infl_gt, 1),
        imae_lag2     = lag(imae, 2)
    )

# ============================================================
# 6) Definir muestra:
#    - se usará 2010 solo para construir var/diffs
#    - Modelo desde 2011
# ============================================================

df_model <- df_2 |>
    filter(fecha >= ymd("2011-01-01")) |>
    # Mantener solo columnas relevantes para el baseline FMI
    select(
        fecha,
        exp_inf_12m, exp_inf_24m,
        infl_gt, tasa_pol, imae, deficit_fiscal, tc, effr, infl_us, oil,
        # Variables construidas
        d_tasa_pol, imae_lag2, tcdep_var, d_effr, d_infl_us, oil_var,
        exp24_lag1, infl_gt_lag1
    ) 
     
key_vars <- c("exp_inf_24m", "infl_gt_lag1", "d_tasa_pol", "exp24_lag1",
              "imae_lag2", "tcdep_var", "deficit_fiscal",
              "d_effr", "d_infl_us", "oil_var")

df_model <- df_model |>
    filter(if_all(all_of(key_vars), ~ !is.na(.x)))

cat("\nObservaciones para modelar (sample efectivo): ", nrow(df_model), "\n")


#gràfica para analizar comportamiento de datos

ggplot(df_model, aes(x = fecha, y = tasa_pol)) +
    geom_line(color = "#0072B2", linewidth = 1.3) +
    labs(
        title = "Serie",
        x = "Fecha",
        y = "Valor"
    ) +
    scale_x_date(
        limits = c(min(df_model$fecha), max(df_model$fecha)),
        date_breaks = "1 year",
        date_labels = "%Y",
        expand = c(0, 0)
    ) +
    theme_minimal(base_size = 13) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

# ============================================================
# 8) Guardar outputs
# ============================================================

write_csv(df_model, "bd_sample_modelo.csv")


