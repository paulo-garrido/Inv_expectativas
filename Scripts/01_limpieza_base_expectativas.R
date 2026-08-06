# ============================================================
# TESINA PES (Guatemala) — Preparación de base para modelo FMI
# Tema: Expectativas de inflación (EEE) — determinantes y anclaje
# Autor: Paulo Garrido (plantilla)
# Fecha: (actualiza)
#
# INSUMO:
#   - Archivo: bd.xlsx
#   - Columnas esperadas (en este orden):
#     fecha, exp_12, exp_24, inf, Tasa_pol, imae, def_fisc, tc, ff_rates, inf_usa, wti_usd_bbl
#
# SALIDA:
#   - data/df_model.csv (base limpia para modelar)
#   - data/df_full_clean.csv (base limpia con transformaciones)
# ============================================================

# 0) Paquetes ------------------------------------------------
# Instala si hace falta:
# install.packages(c("readxl","dplyr","lubridate","janitor","tidyr","stringr"))

library(readxl)
library(dplyr)
library(lubridate)
library(janitor)
library(tidyr)
library(stringr)

# 1) Parámetros ----------------------------------------------
path_in  <- "bd.xlsx"
start_estimation <- as.Date("2011-01-01")   # muestra efectiva para estimar
end_estimation   <- NA                      # opcional: fija un fin, ej as.Date("2024-12-01")

dir.create("data", showWarnings = FALSE)

# 2) Cargar y estandarizar -----------------------------------
raw <- read_excel(path_in) %>%
  clean_names()

# Nota: en tu archivo la fecha viene como "unnamed_0" (columna 1)
df <- raw %>%
  rename(
    fecha        = unnamed_0,
    exp_12m      = exp_12,
    exp_24m      = exp_24,
    infl_gt      = inf,
    mpr_gt       = tasa_pol,
    imae_idx     = imae,
    def_fisc     = def_fisc,
    fx_ref       = tc,
    effr         = ff_rates,
    infl_usa     = inf_usa,
    wti_usd_bbl  = wti_usd_bbl
  ) %>%
  mutate(
    # Asegura que "fecha" sea Date. Si fuera numérica (Excel), descomenta:
    # fecha = as.Date(fecha, origin = "1899-12-30")
    fecha = as.Date(fecha)
  ) %>%
  arrange(fecha)

# 3) Chequeos rápidos de integridad ---------------------------
# 3.1 Duplicados de fecha
if (nrow(df) != n_distinct(df$fecha)) {
  stop("Hay fechas duplicadas en la base. Revisa bd.xlsx.")
}

# 3.2 Secuencia mensual (no es fatal si faltan meses, pero lo reporta)
seq_m <- seq.Date(min(df$fecha, na.rm = TRUE), max(df$fecha, na.rm = TRUE), by = "month")
missing_months <- setdiff(seq_m, df$fecha)
if (length(missing_months) > 0) {
  message("OJO: Faltan meses en la base: ", paste(format(missing_months, "%Y-%m"), collapse = ", "))
}

# 4) Transformaciones (alineadas al FMI) ----------------------
# Convenciones:
# - Tasas de interés: cambios en nivel (puntos porcentuales)
# - Tipo de cambio y petróleo: variación interanual en log (%) aprox
# - Inflación USA: el archivo trae la TASA (no el índice), por eso usamos su cambio mensual

df_full <- df %>%
  mutate(
    # --- Lags de expectativas
    exp_24m_l1 = lag(exp_24m, 1),
    exp_12m_l1 = lag(exp_12m, 1),

    # --- Inflación Guatemala (el FMI usa infl rezagada)
    infl_gt_l1 = lag(infl_gt, 1),

    # --- Política monetaria (Δ MPR)
    d_mpr = mpr_gt - lag(mpr_gt, 1),

    # --- Actividad: IMAE "latest print" (FMI usa t-2)
    # Recomendación: trabajar con crecimiento interanual del IMAE (más comparable en unidades)
    imae_yoy = 100 * (log(imae_idx) - lag(log(imae_idx), 12)),
    imae_yoy_l2 = lag(imae_yoy, 2),

    # Alternativa (si deseas replicar literalmente con el índice): nivel del IMAE con rezago 2
    imae_idx_l2 = lag(imae_idx, 2),

    # --- Tipo de cambio: depreciación interanual (% aprox)
    fxdep_yoy = 100 * (log(fx_ref) - lag(log(fx_ref), 12)),
    # Alternativa exacta (aritmética) si deseas:
    fxdep_yoy_simple = 100 * (fx_ref / lag(fx_ref, 12) - 1),

    # --- EFFR: cambio mensual en puntos porcentuales
    d_effr = effr - lag(effr, 1),

    # --- Inflación USA: el archivo ya trae "inflation rate" (p.ej. 2.6)
    # El FMI usa el cambio mensual en esa tasa:
    d_infl_usa = infl_usa - lag(infl_usa, 1),

    # --- Petróleo (WTI): variación interanual (% aprox)
    oil_yoy = 100 * (log(wti_usd_bbl) - lag(log(wti_usd_bbl), 12))
  )

# 5) Definir muestra efectiva ---------------------------------
df_est <- df_full %>%
  filter(fecha >= start_estimation) %>%
  { if (!is.na(end_estimation)) filter(., fecha <= end_estimation) else . }

# 6) Base final para modelar (baseline FMI) -------------------
# Ecuación FMI (versión replicable con tus variables):
# exp_24m_t ~ infl_gt_{t-1} + d_mpr_t + exp_24m_{t-1} + imae_{t-2} + fxdep_yoy_t +
#             def_fisc_t + d_effr_t + d_infl_usa_t + oil_yoy_t
#
# Nota: aquí uso imae_yoy_l2 como default; si quieres el índice, cambia por imae_idx_l2.

df_model <- df_est %>%
  transmute(
    fecha,
    exp_12m,
    exp_24m,

    # RHS (baseline FMI)
    infl_gt_l1,
    d_mpr,
    exp_24m_l1,
    imae_yoy_l2,
    fxdep_yoy,
    def_fisc,
    d_effr,
    d_infl_usa,
    oil_yoy,

    # (Opcional) dejar variables en nivel para inspección
    infl_gt,
    mpr_gt,
    imae_idx,
    fx_ref,
    effr,
    infl_usa,
    wti_usd_bbl
  ) %>%
  drop_na(exp_24m, infl_gt_l1, d_mpr, exp_24m_l1, imae_yoy_l2, fxdep_yoy, def_fisc, d_effr, d_infl_usa, oil_yoy)

# 7) Resumen de la muestra final ------------------------------
message("Muestra final para modelar: ",
        format(min(df_model$fecha), "%Y-%m"), " a ",
        format(max(df_model$fecha), "%Y-%m"),
        " | N = ", nrow(df_model))

# 8) Guardar salidas ------------------------------------------
write.csv(df_full,  "data/df_full_clean.csv", row.names = FALSE)
write.csv(df_model, "data/df_model.csv",      row.names = FALSE)

# 9) (Opcional) gráficos rápidos de sanity check --------------
# Descomenta si deseas:
# plot(df_model$fecha, df_model$exp_24m, type="l", main="Expectativa 24m", xlab="", ylab="%")
# plot(df_model$fecha, df_model$infl_gt, type="l", main="Inflación GT", xlab="", ylab="%")
# plot(df_model$fecha, df_model$fxdep_yoy, type="l", main="Depreciación yoy (log) TC", xlab="", ylab="%")
