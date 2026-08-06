# ============================================================
#PREPARACIÓN PARA PIB TRIMESTRAL
# ============================================================

# ------------------------------------------------------------
# 0. Limpieza
# ------------------------------------------------------------
rm(list = ls())
gc()

# ------------------------------------------------------------
# 1. Paquetes
# ------------------------------------------------------------
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(zoo)

# ------------------------------------------------------------
# 2. Archivo
# ------------------------------------------------------------

archivo <- "Empalme_CNT_1T_2001_3T_2025.xlsx"

# ------------------------------------------------------------
# 3. Leer hoja 3
#    skip = 5 porque los encabezados útiles empiezan en la fila 6
# ------------------------------------------------------------
pib_raw <- read_excel(
    path  = archivo,
    sheet = 3,
    skip  = 5
)

# Revisar nombres de columnas
names(pib_raw)

# ------------------------------------------------------------
# 4. Seleccionar columnas relevantes
#    - Período
#    - PIB Trimestral
# ------------------------------------------------------------
pib_df <- pib_raw %>%
    select(
        periodo = `Período`,
        pib = `PIB \r\nTrimestral`
    )

# ------------------------------------------------------------
# 5. Limpiar y reconstruir año-trimestre
# ------------------------------------------------------------
pib_largo <- pib_df %>%
    mutate(
        periodo = as.character(periodo),
        
        # identificar filas que contienen un año
        anio = if_else(
            str_detect(periodo, "^[0-9]{4}$"),
            as.integer(periodo),
            NA_integer_
        )
    ) %>%
    fill(anio, .direction = "down") %>%   # arrastrar el año a sus trimestres
    filter(periodo %in% c("I", "II", "III", "IV")) %>%   # quedarnos solo con trimestres
    mutate(
        trimestre = case_when(
            periodo == "I"   ~ 1L,
            periodo == "II"  ~ 2L,
            periodo == "III" ~ 3L,
            periodo == "IV"  ~ 4L
        )
    ) %>%
    arrange(anio, trimestre)

# ------------------------------------------------------------
# 6. Crear variable fecha trimestral
# ------------------------------------------------------------
pib_largo <- pib_largo %>%
    mutate(
        fecha = as.yearqtr(paste(anio, trimestre), format = "%Y %q"),
        fecha = as.Date(fecha)
    )

# ------------------------------------------------------------
# 7. Calcular tasas de variación
#    a) trimestral: respecto al trimestre anterior
#    b) interanual: respecto al mismo trimestre del año previo
# ------------------------------------------------------------
pib_largo <- pib_largo %>%
    mutate(
        var_trimestral = (pib / lag(pib) - 1) * 100,
        var_interanual = (pib / lag(pib, 4) - 1) * 100
    )

# ------------------------------------------------------------
# 8. Orden final de columnas
# ------------------------------------------------------------
pib_largo <- pib_largo %>%
    select(fecha, anio, trimestre, pib, var_trimestral, var_interanual)

# ------------------------------------------------------------
# 9. Ver resultado
# ------------------------------------------------------------
print(pib_largo)

# ------------------------------------------------------------
# 10. Guardar resultado
# ------------------------------------------------------------
write.csv(pib_largo, "pib_trimestral_2001-2025.csv", row.names = FALSE)
