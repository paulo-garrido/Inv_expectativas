# ============================================================
# CURVA DE PHILLIPS INICIAL - RÉPLICA TRIMESTRAL
# Guatemala | 2011-2025

# ============================================================

# ------------------------------------------------------------
# 0. Limpieza
# ------------------------------------------------------------
rm(list = ls())
gc()

# ------------------------------------------------------------
# 1. Paquetes
# ------------------------------------------------------------
packages <- c(
    "readxl", "dplyr", "tidyr", "stringr", "lubridate",
    "zoo", "readr", "broom", "ggplot2", "rvest", "xml2"
)

installed <- rownames(installed.packages())
for (p in packages) {
    if (!(p %in% installed)) install.packages(p)
}
invisible(lapply(packages, library, character.only = TRUE))

# ------------------------------------------------------------
# 2. Archivos
# ------------------------------------------------------------
archivo_mensual <- "bd.xlsx"
archivo_pib     <- "pib_trimestral_2001-2025.csv"
archivo_itcer   <- "itcer.xlsx"

# ------------------------------------------------------------
# 3. Leer base mensual principal
# ------------------------------------------------------------
bd_m <- read_excel(archivo_mensual, sheet = "BD") %>%
    mutate(
        fecha = as.Date(fecha),
        yq = as.yearqtr(fecha)
    ) %>%
    arrange(fecha)

# ------------------------------------------------------------
# 4. Trimestralizar base mensual
# Regla: último mes del trimestre
# ------------------------------------------------------------
bd_q <- bd_m %>%
    group_by(yq) %>%
    arrange(fecha, .by_group = TRUE) %>%
    summarise(
        fecha_trim = max(fecha),
        inf        = last(inf),       # inflación interanual
        exp_12     = last(exp_12),
        exp_24     = last(exp_24),
        tasa_pol   = last(Tasa_pol),
        .groups = "drop"
    ) %>%
    mutate(
        anio = year(fecha_trim),
        trimestre = quarter(fecha_trim)
    )

# ------------------------------------------------------------
# 5. Leer archivo ITCER
# Ojo: el .xls realmente viene como tabla HTML
# ------------------------------------------------------------
itcer_raw <- read_excel("itcer.xlsx")
names(itcer_raw)
# Aplanar nombres de columnas
colnames(itcer_raw) <- c(
    "fecha_txt",
    "itcer_usa",
    "itcer_global",
    "itcer_capard",
    "itcer_resto_mundo",
    "itcer_usa_yoy",
    "itcer_global_yoy",
    "itcer_capard_yoy",
    "itcer_resto_mundo_yoy",
    "itcer_usa_base",
    "itcer_global_base",
    "itcer_capard_base",
    "itcer_resto_mundo_base"
)

# ------------------------------------------------------------
# 6. Convertir fecha mensual del ITCER
# ------------------------------------------------------------
meses <- c(
    "Enero" = "01", "Febrero" = "02", "Marzo" = "03", "Abril" = "04",
    "Mayo" = "05", "Junio" = "06", "Julio" = "07", "Agosto" = "08",
    "Septiembre" = "09", "Octubre" = "10", "Noviembre" = "11", "Diciembre" = "12"
)

itcer_m <- itcer_raw %>%
    mutate(
        fecha_txt = as.character(fecha_txt),
        fecha_txt = str_squish(fecha_txt)
    ) %>%
    # quedarnos solo con filas tipo 2000-Enero
    filter(str_detect(fecha_txt, "^\\d{4}-[A-Za-zÁÉÍÓÚáéíóú]+$")) %>%
    mutate(
        anio = str_extract(fecha_txt, "^\\d{4}"),
        mes_nombre = str_replace(fecha_txt, "^\\d{4}-", ""),
        mes_num = unname(meses[mes_nombre]),
        fecha = as.Date(paste0(anio, "-", mes_num, "-01")),
        itcer_global = as.numeric(itcer_global)
    ) %>%
    select(fecha, itcer_global) %>%
    filter(!is.na(fecha), !is.na(itcer_global)) %>%
    arrange(fecha)

head(itcer_m)
tail(itcer_m)

# ------------------------------------------------------------
# 7. Pasar ITCER a trimestral
# Regla: último mes del trimestre
# ------------------------------------------------------------
itcer_q <- itcer_m %>%
    mutate(yq = as.yearqtr(fecha)) %>%
    group_by(yq) %>%
    arrange(fecha, .by_group = TRUE) %>%
    summarise(
        fecha = max(fecha),
        itcer_global = last(itcer_global),
        .groups = "drop"
    ) %>%
    mutate(
        anio = year(fecha),
        trimestre = quarter(fecha),
        
        # inverso del ITCER
        inv_itcer = 1 / itcer_global
    ) %>%
    arrange(fecha) %>%
    mutate(
        d_inv_itcer_yoy = (inv_itcer / lag(inv_itcer, 4) - 1) * 100
    ) %>%
    select(fecha, anio, trimestre, itcer_global, d_inv_itcer_yoy)

# ------------------------------------------------------------
# 8. Leer PIB trimestral
# ------------------------------------------------------------
pib_q <- read_csv(archivo_pib, show_col_types = FALSE) %>%
    mutate(fecha = as.Date(fecha)) %>%
    select(
        fecha,
        anio,
        trimestre,
        pib,
        g_pib_yoy = var_interanual
    )

# ------------------------------------------------------------
# 9. Unir todas las bases
# ------------------------------------------------------------
# ------------------------------------------------------------
# 9. Unir todas las bases
# ------------------------------------------------------------
df_q <- bd_q %>%
    rename(fecha = fecha_trim) %>%
    left_join(
        pib_q %>% select(anio, trimestre, pib, g_pib_yoy),
        by = c("anio", "trimestre")
    ) %>%
    left_join(
        itcer_q %>% select(anio, trimestre, itcer_global, d_inv_itcer_yoy),
        by = c("anio", "trimestre")
    ) %>%
    arrange(anio, trimestre)

# ------------------------------------------------------------
# 10. Valores de referencia
# ------------------------------------------------------------
meta_inflacion <- 4.0

g_potencial <- mean(df_q$g_pib_yoy, na.rm = TRUE)
dz_equilibrio <- mean(df_q$d_inv_itcer_yoy, na.rm = TRUE)

# ------------------------------------------------------------
# 11. Construcción de variables del modelo
# ------------------------------------------------------------
df_model <- df_q %>%
    mutate(
        pi_gap    = inf - meta_inflacion,
        y_gap     = g_pib_yoy - g_potencial,
        dz_gap    = d_inv_itcer_yoy - dz_equilibrio,
        pi_gap_l1 = lag(pi_gap, 1)
    ) %>%
    filter(fecha >= as.Date("2011-01-01"),
           fecha <= as.Date("2025-12-31")) %>%
    select(
        fecha, anio, trimestre,
        inf, exp_12, exp_24, tasa_pol,
        g_pib_yoy, itcer_global, d_inv_itcer_yoy,
        pi_gap, pi_gap_l1, y_gap, dz_gap
    ) %>%
    drop_na(pi_gap, pi_gap_l1, y_gap, dz_gap)

# ------------------------------------------------------------
# 12. Estimación de la Curva de Phillips inicial
# ------------------------------------------------------------
modelo_cp <- lm(pi_gap ~ pi_gap_l1 + y_gap + dz_gap, data = df_model)

summary(modelo_cp)
broom::tidy(modelo_cp)
broom::glance(modelo_cp)

cor_trim_imae_prom <- cor(bd_q$inf, bd_q$imae_prom, use = "complete.obs", method = "pearson")


# ------------------------------------------------------------
# 4. Trimestralizar base mensual
# Regla:
# - inflación: último mes del trimestre
# - IMAE: promedio del trimestre y fin de trimestre
# ------------------------------------------------------------
bd_q <- bd_m %>%
    group_by(yq) %>%
    arrange(fecha, .by_group = TRUE) %>%
    summarise(
        fecha_trim = max(fecha),
        inf        = last(inf),        # inflación interanual
        exp_12     = last(exp_12),
        exp_24     = last(exp_24),
        tasa_pol   = last(Tasa_pol),
        imae_prom  = mean(imae, na.rm = TRUE),
        imae_fin   = last(imae),
        .groups = "drop"
    ) %>%
    mutate(
        anio = year(fecha_trim),
        trimestre = quarter(fecha_trim)
    )

# Correlación trimestral inflación vs IMAE promedio
cor_trim_imae_prom <- cor(bd_q$inf, bd_q$imae_prom,
                          use = "complete.obs", method = "pearson")
cor_trim_imae_prom

# Correlación trimestral inflación vs IMAE fin de trimestre
cor_trim_imae_fin <- cor(bd_q$inf, bd_q$imae_fin,
                         use = "complete.obs", method = "pearson")
cor_trim_imae_fin

# Correlación trimestral inflación vs crecimiento del PIB
cor_trim_pib <- cor(df_q$inf, df_q$g_pib_yoy,
                    use = "complete.obs", method = "pearson")
cor_trim_pib
# ------------------------------------------------------------
# 13. Guardar salidas
# ------------------------------------------------------------
write_csv(df_model, "df_curva_phillips_itcer_2011_2025.csv")
write_csv(broom::tidy(modelo_cp), "coeficientes_cp_itcer.csv")

# ------------------------------------------------------------
# 14. Gráficas rápidas
# ------------------------------------------------------------
ggplot(df_model, aes(x = fecha, y = inf)) +
    geom_line(linewidth = 1) +
    labs(
        title = "Inflación interanual trimestral",
        subtitle = "Último mes de cada trimestre",
        x = "Fecha",
        y = "Inflación interanual (%)"
    ) +
    theme_minimal(base_size = 12)

ggplot(df_model, aes(x = fecha, y = d_inv_itcer_yoy)) +
    geom_line(linewidth = 1) +
    labs(
        title = "Variación interanual del inverso del ITCER global",
        x = "Fecha",
        y = "Porcentaje"
    ) +
    theme_minimal(base_size = 12)