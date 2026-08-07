# =========================================================
# HECHOS ESTILIZADOS PARA EL WORKING PAPER
# Tema: Evolución de la inflación y las expectativas de inflación
# Autor: Paulo Garrido
# =========================================================

# ---------------------------
# 0. Paquetes
# ---------------------------
paquetes <- c(
  "readxl", "dplyr", "ggplot2", "tidyr", "lubridate",
  "scales", "zoo", "patchwork", "stringr"
)

instalar_si_falta <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}

invisible(lapply(paquetes, instalar_si_falta))

getwd()

# ---------------------------
# 1. Rutas de trabajo
# ---------------------------
ruta_base    <- "bd.xlsx"     # ajusta si fuera necesario
ruta_salida  <- "output"
ruta_figuras <- file.path(ruta_salida, "figuras")
ruta_tablas  <- file.path(ruta_salida, "tablas")

dir.create(ruta_salida, showWarnings = FALSE)
dir.create(ruta_figuras, showWarnings = FALSE)
dir.create(ruta_tablas, showWarnings = FALSE)

# ---------------------------
# 2. Cargar base
# ---------------------------
# Se asume que la hoja relevante es la primera
bd_raw <- read_excel(ruta_base, sheet = 1)

# Verificar nombres originales
names(bd_raw)

# ---------------------------
# 3. Limpieza básica
# ---------------------------
# Se asume que:
# - la primera columna contiene la fecha en formato serial de Excel
# - el resto contiene variables macro y de expectativas

bd <- bd_raw %>%
  rename(fecha_excel = 1)

# Si la fecha viene como número serial de Excel
if (is.numeric(bd$fecha_excel)) {
  bd <- bd %>%
    mutate(fecha = as.Date(fecha_excel, origin = "1899-12-30"))
} else {
  # Si viene como texto
  bd <- bd %>%
    mutate(fecha = suppressWarnings(as.Date(fecha_excel)))
}

bd <- bd %>%
  select(-fecha_excel) %>%
  relocate(fecha) %>%
  mutate(across(-fecha, ~ suppressWarnings(as.numeric(.x)))) %>%
  filter(!is.na(fecha)) %>%
  arrange(fecha)

# Convertir todo lo que se pueda a numérico
bd <- bd %>%
  mutate(across(-fecha, ~ suppressWarnings(as.numeric(.x))))

# Eliminar filas completamente vacías o con fecha faltante
bd <- bd %>%
  filter(!is.na(fecha)) %>%
  filter(if_any(-fecha, ~ !is.na(.x)))

# Ordenar por fecha
bd <- bd %>%
  arrange(fecha)

# Revisar estructura
glimpse(bd)

# ---------------------------
# 4. Definir meta de inflación
# ---------------------------
# Meta central histórica del Banco de Guatemala:
# 2010-2011: 5.0%
# 2012:      4.5%
# Desde 2013: 4.0%
# Si luego deseas una meta variable por año, esta sección se sustituye

bd <- bd %>%
  mutate(
    meta_central = case_when(
      year(fecha) %in% c(2010,2011) ~ 5.0,
      year(fecha) == 2012           ~ 4.5,
      year(fecha) >= 2013           ~ 4.0,
      TRUE                          ~NA_real_    
    ),
    banda_inf = meta_central -1,
    banda_sup = meta_central +1
  )
  

# ---------------------------
# 5. Construcción de variables estilizadas
# ---------------------------
# Se asume nomenclatura de la base:
# exp_12 = expectativas a 12 meses
# exp_24 = expectativas a 24 meses
# inf    = inflación observada
# Tasa_pol = tasa de política monetaria (si existe)

bd <- bd %>%
  mutate(
    gap_inf     = inf - meta,
    gap_exp_12  = exp_12 - meta,
    gap_exp_24  = exp_24 - meta,
    abs_gap_inf    = abs(gap_inf),
    abs_gap_exp_12 = abs(gap_exp_12),
    abs_gap_exp_24 = abs(gap_exp_24),
    spread_exp = exp_24 - exp_12
  )

# Volatilidad móvil de la inflación (12 meses)
bd <- bd %>%
  mutate(
    vol_inf_12m = zoo::rollapply(
      inf,
      width = 12,
      FUN = sd,
      align = "right",
      fill = NA,
      na.rm = TRUE
    ),
    mean_abs_gap_inf_12m = zoo::rollapply(
      abs_gap_inf,
      width = 12,
      FUN = mean,
      align = "right",
      fill = NA,
      na.rm = TRUE
    )
  )

# Indicador simple de "dentro de banda"
bd <- bd %>%
  mutate(
    inf_dentro_banda    = ifelse(!is.na(inf) & inf >= meta_inf & inf <= meta_sup, 1, 0),
    exp12_dentro_banda  = ifelse(!is.na(exp_12) & exp_12 >= meta_inf & exp_12 <= meta_sup, 1, 0),
    exp24_dentro_banda  = ifelse(!is.na(exp_24) & exp_24 >= meta_inf & exp_24 <= meta_sup, 1, 0)
  )

# ---------------------------
# 6. Submuestra útil para antecedentes
# ---------------------------
# Si quieres restringir al período EMEI pleno o a un período particular, ajusta aquí:
# ejemplo: desde 2011
# bd_plot <- bd %>% filter(fecha >= as.Date("2011-01-01"))

bd_plot <- bd

# ---------------------------
# 7. Tema gráfico
# ---------------------------
tema_wp <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 10.5, color = "gray30", hjust = 0),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "gray20"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.4),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# ---------------------------
# 8. GRÁFICAS
# ---------------------------

# 8.1 Inflación observada, meta y banda
g1 <- ggplot(bd_plot, aes(x = fecha)) +
  geom_ribbon(aes(ymin = meta_inf, ymax = meta_sup), alpha = 0.20) +
  geom_line(aes(y = inf, color = "Inflación observada"), linewidth = 1) +
  geom_line(aes(y = meta, color = "Meta central"), linetype = "dashed", linewidth = 0.9) +
  scale_color_manual(values = c(
    "Inflación observada" = "#2C3E50",
    "Meta central" = "#16A085"
  )) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  scale_y_continuous(labels = label_number(suffix = "%")) +
  labs(
    title = "Inflación observada, meta central y banda de tolerancia",
    x = NULL,
    y = "Porcentaje"
  ) +
  tema_wp +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(ruta_figuras, "fig_inflacion_meta.png"),
  plot = g1, width = 10, height = 5, dpi = 300
)

g1
# 8.2 Brecha de inflación respecto de la meta
g2 <- ggplot(bd_plot, aes(x = fecha, y = gap_inf)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col() +
  scale_y_continuous(labels = label_number(suffix = " p.p.")) +
  labs(
    title = "Brecha de inflación respecto de la meta",
    x = NULL, y = "Puntos porcentuales"
  ) +
  tema_wp

ggsave(
  filename = file.path(ruta_figuras, "fig_gap_inflacion.png"),
  plot = g2, width = 10, height = 4.5, dpi = 300
)

g2
# 8.3 Volatilidad móvil de la inflación
g3 <- ggplot(bd_plot, aes(x = fecha, y = vol_inf_12m)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  labs(
    title = "Volatilidad móvil de la inflación (desviación estándar de 12 meses)",
    x = NULL, y = "Desviación estándar"
  ) +
  tema_wp
g3
ggsave(
  filename = file.path(ruta_figuras, "fig_volatilidad_inflacion.png"),
  plot = g3, width = 10, height = 4.5, dpi = 300
)
g3

# 8.4 Expectativas a 12 y 24 meses vs meta
g4 <- ggplot(bd_plot, aes(x = fecha)) +
  geom_ribbon(aes(ymin = meta_inf, ymax = meta_sup), alpha = 0.20) +
  geom_line(aes(y = exp_12, color = "Expectativa 12 meses"), linewidth = 1) +
  geom_line(aes(y = exp_24, color = "Expectativa 24 meses"), linewidth = 1) +
  geom_line(aes(y = meta, color = "Meta central"), linetype = "dashed", linewidth = 0.9) +
  scale_y_continuous(labels = label_number(suffix = "%")) +
  labs(
    title = "Expectativas de inflación y meta",
    x = NULL, y = "Porcentaje"
  ) +
  tema_wp
g4
ggsave(
  filename = file.path(ruta_figuras, "fig_expectativas_meta.png"),
  plot = g4, width = 10, height = 5, dpi = 300
)

# 8.5 Brecha de expectativas respecto a la meta
g5 <- bd_plot %>%
  select(fecha, gap_exp_12, gap_exp_24) %>%
  pivot_longer(
    cols = c(gap_exp_12, gap_exp_24),
    names_to = "serie",
    values_to = "brecha"
  ) %>%
  mutate(
    serie = recode(
      serie,
      gap_exp_12 = "Brecha expectativa 12 meses",
      gap_exp_24 = "Brecha expectativa 24 meses"
    )
  ) %>%
  ggplot(aes(x = fecha, y = brecha, color = serie)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = label_number(suffix = " p.p.")) +
  labs(
    title = "Brechas de expectativas respecto de la meta central",
    x = NULL, y = "Puntos porcentuales"
  ) +
  tema_wp

g5

ggsave(
  filename = file.path(ruta_figuras, "fig_gap_expectativas.png"),
  plot = g5, width = 10, height = 4.8, dpi = 300
)

# 8.6 Distancia absoluta de expectativas a la meta
g6 <- bd_plot %>%
  select(fecha, abs_gap_exp_12, abs_gap_exp_24) %>%
  pivot_longer(
    cols = c(abs_gap_exp_12, abs_gap_exp_24),
    names_to = "serie",
    values_to = "dist_abs"
  ) %>%
  mutate(
    serie = recode(
      serie,
      abs_gap_exp_12 = "Distancia absoluta 12 meses",
      abs_gap_exp_24 = "Distancia absoluta 24 meses"
    )
  ) %>%
  ggplot(aes(x = fecha, y = dist_abs, color = serie)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = label_number(suffix = " p.p.")) +
  labs(
    title = "Distancia absoluta de las expectativas respecto de la meta",
    x = NULL, y = "Puntos porcentuales"
  ) +
  tema_wp

g6
ggsave(
  filename = file.path(ruta_figuras, "fig_distancia_expectativas_meta.png"),
  plot = g6, width = 10, height = 4.8, dpi = 300
)

# 8.7 Pendiente de la estructura temporal de expectativas
# (exp_24 - exp_12)
g7 <- ggplot(bd_plot, aes(x = fecha, y = spread_exp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = label_number(suffix = " p.p.")) +
  labs(
    title = "Pendiente de la estructura temporal de expectativas",
    subtitle = "Diferencia entre expectativas a 24 meses y 12 meses",
    x = NULL, y = "Puntos porcentuales"
  ) +
  tema_wp

g7
ggsave(
  filename = file.path(ruta_figuras, "fig_spread_expectativas.png"),
  plot = g7, width = 10, height = 4.5, dpi = 300
)

# 8.8 Inflación observada y expectativas en una sola figura
g8 <- bd_plot %>%
  select(fecha, inf, exp_12, exp_24) %>%
  pivot_longer(
    cols = c(inf, exp_12, exp_24),
    names_to = "serie",
    values_to = "valor"
  ) %>%
  mutate(
    serie = recode(
      serie,
      inf = "Inflación observada",
      exp_12 = "Expectativa 12 meses",
      exp_24 = "Expectativa 24 meses"
    )
  ) %>%
  ggplot(aes(x = fecha, y = valor, color = serie)) +
  geom_line(linewidth = 1) +
  geom_line(data = bd_plot, aes(x = fecha, y = meta),
            inherit.aes = FALSE, linetype = "dashed", linewidth = 0.9) +
  scale_y_continuous(labels = label_number(suffix = "%")) +
  labs(
    title = "Inflación observada y expectativas de inflación",
    x = NULL, y = "Porcentaje"
  ) +
  tema_wp

ggsave(
  filename = file.path(ruta_figuras, "fig_inflacion_y_expectativas.png"),
  plot = g8, width = 10, height = 5, dpi = 300
)

# 8.9 Gráfico opcional: tasa de política e inflación
if ("Tasa_pol" %in% names(bd_plot)) {
  g9 <- ggplot(bd_plot, aes(x = fecha)) +
    geom_line(aes(y = inf, color = "Inflación"), linewidth = 1) +
    geom_line(aes(y = Tasa_pol, color = "Tasa de política"), linewidth = 1) +
    geom_line(aes(y = meta, color = "Meta"), linetype = "dashed", linewidth = 0.9) +
    scale_y_continuous(labels = label_number(suffix = "%")) +
    labs(
      title = "Inflación, meta y tasa de política monetaria",
      x = NULL, y = "Porcentaje"
    ) +
    tema_wp
  
  ggsave(
    filename = file.path(ruta_figuras, "fig_inflacion_tasa_politica.png"),
    plot = g9, width = 10, height = 5, dpi = 300
  )
}

# ---------------------------
# 9. CUADROS DESCRIPTIVOS
# ---------------------------

# 9.1 Estadísticos descriptivos generales
tabla_desc_general <- bd_plot %>%
  summarise(
    inflacion_promedio = mean(inf, na.rm = TRUE),
    inflacion_sd = sd(inf, na.rm = TRUE),
    exp12_promedio = mean(exp_12, na.rm = TRUE),
    exp12_sd = sd(exp_12, na.rm = TRUE),
    exp24_promedio = mean(exp_24, na.rm = TRUE),
    exp24_sd = sd(exp_24, na.rm = TRUE),
    gap_inf_abs_prom = mean(abs_gap_inf, na.rm = TRUE),
    gap_exp12_abs_prom = mean(abs_gap_exp_12, na.rm = TRUE),
    gap_exp24_abs_prom = mean(abs_gap_exp_24, na.rm = TRUE),
    pct_inf_dentro_banda = 100 * mean(inf_dentro_banda, na.rm = TRUE),
    pct_exp12_dentro_banda = 100 * mean(exp12_dentro_banda, na.rm = TRUE),
    pct_exp24_dentro_banda = 100 * mean(exp24_dentro_banda, na.rm = TRUE)
  )

write.csv(
  tabla_desc_general,
  file = file.path(ruta_tablas, "tabla_descriptivos_generales.csv"),
  row.names = FALSE
)

# 9.2 Estadísticos por año
tabla_anual <- bd_plot %>%
  mutate(anio = year(fecha)) %>%
  group_by(anio) %>%
  summarise(
    inflacion_promedio = mean(inf, na.rm = TRUE),
    exp12_promedio = mean(exp_12, na.rm = TRUE),
    exp24_promedio = mean(exp_24, na.rm = TRUE),
    gap_inf_abs_prom = mean(abs_gap_inf, na.rm = TRUE),
    gap_exp12_abs_prom = mean(abs_gap_exp_12, na.rm = TRUE),
    gap_exp24_abs_prom = mean(abs_gap_exp_24, na.rm = TRUE),
    pct_inf_dentro_banda = 100 * mean(inf_dentro_banda, na.rm = TRUE),
    pct_exp12_dentro_banda = 100 * mean(exp12_dentro_banda, na.rm = TRUE),
    pct_exp24_dentro_banda = 100 * mean(exp24_dentro_banda, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  tabla_anual,
  file = file.path(ruta_tablas, "tabla_anual_hechos_estilizados.csv"),
  row.names = FALSE
)

# 9.3 Correlaciones simples
tabla_cor <- tibble(
  correlacion_inf_exp12 = cor(bd_plot$inf, bd_plot$exp_12, use = "complete.obs"),
  correlacion_inf_exp24 = cor(bd_plot$inf, bd_plot$exp_24, use = "complete.obs"),
  correlacion_exp12_exp24 = cor(bd_plot$exp_12, bd_plot$exp_24, use = "complete.obs")
)

write.csv(
  tabla_cor,
  file = file.path(ruta_tablas, "tabla_correlaciones.csv"),
  row.names = FALSE
)

# ---------------------------
# 10. Guardar base procesada
# ---------------------------
saveRDS(
  bd_plot,
  file = file.path(ruta_salida, "bd_hechos_estilizados.rds")
)

# ---------------------------
# 11. Mensaje final
# ---------------------------
cat("\nProceso completado.\n")
cat("Base procesada guardada en: ", file.path(ruta_salida, "bd_hechos_estilizados.rds"), "\n")
cat("Figuras guardadas en: ", ruta_figuras, "\n")
cat("Tablas guardadas en: ", ruta_tablas, "\n")


# =========================================================
# GRÁFICA HISTÓRICA DE INFLACIÓN CON MARCA DEL EMEI
# =========================================================

# Base histórica de inflación
ruta_inflacion_hist <- "inflacion_bd.xlsx"

infl_hist <- readxl::read_excel(ruta_inflacion_hist, sheet = 1)

# Construir fecha mensual
infl_hist <- infl_hist %>%
  dplyr::mutate(
    fecha = as.Date(paste(anio, mes_num, "01", sep = "-"))
  ) %>%
  dplyr::arrange(fecha)

# Unir metas ya creadas en tu script principal
infl_hist_plot <- infl_hist %>%
  dplyr::left_join(
    bd %>% dplyr::select(fecha, meta, meta_inf, meta_sup),
    by = "fecha"
  )

# Fecha de inicio del EMEI
fecha_emei <- as.Date("2005-01-01")

# Gráfica
g_hist_inf <- ggplot(infl_hist_plot, aes(x = fecha, y = inflacion)) +
  
  # Banda meta sombreada desde el inicio del EMEI
  geom_ribbon(
    data = infl_hist_plot %>% dplyr::filter(fecha >= fecha_emei),
    aes(x = fecha, ymin = meta_inf, ymax = meta_sup),
    inherit.aes = FALSE,
    alpha = 0.18,
    fill = "#A9DFBF"
  ) +
  
  # Bordes superior e inferior de la banda
  geom_line(
    data = infl_hist_plot %>% dplyr::filter(fecha >= fecha_emei),
    aes(x = fecha, y = meta_inf),
    inherit.aes = FALSE,
    color = "#0047AB",
    linetype = "dotted",
    linewidth = 0.8
  ) +
  geom_line(
    data = infl_hist_plot %>% dplyr::filter(fecha >= fecha_emei),
    aes(x = fecha, y = meta_sup),
    inherit.aes = FALSE,
    color = "#0047AB",
    linetype = "dotted",
    linewidth = 0.8
  ) +
  
  # Meta central
  geom_line(
    data = infl_hist_plot %>% dplyr::filter(fecha >= fecha_emei),
    aes(x = fecha, y = meta),
    inherit.aes = FALSE,
    color = "#0047AB",
    linetype = "dashed",
    linewidth = 0.9
  ) +
  
  # Inflación observada
  geom_line(
    color = "#2C3E50",
    linewidth = 1
  ) +
  
  # Línea vertical de inicio del EMEI
  geom_vline(
    xintercept = fecha_emei,
    linetype = "dashed",
    linewidth = 0.9,
    color = "#B03A2E"
  ) +
  
  # Etiqueta del inicio del EMEI
  geom_label(
    data = data.frame(
      fecha = as.Date("2005-06-01"),
      y = max(infl_hist_plot$inflacion, na.rm = TRUE) * 0.98,
      texto = "Inicio del EMEI"
    ),
    aes(x = fecha, y = y, label = texto),
    inherit.aes = FALSE,
    size = 4,
    color = "#B03A2E",
    fill = "white",
    label.size = 0.2,
    fontface = "bold"
  ) +
  
  # Etiqueta del rango meta
  geom_label(
    data = data.frame(
      fecha = as.Date("2015-01-01"),
      y = 8.25,
      texto = "Rango meta de inflación"
    ),
    aes(x = fecha, y = y, label = texto),
    inherit.aes = FALSE,
    size = 3.8,
    color = "#0047AB",
    fill = "white",
    label.size = 0.2,
    fontface = "bold"
  ) +
  
  # Eje horizontal y vertical
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y"
  ) +
  scale_y_continuous(
    labels = scales::label_number(suffix = "%")
  ) +
  
  # Títulos
  labs(
    x = NULL,
    y = "Porcentaje"
  ) +
  
  # Tema
  tema_wp +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 

# Mostrar
g_hist_inf

# Guardar
ggsave(
  filename = file.path(ruta_figuras, "fig_inflacion_historica_1997_banda.png"),
  plot = g_hist_inf,
  width = 10,
  height = 5,
  dpi = 300
)


file.exists("output/figuras/fig_inflacion_historica_1997_banda.png")

list.files("output/figuras", full.names = TRUE)

list.files(recursive = TRUE, pattern = "fig_inflacion_historica_1997_banda.png")

#Expectativas de inflación a 12 y 24 meses frente a la meta, con banda

g3 <- ggplot(bd_plot, aes(x = fecha)) +
  # Banda meta
  geom_ribbon(
    aes(ymin = meta_inf, ymax = meta_sup),
    alpha = 0.22,
    fill = "#E8F1FC"
  ) +
  
  # Bordes de la banda
  geom_line(
    aes(y = meta_inf),
    color = "#4F81BD",
    linetype = "dotted",
    linewidth = 0.55
  ) +
  geom_line(
    aes(y = meta_sup),
    color = "#4F81BD",
    linetype = "dotted",
    linewidth = 0.55
  ) +
  
  # Meta central
  geom_line(
    aes(y = meta),
    color = "#0047AB",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  # Expectativas
  geom_line(
    aes(y = exp_12, color = "Expectativa a 12 meses"),
    linewidth = 1.05
  ) +
  geom_line(
    aes(y = exp_24, color = "Expectativa a 24 meses"),
    linewidth = 0.9
  ) +
  
  # Colores
  scale_color_manual(values = c(
    "Expectativa a 12 meses" = "#1F3A5F",
    "Expectativa a 24 meses" = "#74859A",
    "Meta central" = "#0047AB"
  )) +
  
  # Escalas
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = scales::label_number(suffix = "%"),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  
  labs(
    x = NULL,
    y = "Porcentaje"
  ) +
  geom_label(
    data = data.frame(
      fecha = as.Date("2018-01-01"),
      y = 3.86,
      texto = "Meta central"
    ),
    aes(x = fecha, y = y, label = texto),
    inherit.aes = FALSE,
    size = 2,
    color = "#0047AB",
    fill = "white",
    label.size = 0.2,
    fontface = "bold"
  ) +
  
  tema_wp +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.margin = margin(t = 6, b = 6),
    plot.margin = margin(10, 10, 25, 10)
    )
g3
ggsave(
  filename = file.path(ruta_figuras, "3.fig_expectat_meta.png"),
  plot = g3,
  width = 10,
  height = 5,
  dpi = 300
)



bd_gap_exp <- bd_plot %>%
  select(fecha, gap_exp_12, gap_exp_24) %>%
  pivot_longer(
    cols = c(gap_exp_12, gap_exp_24),
    names_to = "serie",
    values_to = "brecha"
  ) %>%
  mutate(
    serie = recode(
      serie,
      gap_exp_12 = "Brecha expectativa 12 meses",
      gap_exp_24 = "Brecha expectativa 24 meses"
    )
  )

bd_dist_abs <- bd_plot %>%
  dplyr::select(fecha, abs_gap_exp_12, abs_gap_exp_24) %>%
  tidyr::pivot_longer(
    cols = c(abs_gap_exp_12, abs_gap_exp_24),
    names_to = "serie",
    values_to = "dist_abs"
  ) %>%
  dplyr::mutate(
    serie = dplyr::recode(
      serie,
      abs_gap_exp_12 = "Distancia absoluta 12 meses",
      abs_gap_exp_24 = "Distancia absoluta 24 meses"
    )
  )

g5 <- ggplot(bd_dist_abs, aes(x = fecha, y = dist_abs, color = serie)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(
    "Distancia absoluta 12 meses" = "#1F3A5F",
    "Distancia absoluta 24 meses" = "#74859A"
  )) +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = scales::label_number(suffix = " p.p."),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Puntos porcentuales"
  ) +
  tema_wp +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.3, "cm"),
    plot.margin = margin(10, 10, 20, 10)
  )

g5

ggsave(
  filename = file.path(ruta_figuras, "4.distancia_abs_expectativas_meta.png"),
  plot = g5,
  width = 10,
  height = 5.6,
  dpi = 300
)


# =========================================================
# CUADRO 1. ESTADÍSTICOS DESCRIPTIVOS DE INFLACIÓN Y EXPECTATIVAS
# =========================================================

# ---------------------------
# A. Definir subperíodos
# ---------------------------
bd_cuadro <- bd_plot %>%
  mutate(
    subperiodo = case_when(
      fecha >= as.Date("2010-01-01") & fecha <= as.Date("2015-12-01") ~ "2010-2015",
      fecha >= as.Date("2016-01-01") & fecha <= as.Date("2020-12-01") ~ "2016-2020",
      fecha >= as.Date("2021-01-01") & fecha <= as.Date("2026-12-01") ~ "2021-2026",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(subperiodo))

# ---------------------------
# B. Función auxiliar
# ---------------------------
resumen_vars <- function(data, etiqueta_periodo) {
  tibble(
    Periodo = etiqueta_periodo,
    
    inflacion_media = mean(data$inf, na.rm = TRUE),
    inflacion_sd    = sd(data$inf, na.rm = TRUE),
    
    exp12_media = mean(data$exp_12, na.rm = TRUE),
    exp12_sd    = sd(data$exp_12, na.rm = TRUE),
    
    exp24_media = mean(data$exp_24, na.rm = TRUE),
    exp24_sd    = sd(data$exp_24, na.rm = TRUE),
    
    dist12_media = mean(data$abs_gap_exp_12, na.rm = TRUE),
    dist12_sd    = sd(data$abs_gap_exp_12, na.rm = TRUE),
    
    dist24_media = mean(data$abs_gap_exp_24, na.rm = TRUE),
    dist24_sd    = sd(data$abs_gap_exp_24, na.rm = TRUE),
    
    exp12_dentro_banda = 100 * mean(data$exp12_dentro_banda, na.rm = TRUE),
    exp24_dentro_banda = 100 * mean(data$exp24_dentro_banda, na.rm = TRUE)
  )
}

# ---------------------------
# C. Construir cuadro ancho
# ---------------------------
cuadro1_wide <- bind_rows(
  resumen_vars(
    bd_cuadro %>% filter(fecha >= as.Date("2010-01-01") & fecha <= as.Date("2026-12-01")),
    "2010-2026"
  ),
  resumen_vars(
    bd_cuadro %>% filter(subperiodo == "2010-2015"),
    "2010-2015"
  ),
  resumen_vars(
    bd_cuadro %>% filter(subperiodo == "2016-2020"),
    "2016-2020"
  ),
  resumen_vars(
    bd_cuadro %>% filter(subperiodo == "2021-2026"),
    "2021-2026"
  )
)

# ---------------------------
# D. Reorganizar para presentación
# ---------------------------
cuadro1_long <- cuadro1_wide %>%
  pivot_longer(
    cols = -Periodo,
    names_to = "estadistico",
    values_to = "valor"
  ) %>%
  mutate(
    bloque = case_when(
      estadistico %in% c("inflacion_media", "inflacion_sd") ~ "Inflación observada",
      estadistico %in% c("exp12_media", "exp12_sd", "exp12_dentro_banda") ~ "Expectativa a 12 meses",
      estadistico %in% c("exp24_media", "exp24_sd", "exp24_dentro_banda") ~ "Expectativa a 24 meses",
      estadistico %in% c("dist12_media", "dist12_sd") ~ "Distancia absoluta a la meta (12 meses)",
      estadistico %in% c("dist24_media", "dist24_sd") ~ "Distancia absoluta a la meta (24 meses)",
      TRUE ~ "Otros"
    ),
    medida = case_when(
      stringr::str_detect(estadistico, "_media$") ~ "Media",
      stringr::str_detect(estadistico, "_sd$") ~ "Desviación estándar",
      stringr::str_detect(estadistico, "dentro_banda") ~ "% dentro de la banda",
      TRUE ~ estadistico
    )
  ) %>%
  select(bloque, medida, Periodo, valor)

cuadro1_final <- cuadro1_long %>%
  pivot_wider(
    names_from = Periodo,
    values_from = valor
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 2))
  )

# ---------------------------
# E. Guardar cuadro
# ---------------------------
write.csv(
  cuadro1_final,
  file = file.path(ruta_tablas, "cuadro1_estadisticos_expectativas.csv"),
  row.names = FALSE
)

saveRDS(
  cuadro1_final,
  file = file.path(ruta_tablas, "cuadro1_estadisticos_expectativas.rds")
)

#| label: tbl-estadisticos-expectativas
#| tbl-cap: "Estadísticos descriptivos de inflación y expectativas de inflación por subperíodos."
knitr::kable(
  cuadro1_final,
  align = c("l", "l", "c", "c", "c", "c")
)
# Ver en consola
cuadro1_final

