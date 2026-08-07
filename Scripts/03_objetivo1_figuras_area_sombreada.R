# =========================================================
# OBJETIVO 1 - FIGURAS Y CUADRO DESCRIPTIVO PARA EL WORKING PAPER
# Evolución, proximidad y estabilidad de las expectativas de inflación
# a 12 y 24 meses respecto de la meta vigente
# =========================================================

# ------------------------------------------------------------
# 0) Limpieza de entorno
# ------------------------------------------------------------
rm(list = ls())
cat("\014")
if (!is.null(dev.list())) dev.off()

# ------------------------------------------------------------
# 1) Paquetes
# ------------------------------------------------------------
paquetes <- c(
  "dplyr", "ggplot2", "tidyr", "lubridate", "scales",
  "zoo", "patchwork", "stringr", "readxl", "tibble"
)

instalar_si_falta <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}

invisible(lapply(paquetes, instalar_si_falta))

# ------------------------------------------------------------
# 2) Rutas
# ------------------------------------------------------------
ruta_salida  <- "figuras y tablas"
ruta_figuras <- file.path(ruta_salida, "figuras")
ruta_tablas  <- file.path(ruta_salida, "tablas")
ruta_bd_rds  <- file.path(ruta_salida, "bd_hechos_estilizados.rds")

dir.create(ruta_salida,  showWarnings = FALSE)
dir.create(ruta_figuras, showWarnings = FALSE)
dir.create(ruta_tablas,  showWarnings = FALSE)

# ------------------------------------------------------------
# 3) Cargar base procesada
# ------------------------------------------------------------
if (!file.exists(ruta_bd_rds)) {
  stop(
    paste0(
      "No se encontró la base procesada en: ", ruta_bd_rds,
      "\nEjecuta primero el script generador de la base."
    )
  )
}

bd <- readRDS(ruta_bd_rds) %>%
  dplyr::arrange(fecha)

# ------------------------------------------------------------
# 4) Función de meta histórica
# ------------------------------------------------------------
meta_inflacion_historica <- function(fecha) {
  anio <- lubridate::year(fecha)

  dplyr::case_when(
    anio == 2005 ~ 5.0,
    anio == 2006 ~ 6.0,
    anio == 2007 ~ 5.0,
    anio %in% c(2008, 2009) ~ 5.5,
    anio %in% c(2010, 2011) ~ 5.0,
    anio == 2012 ~ 4.5,
    anio >= 2013 ~ 4.0,
    TRUE ~ NA_real_
  )
}

# ------------------------------------------------------------
# 5) Verificar variables básicas y reconstruir meta/brechas
# ------------------------------------------------------------
requeridas <- c("fecha", "inf", "exp_12", "exp_24")
faltantes <- setdiff(requeridas, names(bd))

if (length(faltantes) > 0) {
  stop(
    paste0(
      "Faltan variables en la base: ",
      paste(faltantes, collapse = ", ")
    )
  )
}

bd <- bd %>%
  dplyr::mutate(
    meta_central = meta_inflacion_historica(fecha),
    banda_inf    = meta_central - 1,
    banda_sup    = meta_central + 1,

    gap_exp_12 = exp_12 - meta_central,
    gap_exp_24 = exp_24 - meta_central,

    abs_gap_exp_12 = abs(gap_exp_12),
    abs_gap_exp_24 = abs(gap_exp_24),

    d_exp_12 = exp_12 - dplyr::lag(exp_12),
    d_exp_24 = exp_24 - dplyr::lag(exp_24),

    exp12_dentro_banda = dplyr::if_else(
      !is.na(exp_12) & !is.na(banda_inf) &
        exp_12 >= banda_inf & exp_12 <= banda_sup,
      1, 0, missing = NA_real_
    ),
    exp24_dentro_banda = dplyr::if_else(
      !is.na(exp_24) & !is.na(banda_inf) &
        exp_24 >= banda_inf & exp_24 <= banda_sup,
      1, 0, missing = NA_real_
    )
  )

# Verificaciones automáticas
stopifnot(
  isTRUE(all.equal(bd$gap_exp_12, bd$exp_12 - bd$meta_central)),
  isTRUE(all.equal(bd$gap_exp_24, bd$exp_24 - bd$meta_central)),
  isTRUE(all.equal(bd$abs_gap_exp_12, abs(bd$gap_exp_12))),
  isTRUE(all.equal(bd$abs_gap_exp_24, abs(bd$gap_exp_24)))
)

# Indicador de estabilidad
bd <- bd %>%
  dplyr::mutate(
    estabilidad_exp12_12m = zoo::rollapply(
      abs(d_exp_12),
      width = 12,
      FUN = mean,
      align = "right",
      fill = NA,
      na.rm = TRUE
    ),
    estabilidad_exp24_12m = zoo::rollapply(
      abs(d_exp_24),
      width = 12,
      FUN = mean,
      align = "right",
      fill = NA,
      na.rm = TRUE
    )
  )

bd_plot <- bd %>%
  dplyr::filter(fecha >= as.Date("2010-01-01"))

# ------------------------------------------------------------
# 6) Sistema visual del working paper
# ------------------------------------------------------------
col_exp12 <- "#1F3A5F"
col_exp24 <- "#7A8798"
col_meta  <- "#0D47A1"
col_banda <- "#DCE6F1"   # sombra más visible
col_borde <- "#4F86D9"
col_inf   <- "#2C3E50"
col_emei  <- "#B03A2E"
col_cero  <- "#7F8C8D"

tema_wp <- ggplot2::theme_minimal(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(
      face = "plain",
      size = 10.5,
      color = "#333333",
      margin = ggplot2::margin(r = 8)
    ),
    axis.text = ggplot2::element_text(
      size = 9.5,
      color = "#444444"
    ),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E8E8E8",
      linewidth = 0.35
    ),
    panel.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    ),
    plot.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    ),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 9.5),
    legend.key.width = grid::unit(1.35, "cm"),
    legend.margin = ggplot2::margin(t = 5),
    plot.margin = ggplot2::margin(8, 12, 10, 8)
  )

escala_x_wp <- ggplot2::scale_x_date(
  date_breaks = "2 years",
  date_labels = "%Y",
  expand = ggplot2::expansion(mult = c(0.01, 0.015))
)

guardar_figura <- function(plot, nombre, width, height) {
  ggplot2::ggsave(
    filename = file.path(ruta_figuras, paste0(nombre, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = file.path(ruta_figuras, paste0(nombre, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )
}

# ============================================================
# FIGURA DE CONTEXTO
# Evolución histórica de la inflación y metas bajo el EMEI
# ============================================================
ruta_inflacion_hist <- "inflacion_bd.xlsx"

if (file.exists(ruta_inflacion_hist)) {

  infl_hist <- readxl::read_excel(ruta_inflacion_hist, sheet = 1) %>%
    dplyr::mutate(
      fecha = as.Date(paste(anio, mes_num, "01", sep = "-")),
      meta_central = meta_inflacion_historica(fecha),
      banda_inf = meta_central - 1,
      banda_sup = meta_central + 1
    ) %>%
    dplyr::arrange(fecha)

  fecha_emei <- as.Date("2005-01-01")

  g_contexto <- ggplot2::ggplot(
    infl_hist,
    ggplot2::aes(x = fecha, y = inflacion)
  ) +
    ggplot2::geom_ribbon(
      data = infl_hist %>%
        dplyr::filter(fecha >= fecha_emei, !is.na(meta_central)),
      ggplot2::aes(
        ymin = banda_inf,
        ymax = banda_sup,
        group = 1
      ),
      inherit.aes = TRUE,
      fill = col_banda,
      alpha = 0.55,
      na.rm = TRUE
    ) +
    ggplot2::geom_step(
      data = infl_hist %>%
        dplyr::filter(fecha >= fecha_emei, !is.na(meta_central)),
      ggplot2::aes(y = banda_inf),
      inherit.aes = TRUE,
      color = col_borde,
      linetype = "dotted",
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    ggplot2::geom_step(
      data = infl_hist %>%
        dplyr::filter(fecha >= fecha_emei, !is.na(meta_central)),
      ggplot2::aes(y = banda_sup),
      inherit.aes = TRUE,
      color = col_borde,
      linetype = "dotted",
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    ggplot2::geom_step(
      data = infl_hist %>%
        dplyr::filter(fecha >= fecha_emei, !is.na(meta_central)),
      ggplot2::aes(y = meta_central),
      inherit.aes = TRUE,
      color = col_meta,
      linetype = "dashed",
      linewidth = 0.95,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      color = col_inf,
      linewidth = 0.9
    ) +
    ggplot2::geom_vline(
      xintercept = fecha_emei,
      color = col_emei,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    ggplot2::geom_label(
      data = data.frame(
        fecha = as.Date("2005-06-01"),
        y = max(infl_hist$inflacion, na.rm = TRUE) * 0.965,
        texto = "Inicio del EMEI"
      ),
      ggplot2::aes(x = fecha, y = y, label = texto),
      inherit.aes = FALSE,
      size = 3.3,
      color = col_emei,
      fill = "white",
      label.size = 0.15,
      fontface = "bold"
    ) +
    ggplot2::geom_label(
      data = data.frame(
        fecha = as.Date("2016-01-01"),
        y = 8.8,
        texto = "Rango meta de inflación"
      ),
      ggplot2::aes(x = fecha, y = y, label = texto),
      inherit.aes = FALSE,
      size = 3.1,
      color = col_meta,
      fill = "white",
      label.size = 0.15,
      fontface = "bold"
    ) +
    ggplot2::scale_x_date(
      date_breaks = "2 years",
      date_labels = "%Y",
      expand = ggplot2::expansion(mult = c(0.005, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(suffix = "%"),
      expand = ggplot2::expansion(mult = c(0.02, 0.07))
    ) +
    ggplot2::labs(x = NULL, y = "Porcentaje") +
    tema_wp +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  guardar_figura(
    g_contexto,
    "fig_contexto_inflacion_historica_emei",
    width = 10,
    height = 5
  )
}

# ============================================================
# FIGURA 1 DEL OBJETIVO 1
# Expectativas a 12 y 24 meses frente a la meta y rango
# IMPORTANTE: esta versión fuerza el área sombreada y utiliza
# geom_step para que la meta aparezca como escalón.
# ============================================================
g_obj1_expectativas <- ggplot2::ggplot(
  bd_plot,
  ggplot2::aes(x = fecha)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = banda_inf,
      ymax = banda_sup,
      group = 1
    ),
    inherit.aes = TRUE,
    fill = col_banda,
    alpha = 0.55,
    na.rm = TRUE
  ) +
  ggplot2::geom_step(
    ggplot2::aes(y = banda_inf),
    color = col_borde,
    linetype = "dotted",
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  ggplot2::geom_step(
    ggplot2::aes(y = banda_sup),
    color = col_borde,
    linetype = "dotted",
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  ggplot2::geom_step(
    ggplot2::aes(y = meta_central),
    color = col_meta,
    linetype = "dashed",
    linewidth = 0.95,
    na.rm = TRUE
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      y = exp_12,
      color = "Expectativa a 12 meses"
    ),
    linewidth = 1.05,
    na.rm = TRUE
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      y = exp_24,
      color = "Expectativa a 24 meses"
    ),
    linewidth = 1.05,
    na.rm = TRUE
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Expectativa a 12 meses" = col_exp12,
      "Expectativa a 24 meses" = col_exp24
    )
  ) +
  ggplot2::scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y",
    expand = ggplot2::expansion(mult = c(0.01, 0.015))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      suffix = "%",
      accuracy = 0.1
    ),
    limits = c(2.9, 7.0),
    breaks = c(3, 4, 5, 6, 7),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::labs(y = "Porcentaje") +
  ggplot2::annotate(
    "text",
    x = as.Date("2017-08-01"),
    y = 3.78,
    label = "Meta central",
    size = 4.0,
    color = col_meta,
    fontface = "bold"
  ) +
  tema_wp +
  ggplot2::theme(
    legend.position = "bottom"
  )

g_obj1_expectativas

guardar_figura(
  g_obj1_expectativas,
  "fig_obj1_expectativas_meta",
  width = 10,
  height = 5
)

# ============================================================
# FIGURA 2 DEL OBJETIVO 1
# Desviación con signo y distancia absoluta respecto de la meta
# ============================================================
bd_fig2 <- bd_plot %>%
  dplyr::transmute(
    fecha,
    gap_12 = exp_12 - meta_central,
    gap_24 = exp_24 - meta_central,
    distancia_12 = abs(exp_12 - meta_central),
    distancia_24 = abs(exp_24 - meta_central)
  )

bd_gap_exp <- bd_fig2 %>%
  dplyr::select(fecha, gap_12, gap_24) %>%
  tidyr::pivot_longer(
    cols = c(gap_12, gap_24),
    names_to = "serie",
    values_to = "brecha"
  ) %>%
  dplyr::mutate(
    serie = dplyr::recode(
      serie,
      gap_12 = "Expectativa a 12 meses",
      gap_24 = "Expectativa a 24 meses"
    )
  )

bd_dist_abs <- bd_fig2 %>%
  dplyr::select(fecha, distancia_12, distancia_24) %>%
  tidyr::pivot_longer(
    cols = c(distancia_12, distancia_24),
    names_to = "serie",
    values_to = "distancia"
  ) %>%
  dplyr::mutate(
    serie = dplyr::recode(
      serie,
      distancia_12 = "Expectativa a 12 meses",
      distancia_24 = "Expectativa a 24 meses"
    )
  )

g_gap <- ggplot2::ggplot(
  bd_gap_exp,
  ggplot2::aes(
    x = fecha,
    y = brecha,
    color = serie
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = col_cero,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::scale_color_manual(
    values = c(
      "Expectativa a 12 meses" = col_exp12,
      "Expectativa a 24 meses" = col_exp24
    )
  ) +
  escala_x_wp +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      suffix = " p.p.",
      accuracy = 0.1
    ),
    expand = ggplot2::expansion(mult = c(0.05, 0.07))
  ) +
  ggplot2::labs(y = "Desviación respecto de la meta") +
  tema_wp +
  ggplot2::theme(
    legend.position = "none",
    plot.margin = ggplot2::margin(7, 12, 2, 8)
  )

g_distancia <- ggplot2::ggplot(
  bd_dist_abs,
  ggplot2::aes(
    x = fecha,
    y = distancia,
    color = serie
  )
) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::scale_color_manual(
    values = c(
      "Expectativa a 12 meses" = col_exp12,
      "Expectativa a 24 meses" = col_exp24
    )
  ) +
  escala_x_wp +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      suffix = " p.p.",
      accuracy = 0.1
    ),
    limits = c(0, NA),
    expand = ggplot2::expansion(mult = c(0, 0.07))
  ) +
  ggplot2::labs(y = "Distancia respecto de la meta") +
  tema_wp +
  ggplot2::theme(
    plot.margin = ggplot2::margin(2, 12, 7, 8)
  )

g_obj1_desviaciones <- (
  g_gap /
    g_distancia +
    patchwork::plot_layout(
      heights = c(1, 1),
      guides = "collect"
    ) +
    patchwork::plot_annotation(tag_levels = "A")
) &
  ggplot2::theme(
    legend.position = "bottom",
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 10.5,
      color = "#333333"
    )
  )

g_obj1_desviaciones

guardar_figura(
  g_obj1_desviaciones,
  "fig_obj1_desviacion_distancia_meta",
  width = 10,
  height = 7.2
)

# ============================================================
# FIGURA 3 DEL OBJETIVO 1
# Estabilidad temporal de las expectativas
# ============================================================
bd_estabilidad <- bd_plot %>%
  dplyr::select(
    fecha,
    estabilidad_exp12_12m,
    estabilidad_exp24_12m
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      estabilidad_exp12_12m,
      estabilidad_exp24_12m
    ),
    names_to = "serie",
    values_to = "variacion_abs_media"
  ) %>%
  dplyr::mutate(
    serie = dplyr::recode(
      serie,
      estabilidad_exp12_12m = "Expectativa a 12 meses",
      estabilidad_exp24_12m = "Expectativa a 24 meses"
    )
  )

g_obj1_estabilidad <- ggplot2::ggplot(
  bd_estabilidad,
  ggplot2::aes(
    x = fecha,
    y = variacion_abs_media,
    color = serie
  )
) +
  ggplot2::geom_line(linewidth = 0.95) +
  ggplot2::scale_color_manual(
    values = c(
      "Expectativa a 12 meses" = col_exp12,
      "Expectativa a 24 meses" = col_exp24
    )
  ) +
  escala_x_wp +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      suffix = " p.p.",
      accuracy = 0.01
    ),
    expand = ggplot2::expansion(mult = c(0.03, 0.07))
  ) +
  ggplot2::labs(y = "Variación mensual absoluta promedio") +
  tema_wp

g_obj1_estabilidad

guardar_figura(
  g_obj1_estabilidad,
  "fig_obj1_estabilidad_expectativas",
  width = 10,
  height = 4.8
)

# ============================================================
# CUADRO 1 DEL OBJETIVO 1
# Estadísticos descriptivos
# ============================================================
bd_cuadro <- bd_plot %>%
  dplyr::mutate(
    subperiodo = dplyr::case_when(
      fecha >= as.Date("2010-01-01") &
        fecha <= as.Date("2015-12-01") ~ "2010-2015",
      fecha >= as.Date("2016-01-01") &
        fecha <= as.Date("2020-12-01") ~ "2016-2020",
      fecha >= as.Date("2021-01-01") ~ "2021-2026",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(subperiodo))

resumen_obj1 <- function(data, etiqueta_periodo) {
  tibble::tibble(
    Periodo = etiqueta_periodo,
    exp12_media = mean(data$exp_12, na.rm = TRUE),
    exp24_media = mean(data$exp_24, na.rm = TRUE),
    exp12_sd = stats::sd(data$exp_12, na.rm = TRUE),
    exp24_sd = stats::sd(data$exp_24, na.rm = TRUE),
    dist12_media = mean(data$abs_gap_exp_12, na.rm = TRUE),
    dist24_media = mean(data$abs_gap_exp_24, na.rm = TRUE),
    cambio_abs12_media = mean(abs(data$d_exp_12), na.rm = TRUE),
    cambio_abs24_media = mean(abs(data$d_exp_24), na.rm = TRUE),
    exp12_dentro_banda = 100 * mean(data$exp12_dentro_banda, na.rm = TRUE),
    exp24_dentro_banda = 100 * mean(data$exp24_dentro_banda, na.rm = TRUE)
  )
}

cuadro1_wide <- dplyr::bind_rows(
  resumen_obj1(bd_cuadro, "2010-2026"),
  resumen_obj1(
    bd_cuadro %>% dplyr::filter(subperiodo == "2010-2015"),
    "2010-2015"
  ),
  resumen_obj1(
    bd_cuadro %>% dplyr::filter(subperiodo == "2016-2020"),
    "2016-2020"
  ),
  resumen_obj1(
    bd_cuadro %>% dplyr::filter(subperiodo == "2021-2026"),
    "2021-2026"
  )
)

cuadro1_long <- cuadro1_wide %>%
  tidyr::pivot_longer(
    cols = -Periodo,
    names_to = "estadistico",
    values_to = "valor"
  ) %>%
  dplyr::mutate(
    bloque = dplyr::case_when(
      estadistico %in% c("exp12_media", "exp24_media") ~
        "Nivel promedio",
      estadistico %in% c("exp12_sd", "exp24_sd") ~
        "Desviación estándar",
      estadistico %in% c("dist12_media", "dist24_media") ~
        "Distancia absoluta media respecto de la meta",
      estadistico %in% c("cambio_abs12_media", "cambio_abs24_media") ~
        "Variación mensual absoluta media",
      estadistico %in% c("exp12_dentro_banda", "exp24_dentro_banda") ~
        "Observaciones dentro del rango meta (%)",
      TRUE ~ "Otros"
    ),
    horizonte = dplyr::case_when(
      stringr::str_detect(estadistico, "12") ~ "12 meses",
      stringr::str_detect(estadistico, "24") ~ "24 meses",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(
    bloque,
    horizonte,
    Periodo,
    valor
  )

cuadro1_final <- cuadro1_long %>%
  tidyr::pivot_wider(
    names_from = Periodo,
    values_from = valor
  ) %>%
  dplyr::mutate(
    dplyr::across(
      where(is.numeric),
      ~ round(.x, 2)
    )
  )

write.csv(
  cuadro1_final,
  file = file.path(
    ruta_tablas,
    "cuadro_obj1_expectativas.csv"
  ),
  row.names = FALSE
)

saveRDS(
  cuadro1_final,
  file = file.path(
    ruta_tablas,
    "cuadro_obj1_expectativas.rds"
  )
)

print(cuadro1_final)

# ------------------------------------------------------------
# 7) Guardar base enriquecida
# ------------------------------------------------------------
saveRDS(
  bd_plot,
  file = file.path(
    ruta_salida,
    "bd_objetivo1_figuras.rds"
  )
)

# ------------------------------------------------------------
# 8) Resumen en consola
# ------------------------------------------------------------
cat("\n============================================================\n")
cat("OBJETIVO 1 - FIGURAS Y CUADRO GENERADOS\n")
cat("============================================================\n\n")
cat("La FIGURA 1 y la FIGURA DE CONTEXTO se generan con área sombreada.\n\n")
cat("FIGURAS PRINCIPALES:\n")
cat("1. fig_obj1_expectativas_meta\n")
cat("2. fig_obj1_desviacion_distancia_meta\n")
cat("3. fig_obj1_estabilidad_expectativas\n\n")
cat("FIGURA DE CONTEXTO:\n")
cat("- fig_contexto_inflacion_historica_emei\n\n")
cat("CUADRO PRINCIPAL:\n")
cat("- cuadro_obj1_expectativas.csv\n\n")
cat("Directorio de figuras: ", ruta_figuras, "\n")
cat("Directorio de tablas: ", ruta_tablas, "\n")
