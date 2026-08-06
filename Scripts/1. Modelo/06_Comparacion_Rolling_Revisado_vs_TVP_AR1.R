# ==============================================================================
# OBJETIVO 2 - SCRIPT 06 CORREGIDO
# COMPARACIÓN ENTRE ROLLING REGRESSIONS REVISADAS Y TVP-KALMAN AR(1)
#
# ROLLING PRINCIPAL (ventana de 36 meses):
#
#   (e24_t - meta_t) = alpha_w
#                    + rho_w  (e24_(t-1) - meta_(t-1))
#                    + beta_w (pi_(t-1)  - meta_(t-1))
#                    + gamma_w q_(t-1)
#                    + u_t
#
# Robustez rolling: ventana de 30 meses.
# Inferencia rolling: errores estándar Newey-West HAC, rezago 3.
#
# Medidas rolling:
#   - beta_w: sensibilidad de corto plazo a la brecha de inflación.
#   - rho_w: persistencia de la brecha de expectativas.
#   - gamma_w: respuesta a la depreciación cambiaria.
#   - beta_LP,w = beta_w / (1 - rho_w): sensibilidad de largo plazo.
#
# TVP PRINCIPAL:
#
#   Delta e24_t = alpha
#                 + rho Delta e24_(t-1)
#                 + beta_t shock_inf_t
#                 + epsilon_t
#
# TVP DE ROBUSTEZ:
#   Misma ecuación con pulso en septiembre de 2012.
#
# ADVERTENCIA:
# Los coeficientes rolling y TVP no estiman exactamente el mismo parámetro:
# cambian la variable dependiente, la definición del impulso inflacionario y
# la forma de incorporar la dinámica. Por eso la comparación principal se
# concentra en trayectoria, tendencia, comovimiento, primeras diferencias,
# promedios anuales e índices estandarizados; no en igualdad de magnitudes.
# ==============================================================================


# 0. PAQUETES ------------------------------------------------------------------

paquetes <- c(
  "dplyr",
  "tidyr",
  "readr",
  "lubridate",
  "ggplot2",
  "purrr",
  "tibble",
  "sandwich",
  "lmtest"
)

faltantes <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(faltantes) > 0) {
  install.packages(faltantes, dependencies = TRUE)
}

invisible(lapply(paquetes, library, character.only = TRUE))


# 1. CONFIGURACIÓN -------------------------------------------------------------

archivos_bd_candidatos <- c(
  "bd_sample_modelo.csv",
  "bd_sample_modelo(1).csv"
)

archivos_tvp_base_candidatos <- c(
  file.path(
    "resultados_objetivo2_constante_vs_tvp_ar1",
    "04_trayectoria_beta_tvp.csv"
  ),
  file.path(
    "resultados_objetivo2_tvp_ar1",
    "03_resultados_tvp_kalman_ar1.csv"
  ),
  "04_trayectoria_beta_tvp.csv"
)

archivos_tvp_pulso_candidatos <- c(
  file.path(
    "resultados_objetivo2_constante_vs_tvp_ar1_pulso_sep2012",
    "04_trayectoria_beta_tvp_pulso.csv"
  ),
  file.path(
    "resultados_objetivo2_tvp_ar1_dummy_sep2012",
    "04_resultados_modelo_dummy.csv"
  ),
  "04_trayectoria_beta_tvp_pulso.csv"
)

archivos_rolling_36_previos <- c(
  file.path("output", "rolling_modelo_principal_36m.csv"),
  "rolling_modelo_principal_36m.csv"
)

archivos_rolling_30_previos <- c(
  file.path("output", "rolling_robustez_30m.csv"),
  "rolling_robustez_30m.csv"
)

carpeta_salida <- file.path(
  getwd(),
  "resultados_objetivo2_comparacion_rolling_revisado_tvp"
)

col_fecha <- "fecha"
col_exp24 <- "exp_inf_24m"
col_inflacion <- "infl_gt"

candidatas_tc <- c(
  "tcdep_var",
  "depreciacion_tc",
  "tc_var",
  "q_t",
  "tipo_cambio_var",
  "dep_tc"
)

ventana_principal <- 36L
ventana_robustez <- 30L
hac_lag <- 3L
nivel_confianza <- 0.95
z_critico <- qnorm(1 - (1 - nivel_confianza) / 2)

dir.create(
  carpeta_salida,
  showWarnings = FALSE,
  recursive = TRUE
)

if (!dir.exists(carpeta_salida)) {
  stop("No fue posible crear la carpeta: ", carpeta_salida)
}

if (file.access(carpeta_salida, 2) != 0) {
  stop("R no tiene permiso de escritura en: ", carpeta_salida)
}


# 2. FUNCIONES GENERALES -------------------------------------------------------

encontrar_archivo <- function(
    candidatos,
    obligatorio = TRUE,
    descripcion = "archivo"
) {

  encontrados <- candidatos[file.exists(candidatos)]

  if (length(encontrados) == 0) {

    if (obligatorio) {
      stop(
        paste0(
          "No se encontró ", descripcion, ".\n",
          "Rutas revisadas:\n",
          paste0("- ", candidatos, collapse = "\n"),
          "\n\nDirectorio actual: ",
          getwd()
        )
      )
    }

    return(NA_character_)
  }

  encontrados[[1]]
}


primera_columna_disponible <- function(datos, candidatas) {

  disponibles <- candidatas[candidatas %in% names(datos)]

  if (length(disponibles) == 0) {
    return(NA_character_)
  }

  disponibles[[1]]
}


verificar_columnas <- function(datos, columnas) {

  faltan <- setdiff(columnas, names(datos))

  if (length(faltan) > 0) {
    stop(
      "No se encontraron estas columnas: ",
      paste(faltan, collapse = ", ")
    )
  }
}


estandarizar <- function(x) {

  desv <- sd(x, na.rm = TRUE)

  if (!is.finite(desv) || desv == 0) {
    return(rep(NA_real_, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) / desv
}


normalizar_01 <- function(x) {

  rango <- range(x, na.rm = TRUE)
  amplitud <- diff(rango)

  if (
    length(amplitud) == 0 ||
      any(!is.finite(rango)) ||
      !is.finite(amplitud) ||
      amplitud == 0
  ) {
    return(rep(NA_real_, length(x)))
  }

  (x - rango[1]) / amplitud
}


correlacion_segura <- function(x, y, metodo = "pearson") {

  validos <- is.finite(x) & is.finite(y)

  if (sum(validos) < 3) {
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      x[validos],
      y[validos],
      method = metodo
    )
  )
}


pendiente_anual <- function(fecha, serie) {

  datos <- tibble(
    fecha = fecha,
    serie = serie
  ) |>
    filter(
      !is.na(fecha),
      is.finite(serie)
    )

  if (nrow(datos) < 3) {
    return(NA_real_)
  }

  modelo <- lm(
    serie ~ as.numeric(fecha),
    data = datos
  )

  unname(coef(modelo)[2]) * 365.25
}


clasificar_tendencia <- function(pendiente) {

  if (!is.finite(pendiente)) {
    return("No disponible")
  }

  if (pendiente < 0) {
    return("Descendente")
  }

  if (pendiente > 0) {
    return("Ascendente")
  }

  "Sin cambio"
}


signo_cambio <- function(x, tolerancia = 1e-10) {

  case_when(
    !is.finite(x) ~ NA_character_,
    x > tolerancia ~ "Aumento de sensibilidad",
    x < -tolerancia ~ "Reducción de sensibilidad",
    TRUE ~ "Sin cambio"
  )
}


# 3. LOCALIZACIÓN DE ARCHIVOS --------------------------------------------------

archivo_bd <- encontrar_archivo(
  archivos_bd_candidatos,
  obligatorio = TRUE,
  descripcion = "la base de datos"
)

archivo_tvp_base <- encontrar_archivo(
  archivos_tvp_base_candidatos,
  obligatorio = TRUE,
  descripcion = "la trayectoria TVP base"
)

archivo_tvp_pulso <- encontrar_archivo(
  archivos_tvp_pulso_candidatos,
  obligatorio = FALSE,
  descripcion = "la trayectoria TVP con pulso"
)

archivo_rolling_36_previo <- encontrar_archivo(
  archivos_rolling_36_previos,
  obligatorio = FALSE,
  descripcion = "el resultado rolling previo de 36 meses"
)

archivo_rolling_30_previo <- encontrar_archivo(
  archivos_rolling_30_previos,
  obligatorio = FALSE,
  descripcion = "el resultado rolling previo de 30 meses"
)


# 4. PREPARACIÓN DE LA BASE ROLLING -------------------------------------------

df_model <- readr::read_csv(
  archivo_bd,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", ".", "null")
)

verificar_columnas(
  df_model,
  c(
    col_fecha,
    col_exp24,
    col_inflacion
  )
)

col_tc <- primera_columna_disponible(
  df_model,
  candidatas_tc
)

if (is.na(col_tc)) {
  stop(
    "No se encontró una variable de depreciación cambiaria. ",
    "Agregue su nombre a 'candidatas_tc'. Columnas disponibles: ",
    paste(names(df_model), collapse = ", ")
  )
}

message(
  "Variable de tipo de cambio seleccionada: ",
  col_tc
)

df_model <- df_model |>
  mutate(
    fecha = lubridate::ymd(.data[[col_fecha]])
  ) |>
  arrange(fecha)

if (anyNA(df_model$fecha)) {
  stop("Existen fechas que no pudieron convertirse al formato Date.")
}

# Se replica exactamente la construcción del script rolling revisado.
if (!"meta_inflacion" %in% names(df_model)) {

  df_model <- df_model |>
    mutate(
      meta_inflacion = case_when(
        fecha < ymd("2012-01-01") ~ 5.0,
        fecha < ymd("2013-01-01") ~ 4.5,
        TRUE ~ 4.0
      )
    )

} else {

  message(
    "Se utilizará la columna existente 'meta_inflacion'."
  )
}

df_est <- df_model |>
  transmute(
    fecha,
    exp_inf_24m = .data[[col_exp24]],
    infl_gt = .data[[col_inflacion]],
    meta_inflacion,
    tcdep = .data[[col_tc]],

    exp_gap_24 =
      exp_inf_24m - meta_inflacion,

    infl_gap =
      infl_gt - meta_inflacion,

    exp_gap_24_lag1 =
      lag(exp_gap_24, 1),

    infl_gap_lag1 =
      lag(infl_gap, 1),

    tcdep_lag1 =
      lag(tcdep, 1)
  ) |>
  tidyr::drop_na(
    fecha,
    exp_gap_24,
    exp_gap_24_lag1,
    infl_gap_lag1,
    tcdep_lag1
  )

if (nrow(df_est) < ventana_principal) {
  stop(
    "La muestra efectiva tiene ",
    nrow(df_est),
    " observaciones, menos que la ventana principal de ",
    ventana_principal,
    " meses."
  )
}

saltos_meses <- diff(
  year(df_est$fecha) * 12 +
    month(df_est$fecha)
)

if (any(saltos_meses != 1)) {
  warning(
    "La muestra efectiva presenta saltos entre meses. ",
    "Revise observaciones faltantes antes de interpretar las ventanas."
  )
}


# 5. ESTIMACIÓN ROLLING REVISADA ----------------------------------------------

formula_principal <- exp_gap_24 ~
  exp_gap_24_lag1 +
  infl_gap_lag1 +
  tcdep_lag1

nombre_rho <- "exp_gap_24_lag1"
nombre_beta <- "infl_gap_lag1"
nombre_gamma <- "tcdep_lag1"


extraer_coeficiente_hac <- function(
    tabla_hac,
    nombre
) {

  if (!nombre %in% rownames(tabla_hac)) {
    return(
      c(
        estimate = NA_real_,
        se = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_
      )
    )
  }

  c(
    estimate =
      tabla_hac[nombre, "Estimate"],

    se =
      tabla_hac[nombre, "Std. Error"],

    statistic =
      tabla_hac[nombre, "t value"],

    p_value =
      tabla_hac[nombre, "Pr(>|t|)"]
  )
}


calcular_largo_plazo <- function(
    coeficientes,
    vcov_hac,
    nombre_rho,
    nombre_beta
) {

  rho <- unname(coeficientes[nombre_rho])
  beta <- unname(coeficientes[nombre_beta])

  denominador <- 1 - rho

  if (
    !is.finite(rho) ||
      !is.finite(beta) ||
      !is.finite(denominador) ||
      abs(denominador) < 0.05
  ) {
    return(
      c(
        efecto_lp = NA_real_,
        se_lp = NA_real_,
        ic_inf_lp = NA_real_,
        ic_sup_lp = NA_real_
      )
    )
  }

  efecto_lp <- beta / denominador

  nombres_necesarios <- c(
    nombre_rho,
    nombre_beta
  )

  if (
    !all(nombres_necesarios %in% rownames(vcov_hac)) ||
      !all(nombres_necesarios %in% colnames(vcov_hac))
  ) {
    return(
      c(
        efecto_lp = efecto_lp,
        se_lp = NA_real_,
        ic_inf_lp = NA_real_,
        ic_sup_lp = NA_real_
      )
    )
  }

  gradiente <- matrix(
    c(
      beta / (1 - rho)^2,
      1 / (1 - rho)
    ),
    ncol = 1
  )

  V_sub <- vcov_hac[
    nombres_necesarios,
    nombres_necesarios,
    drop = FALSE
  ]

  var_lp <- as.numeric(
    t(gradiente) %*%
      V_sub %*%
      gradiente
  )

  se_lp <- ifelse(
    is.finite(var_lp) && var_lp >= 0,
    sqrt(var_lp),
    NA_real_
  )

  c(
    efecto_lp = efecto_lp,
    se_lp = se_lp,
    ic_inf_lp = efecto_lp - z_critico * se_lp,
    ic_sup_lp = efecto_lp + z_critico * se_lp
  )
}


estimar_ventana <- function(
    data_window,
    window_size,
    hac_lag = 3L
) {

  modelo <- lm(
    formula_principal,
    data = data_window
  )

  vcov_hac <- sandwich::NeweyWest(
    modelo,
    lag = hac_lag,
    prewhite = FALSE,
    adjust = TRUE
  )

  tabla_hac <- lmtest::coeftest(
    modelo,
    vcov. = vcov_hac
  )

  coeficientes <- coef(modelo)

  alpha <- extraer_coeficiente_hac(
    tabla_hac,
    "(Intercept)"
  )

  rho <- extraer_coeficiente_hac(
    tabla_hac,
    nombre_rho
  )

  beta <- extraer_coeficiente_hac(
    tabla_hac,
    nombre_beta
  )

  gamma <- extraer_coeficiente_hac(
    tabla_hac,
    nombre_gamma
  )

  largo_plazo <- calcular_largo_plazo(
    coeficientes = coeficientes,
    vcov_hac = vcov_hac,
    nombre_rho = nombre_rho,
    nombre_beta = nombre_beta
  )

  tibble(
    fecha_inicio = min(data_window$fecha),
    fecha_fin = max(data_window$fecha),
    ventana = window_size,
    n_obs = nobs(modelo),

    alpha = unname(alpha["estimate"]),
    alpha_se_hac = unname(alpha["se"]),
    alpha_p_value = unname(alpha["p_value"]),

    rho = unname(rho["estimate"]),
    rho_se_hac = unname(rho["se"]),
    rho_p_value = unname(rho["p_value"]),
    rho_estable = abs(unname(rho["estimate"])) < 1,

    beta_corto_plazo = unname(beta["estimate"]),
    beta_se_hac = unname(beta["se"]),
    beta_t_hac = unname(beta["statistic"]),
    beta_p_value = unname(beta["p_value"]),
    beta_ic_inf = unname(
      beta["estimate"] -
        z_critico * beta["se"]
    ),
    beta_ic_sup = unname(
      beta["estimate"] +
        z_critico * beta["se"]
    ),

    gamma_tc = unname(gamma["estimate"]),
    gamma_se_hac = unname(gamma["se"]),
    gamma_p_value = unname(gamma["p_value"]),

    beta_largo_plazo =
      unname(largo_plazo["efecto_lp"]),
    beta_largo_plazo_se =
      unname(largo_plazo["se_lp"]),
    beta_largo_plazo_ic_inf =
      unname(largo_plazo["ic_inf_lp"]),
    beta_largo_plazo_ic_sup =
      unname(largo_plazo["ic_sup_lp"]),

    r_squared = summary(modelo)$r.squared,
    adj_r_squared = summary(modelo)$adj.r.squared
  )
}


correr_rolling <- function(
    datos,
    window_size,
    hac_lag = 3L
) {

  if (nrow(datos) < window_size) {
    stop(
      "La ventana solicitada excede la muestra disponible."
    )
  }

  purrr::map_dfr(
    seq.int(window_size, nrow(datos)),
    function(i) {

      data_window <- datos[
        (i - window_size + 1):i,
        ,
        drop = FALSE
      ]

      estimar_ventana(
        data_window = data_window,
        window_size = window_size,
        hac_lag = hac_lag
      )
    }
  )
}


rolling_36 <- correr_rolling(
  datos = df_est,
  window_size = ventana_principal,
  hac_lag = hac_lag
)

rolling_30 <- correr_rolling(
  datos = df_est,
  window_size = ventana_robustez,
  hac_lag = hac_lag
)

readr::write_csv(
  rolling_36,
  file.path(
    carpeta_salida,
    "01_rolling_revisado_completo_36m.csv"
  )
)

readr::write_csv(
  rolling_30,
  file.path(
    carpeta_salida,
    "02_rolling_revisado_completo_30m.csv"
  )
)


# 6. VERIFICACIÓN CONTRA LOS RESULTADOS PREVIOS -------------------------------

verificar_rolling_previo <- function(
    archivo_previo,
    rolling_nuevo,
    ventana
) {

  if (is.na(archivo_previo)) {
    return(
      tibble(
        ventana = ventana,
        archivo_previo_encontrado = FALSE,
        archivo_previo = "",
        observaciones_comparadas = 0L,
        diferencia_absoluta_media = NA_real_,
        diferencia_absoluta_maxima = NA_real_,
        coincide_tolerancia_1e_8 = NA
      )
    )
  }

  previo <- readr::read_csv(
    archivo_previo,
    show_col_types = FALSE
  )

  if (
    !all(c("fecha_fin", "estimate") %in% names(previo))
  ) {
    return(
      tibble(
        ventana = ventana,
        archivo_previo_encontrado = TRUE,
        archivo_previo = archivo_previo,
        observaciones_comparadas = 0L,
        diferencia_absoluta_media = NA_real_,
        diferencia_absoluta_maxima = NA_real_,
        coincide_tolerancia_1e_8 = FALSE
      )
    )
  }

  comparacion <- rolling_nuevo |>
    transmute(
      fecha_fin,
      beta_nuevo = beta_corto_plazo
    ) |>
    inner_join(
      previo |>
        transmute(
          fecha_fin = ymd(fecha_fin),
          beta_previo = estimate
        ),
      by = "fecha_fin"
    ) |>
    mutate(
      diferencia_absoluta =
        abs(beta_nuevo - beta_previo)
    )

  if (nrow(comparacion) == 0) {
    return(
      tibble(
        ventana = ventana,
        archivo_previo_encontrado = TRUE,
        archivo_previo = archivo_previo,
        observaciones_comparadas = 0L,
        diferencia_absoluta_media = NA_real_,
        diferencia_absoluta_maxima = NA_real_,
        coincide_tolerancia_1e_8 = FALSE
      )
    )
  }

  diferencia_maxima <- max(
    comparacion$diferencia_absoluta,
    na.rm = TRUE
  )

  tibble(
    ventana = ventana,
    archivo_previo_encontrado = TRUE,
    archivo_previo = archivo_previo,
    observaciones_comparadas = nrow(comparacion),
    diferencia_absoluta_media = mean(
      comparacion$diferencia_absoluta,
      na.rm = TRUE
    ),
    diferencia_absoluta_maxima = diferencia_maxima,
    coincide_tolerancia_1e_8 =
      diferencia_maxima < 1e-8
  )
}


tabla_verificacion <- bind_rows(
  verificar_rolling_previo(
    archivo_previo = archivo_rolling_36_previo,
    rolling_nuevo = rolling_36,
    ventana = 36L
  ),
  verificar_rolling_previo(
    archivo_previo = archivo_rolling_30_previo,
    rolling_nuevo = rolling_30,
    ventana = 30L
  )
)

readr::write_csv(
  tabla_verificacion,
  file.path(
    carpeta_salida,
    "03_verificacion_resultados_rolling_previos.csv"
  )
)


# 7. LECTURA DE TRAYECTORIAS TVP ----------------------------------------------

extraer_tvp <- function(
    archivo,
    nombre_beta,
    nombre_li,
    nombre_ls
) {

  datos <- readr::read_csv(
    archivo,
    show_col_types = FALSE,
    na = c("", "NA", "N/A", ".", "null")
  )

  verificar_columnas(
    datos,
    c("fecha", "beta_t")
  )

  salida <- datos |>
    transmute(
      fecha = ymd(fecha),
      "{nombre_beta}" := beta_t
    )

  if ("beta_li_95" %in% names(datos)) {
    salida[[nombre_li]] <- datos$beta_li_95
  } else {
    salida[[nombre_li]] <- NA_real_
  }

  if ("beta_ls_95" %in% names(datos)) {
    salida[[nombre_ls]] <- datos$beta_ls_95
  } else {
    salida[[nombre_ls]] <- NA_real_
  }

  salida |>
    arrange(fecha)
}


tvp_base <- extraer_tvp(
  archivo = archivo_tvp_base,
  nombre_beta = "beta_tvp_base",
  nombre_li = "beta_tvp_base_li_95",
  nombre_ls = "beta_tvp_base_ls_95"
)

tvp_pulso_disponible <- !is.na(archivo_tvp_pulso)

if (tvp_pulso_disponible) {

  tvp_pulso <- extraer_tvp(
    archivo = archivo_tvp_pulso,
    nombre_beta = "beta_tvp_pulso",
    nombre_li = "beta_tvp_pulso_li_95",
    nombre_ls = "beta_tvp_pulso_ls_95"
  )

} else {

  tvp_pulso <- tibble(
    fecha = as.Date(character()),
    beta_tvp_pulso = numeric(),
    beta_tvp_pulso_li_95 = numeric(),
    beta_tvp_pulso_ls_95 = numeric()
  )
}


# 8. ALINEACIÓN TEMPORAL -------------------------------------------------------

rolling_36_fecha <- rolling_36 |>
  transmute(
    fecha = fecha_fin,
    rho_rolling_36 = rho,
    beta_rolling_36 = beta_corto_plazo,
    beta_rolling_36_li_95 = beta_ic_inf,
    beta_rolling_36_ls_95 = beta_ic_sup,
    beta_lp_rolling_36 = beta_largo_plazo,
    beta_lp_rolling_36_li_95 =
      beta_largo_plazo_ic_inf,
    beta_lp_rolling_36_ls_95 =
      beta_largo_plazo_ic_sup,
    gamma_tc_rolling_36 = gamma_tc
  )

rolling_30_fecha <- rolling_30 |>
  transmute(
    fecha = fecha_fin,
    rho_rolling_30 = rho,
    beta_rolling_30 = beta_corto_plazo,
    beta_lp_rolling_30 = beta_largo_plazo,
    gamma_tc_rolling_30 = gamma_tc
  )

comparacion <- rolling_36_fecha |>
  inner_join(
    tvp_base,
    by = "fecha"
  ) |>
  left_join(
    rolling_30_fecha,
    by = "fecha"
  )

if (tvp_pulso_disponible) {
  comparacion <- comparacion |>
    left_join(
      tvp_pulso,
      by = "fecha"
    )
}

if (nrow(comparacion) < 12) {
  stop(
    "La muestra común entre rolling revisado y TVP ",
    "tiene menos de 12 observaciones."
  )
}


# 9. TRANSFORMACIONES PARA COMPARACIÓN ----------------------------------------

comparacion <- comparacion |>
  arrange(fecha) |>
  mutate(
    z_beta_rolling_36 =
      estandarizar(beta_rolling_36),

    z_beta_lp_rolling_36 =
      estandarizar(beta_lp_rolling_36),

    z_beta_tvp_base =
      estandarizar(beta_tvp_base),

    indice_anclaje_rolling_36 =
      1 - normalizar_01(beta_rolling_36),

    indice_anclaje_lp_rolling_36 =
      1 - normalizar_01(beta_lp_rolling_36),

    indice_anclaje_tvp_base =
      1 - normalizar_01(beta_tvp_base),

    d_beta_rolling_36 =
      beta_rolling_36 -
        lag(beta_rolling_36),

    d_beta_lp_rolling_36 =
      beta_lp_rolling_36 -
        lag(beta_lp_rolling_36),

    d_beta_tvp_base =
      beta_tvp_base -
        lag(beta_tvp_base)
  )

if ("beta_rolling_30" %in% names(comparacion)) {

  comparacion <- comparacion |>
    mutate(
      z_beta_rolling_30 =
        estandarizar(beta_rolling_30),

      d_beta_rolling_30 =
        beta_rolling_30 -
          lag(beta_rolling_30)
    )
}

if (
  tvp_pulso_disponible &&
    "beta_tvp_pulso" %in% names(comparacion)
) {

  comparacion <- comparacion |>
    mutate(
      z_beta_tvp_pulso =
        estandarizar(beta_tvp_pulso),

      indice_anclaje_tvp_pulso =
        1 - normalizar_01(beta_tvp_pulso),

      d_beta_tvp_pulso =
        beta_tvp_pulso -
          lag(beta_tvp_pulso)
    )
}

readr::write_csv(
  comparacion,
  file.path(
    carpeta_salida,
    "04_comparacion_mensual_rolling_revisado_tvp.csv"
  )
)


# 10. CORRELACIONES ------------------------------------------------------------

crear_fila_correlacion <- function(
    serie_1_nombre,
    serie_2_nombre,
    serie_1,
    serie_2,
    diferencia_1,
    diferencia_2,
    tipo_comparacion
) {

  validos_nivel <-
    is.finite(serie_1) &
    is.finite(serie_2)

  validos_cambio <-
    is.finite(diferencia_1) &
    is.finite(diferencia_2)

  tibble(
    tipo_comparacion = tipo_comparacion,
    serie_1 = serie_1_nombre,
    serie_2 = serie_2_nombre,
    observaciones_nivel =
      sum(validos_nivel),
    correlacion_pearson_nivel =
      correlacion_segura(
        serie_1,
        serie_2,
        "pearson"
      ),
    correlacion_spearman_nivel =
      correlacion_segura(
        serie_1,
        serie_2,
        "spearman"
      ),
    observaciones_primeras_diferencias =
      sum(validos_cambio),
    correlacion_pearson_primeras_diferencias =
      correlacion_segura(
        diferencia_1,
        diferencia_2,
        "pearson"
      ),
    correlacion_spearman_primeras_diferencias =
      correlacion_segura(
        diferencia_1,
        diferencia_2,
        "spearman"
      )
  )
}


lista_correlaciones <- list(
  crear_fila_correlacion(
    "Rolling 36m: beta corto plazo",
    "TVP-AR(1) base",
    comparacion$beta_rolling_36,
    comparacion$beta_tvp_base,
    comparacion$d_beta_rolling_36,
    comparacion$d_beta_tvp_base,
    "Comparación metodológica principal"
  ),

  crear_fila_correlacion(
    "Rolling 36m: beta largo plazo",
    "TVP-AR(1) base",
    comparacion$beta_lp_rolling_36,
    comparacion$beta_tvp_base,
    comparacion$d_beta_lp_rolling_36,
    comparacion$d_beta_tvp_base,
    "Comparación complementaria"
  ),

  crear_fila_correlacion(
    "Rolling 36m",
    "Rolling 30m",
    comparacion$beta_rolling_36,
    comparacion$beta_rolling_30,
    comparacion$d_beta_rolling_36,
    comparacion$d_beta_rolling_30,
    "Robustez de ventana rolling"
  )
)

if (
  tvp_pulso_disponible &&
    "beta_tvp_pulso" %in% names(comparacion)
) {

  lista_correlaciones <- c(
    lista_correlaciones,
    list(
      crear_fila_correlacion(
        "Rolling 36m: beta corto plazo",
        "TVP-AR(1) con pulso",
        comparacion$beta_rolling_36,
        comparacion$beta_tvp_pulso,
        comparacion$d_beta_rolling_36,
        comparacion$d_beta_tvp_pulso,
        "Robustez TVP"
      ),

      crear_fila_correlacion(
        "TVP-AR(1) base",
        "TVP-AR(1) con pulso",
        comparacion$beta_tvp_base,
        comparacion$beta_tvp_pulso,
        comparacion$d_beta_tvp_base,
        comparacion$d_beta_tvp_pulso,
        "Robustez interna TVP"
      )
    )
  )
}

tabla_correlaciones <- bind_rows(
  lista_correlaciones
)

readr::write_csv(
  tabla_correlaciones,
  file.path(
    carpeta_salida,
    "05_correlaciones_trayectorias.csv"
  )
)


# 11. PROMEDIOS ANUALES --------------------------------------------------------

comparacion_anual <- comparacion |>
  mutate(
    anio = year(fecha)
  ) |>
  group_by(anio) |>
  summarise(
    observaciones = n(),
    beta_rolling_36 =
      mean(beta_rolling_36, na.rm = TRUE),
    beta_lp_rolling_36 =
      mean(beta_lp_rolling_36, na.rm = TRUE),
    beta_tvp_base =
      mean(beta_tvp_base, na.rm = TRUE),
    beta_rolling_30 =
      mean(beta_rolling_30, na.rm = TRUE),
    .groups = "drop"
  )

if (
  tvp_pulso_disponible &&
    "beta_tvp_pulso" %in% names(comparacion)
) {

  comparacion_anual <- comparacion |>
    mutate(
      anio = year(fecha)
    ) |>
    group_by(anio) |>
    summarise(
      observaciones = n(),
      beta_rolling_36 =
        mean(beta_rolling_36, na.rm = TRUE),
      beta_lp_rolling_36 =
        mean(beta_lp_rolling_36, na.rm = TRUE),
      beta_tvp_base =
        mean(beta_tvp_base, na.rm = TRUE),
      beta_tvp_pulso =
        mean(beta_tvp_pulso, na.rm = TRUE),
      beta_rolling_30 =
        mean(beta_rolling_30, na.rm = TRUE),
      .groups = "drop"
    )
}

readr::write_csv(
  comparacion_anual,
  file.path(
    carpeta_salida,
    "06_promedios_anuales_rolling_tvp.csv"
  )
)

tabla_correlaciones_anuales <- bind_rows(
  tibble(
    serie_1 =
      "Rolling 36m: beta corto plazo",
    serie_2 =
      "TVP-AR(1) base",
    correlacion_pearson =
      correlacion_segura(
        comparacion_anual$beta_rolling_36,
        comparacion_anual$beta_tvp_base,
        "pearson"
      ),
    correlacion_spearman =
      correlacion_segura(
        comparacion_anual$beta_rolling_36,
        comparacion_anual$beta_tvp_base,
        "spearman"
      )
  ),

  tibble(
    serie_1 =
      "Rolling 36m: beta largo plazo",
    serie_2 =
      "TVP-AR(1) base",
    correlacion_pearson =
      correlacion_segura(
        comparacion_anual$beta_lp_rolling_36,
        comparacion_anual$beta_tvp_base,
        "pearson"
      ),
    correlacion_spearman =
      correlacion_segura(
        comparacion_anual$beta_lp_rolling_36,
        comparacion_anual$beta_tvp_base,
        "spearman"
      )
  )
)

readr::write_csv(
  tabla_correlaciones_anuales,
  file.path(
    carpeta_salida,
    "07_correlaciones_promedios_anuales.csv"
  )
)


# 12. RESUMEN DE TENDENCIAS ----------------------------------------------------

resumir_serie <- function(
    fecha,
    serie,
    nombre
) {

  validos <-
    !is.na(fecha) &
    is.finite(serie)

  fecha_valida <- fecha[validos]
  serie_valida <- serie[validos]

  if (length(serie_valida) == 0) {
    return(
      tibble(
        metodo = nombre,
        observaciones = 0L,
        fecha_inicial = as.Date(NA),
        valor_inicial = NA_real_,
        fecha_final = as.Date(NA),
        valor_final = NA_real_,
        cambio_absoluto = NA_real_,
        reduccion_porcentual = NA_real_,
        promedio = NA_real_,
        desviacion_estandar = NA_real_,
        minimo = NA_real_,
        fecha_minimo = as.Date(NA),
        maximo = NA_real_,
        fecha_maximo = as.Date(NA),
        pendiente_anual = NA_real_,
        direccion_tendencia = "No disponible"
      )
    )
  }

  indice_min <- which.min(serie_valida)
  indice_max <- which.max(serie_valida)

  pendiente <- pendiente_anual(
    fecha_valida,
    serie_valida
  )

  tibble(
    metodo = nombre,
    observaciones = length(serie_valida),
    fecha_inicial = first(fecha_valida),
    valor_inicial = first(serie_valida),
    fecha_final = last(fecha_valida),
    valor_final = last(serie_valida),
    cambio_absoluto =
      last(serie_valida) -
        first(serie_valida),
    reduccion_porcentual = ifelse(
      first(serie_valida) != 0,
      100 * (
        1 -
          last(serie_valida) /
            first(serie_valida)
      ),
      NA_real_
    ),
    promedio = mean(serie_valida),
    desviacion_estandar = sd(serie_valida),
    minimo = serie_valida[indice_min],
    fecha_minimo = fecha_valida[indice_min],
    maximo = serie_valida[indice_max],
    fecha_maximo = fecha_valida[indice_max],
    pendiente_anual = pendiente,
    direccion_tendencia =
      clasificar_tendencia(pendiente)
  )
}


lista_resumen <- list(
  resumir_serie(
    comparacion$fecha,
    comparacion$beta_rolling_36,
    "Rolling 36m: beta corto plazo"
  ),

  resumir_serie(
    comparacion$fecha,
    comparacion$beta_lp_rolling_36,
    "Rolling 36m: beta largo plazo"
  ),

  resumir_serie(
    comparacion$fecha,
    comparacion$rho_rolling_36,
    "Rolling 36m: rho"
  ),

  resumir_serie(
    comparacion$fecha,
    comparacion$gamma_tc_rolling_36,
    "Rolling 36m: gamma tipo de cambio"
  ),

  resumir_serie(
    comparacion$fecha,
    comparacion$beta_rolling_30,
    "Rolling 30m: beta corto plazo"
  ),

  resumir_serie(
    comparacion$fecha,
    comparacion$beta_tvp_base,
    "TVP-AR(1) base"
  )
)

if (
  tvp_pulso_disponible &&
    "beta_tvp_pulso" %in% names(comparacion)
) {

  lista_resumen <- c(
    lista_resumen,
    list(
      resumir_serie(
        comparacion$fecha,
        comparacion$beta_tvp_pulso,
        "TVP-AR(1) con pulso"
      )
    )
  )
}

tabla_resumen_tendencias <- bind_rows(
  lista_resumen
)

readr::write_csv(
  tabla_resumen_tendencias,
  file.path(
    carpeta_salida,
    "08_resumen_tendencias_comparadas.csv"
  )
)


# 13. CONCORDANCIA DE DIRECCIÓN ------------------------------------------------

tabla_concordancia <- comparacion |>
  transmute(
    fecha,

    direccion_rolling_36 =
      signo_cambio(d_beta_rolling_36),

    direccion_rolling_lp_36 =
      signo_cambio(d_beta_lp_rolling_36),

    direccion_tvp_base =
      signo_cambio(d_beta_tvp_base),

    coincide_rolling36_tvp =
      direccion_rolling_36 ==
        direccion_tvp_base,

    coincide_rolling_lp_tvp =
      direccion_rolling_lp_36 ==
        direccion_tvp_base
  )

resumen_concordancia <- tabla_concordancia |>
  summarise(
    observaciones_corto_plazo =
      sum(!is.na(coincide_rolling36_tvp)),

    porcentaje_misma_direccion_corto_plazo =
      100 * mean(
        coincide_rolling36_tvp,
        na.rm = TRUE
      ),

    observaciones_largo_plazo =
      sum(!is.na(coincide_rolling_lp_tvp)),

    porcentaje_misma_direccion_largo_plazo =
      100 * mean(
        coincide_rolling_lp_tvp,
        na.rm = TRUE
      )
  )

readr::write_csv(
  tabla_concordancia,
  file.path(
    carpeta_salida,
    "09_concordancia_direccion_mensual.csv"
  )
)

readr::write_csv(
  resumen_concordancia,
  file.path(
    carpeta_salida,
    "10_resumen_concordancia_direccion.csv"
  )
)


# 14. RESUMEN POR SUBPERÍODOS --------------------------------------------------

asignar_subperiodo <- function(fecha) {

  case_when(
    fecha <= as.Date("2016-12-01") ~ "2014-2016",
    fecha <= as.Date("2019-12-01") ~ "2017-2019",
    fecha <= as.Date("2022-12-01") ~ "2020-2022",
    TRUE ~ "2023-2026"
  )
}


tabla_subperiodos <- comparacion |>
  mutate(
    subperiodo = asignar_subperiodo(fecha)
  ) |>
  group_by(subperiodo) |>
  summarise(
    fecha_inicial = min(fecha),
    fecha_final = max(fecha),
    observaciones = n(),

    promedio_beta_rolling_36 =
      mean(beta_rolling_36, na.rm = TRUE),

    promedio_beta_lp_rolling_36 =
      mean(beta_lp_rolling_36, na.rm = TRUE),

    promedio_rho_rolling_36 =
      mean(rho_rolling_36, na.rm = TRUE),

    promedio_gamma_tc_rolling_36 =
      mean(gamma_tc_rolling_36, na.rm = TRUE),

    promedio_beta_tvp_base =
      mean(beta_tvp_base, na.rm = TRUE),

    desviacion_beta_rolling_36 =
      sd(beta_rolling_36, na.rm = TRUE),

    desviacion_beta_lp_rolling_36 =
      sd(beta_lp_rolling_36, na.rm = TRUE),

    desviacion_beta_tvp_base =
      sd(beta_tvp_base, na.rm = TRUE),

    .groups = "drop"
  )

if (
  tvp_pulso_disponible &&
    "beta_tvp_pulso" %in% names(comparacion)
) {

  tabla_subperiodos <- comparacion |>
    mutate(
      subperiodo = asignar_subperiodo(fecha)
    ) |>
    group_by(subperiodo) |>
    summarise(
      fecha_inicial = min(fecha),
      fecha_final = max(fecha),
      observaciones = n(),

      promedio_beta_rolling_36 =
        mean(beta_rolling_36, na.rm = TRUE),

      promedio_beta_lp_rolling_36 =
        mean(beta_lp_rolling_36, na.rm = TRUE),

      promedio_rho_rolling_36 =
        mean(rho_rolling_36, na.rm = TRUE),

      promedio_gamma_tc_rolling_36 =
        mean(gamma_tc_rolling_36, na.rm = TRUE),

      promedio_beta_tvp_base =
        mean(beta_tvp_base, na.rm = TRUE),

      promedio_beta_tvp_pulso =
        mean(beta_tvp_pulso, na.rm = TRUE),

      desviacion_beta_rolling_36 =
        sd(beta_rolling_36, na.rm = TRUE),

      desviacion_beta_lp_rolling_36 =
        sd(beta_lp_rolling_36, na.rm = TRUE),

      desviacion_beta_tvp_base =
        sd(beta_tvp_base, na.rm = TRUE),

      desviacion_beta_tvp_pulso =
        sd(beta_tvp_pulso, na.rm = TRUE),

      .groups = "drop"
    )
}

readr::write_csv(
  tabla_subperiodos,
  file.path(
    carpeta_salida,
    "11_resumen_por_subperiodos.csv"
  )
)


# 15. GRÁFICAS -----------------------------------------------------------------

# 15.1 Comparación principal estandarizada -------------------------------------

datos_principales <- comparacion |>
  select(
    fecha,
    `Rolling 36m: beta corto plazo` =
      z_beta_rolling_36,
    `TVP-AR(1) base` =
      z_beta_tvp_base
  ) |>
  pivot_longer(
    cols = -fecha,
    names_to = "metodo",
    values_to = "sensibilidad_estandarizada"
  )

p_principal <- ggplot(
  datos_principales,
  aes(
    x = fecha,
    y = sensibilidad_estandarizada,
    linetype = metodo
  )
) +
  geom_line(linewidth = 0.9) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  labs(
    title =
      "Sensibilidad inflacionaria: rolling revisado y TVP",
    subtitle =
      "Rolling de 36 meses y TVP-AR(1); series estandarizadas en la muestra común",
    x = NULL,
    y = "Sensibilidad estandarizada",
    linetype = NULL,
    caption = paste(
      "Las magnitudes originales no son equivalentes.",
      "La comparación se refiere a la trayectoria temporal."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    carpeta_salida,
    "12_comparacion_beta_corto_plazo_rolling36_tvp.png"
  ),
  plot = p_principal,
  width = 10,
  height = 6,
  dpi = 300
)


# 15.2 Sensibilidad rolling de largo plazo frente a TVP ------------------------

datos_lp <- comparacion |>
  select(
    fecha,
    `Rolling 36m: beta largo plazo` =
      z_beta_lp_rolling_36,
    `TVP-AR(1) base` =
      z_beta_tvp_base
  ) |>
  pivot_longer(
    cols = -fecha,
    names_to = "metodo",
    values_to = "sensibilidad_estandarizada"
  )

p_lp <- ggplot(
  datos_lp,
  aes(
    x = fecha,
    y = sensibilidad_estandarizada,
    linetype = metodo
  )
) +
  geom_line(linewidth = 0.9) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  labs(
    title =
      "Sensibilidad de largo plazo rolling y trayectoria TVP",
    subtitle =
      "Comparación complementaria mediante series estandarizadas",
    x = NULL,
    y = "Sensibilidad estandarizada",
    linetype = NULL,
    caption =
      "La sensibilidad rolling de largo plazo se calcula como beta/(1-rho)."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    carpeta_salida,
    "13_comparacion_beta_largo_plazo_rolling36_tvp.png"
  ),
  plot = p_lp,
  width = 10,
  height = 6,
  dpi = 300
)


# 15.3 Índices relativos de anclaje --------------------------------------------

datos_indices <- comparacion |>
  select(
    fecha,
    `Rolling 36m: corto plazo` =
      indice_anclaje_rolling_36,
    `Rolling 36m: largo plazo` =
      indice_anclaje_lp_rolling_36,
    `TVP-AR(1)` =
      indice_anclaje_tvp_base
  ) |>
  pivot_longer(
    cols = -fecha,
    names_to = "metodo",
    values_to = "indice"
  )

p_indices <- ggplot(
  datos_indices,
  aes(
    x = fecha,
    y = indice,
    linetype = metodo
  )
) +
  geom_line(linewidth = 0.85) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title =
      "Índices relativos de anclaje por método",
    subtitle =
      "Mayor valor indica menor sensibilidad dentro de la muestra común",
    x = NULL,
    y = "Índice relativo de anclaje",
    linetype = NULL,
    caption = paste(
      "Índice min-max invertido.",
      "No constituye una medida estructural en unidades económicas."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    carpeta_salida,
    "14_indices_relativos_anclaje.png"
  ),
  plot = p_indices,
  width = 10,
  height = 6,
  dpi = 300
)


# 15.4 Robustez de ventanas rolling --------------------------------------------

datos_ventanas <- bind_rows(
  rolling_36 |>
    transmute(
      fecha = fecha_fin,
      ventana = "36 meses",
      beta = beta_corto_plazo
    ),

  rolling_30 |>
    transmute(
      fecha = fecha_fin,
      ventana = "30 meses",
      beta = beta_corto_plazo
    )
)

p_ventanas <- ggplot(
  datos_ventanas,
  aes(
    x = fecha,
    y = beta,
    linetype = ventana
  )
) +
  geom_line(linewidth = 0.85) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title =
      "Robustez del coeficiente rolling a la longitud de ventana",
    subtitle =
      "Especificación revisada con ventanas de 36 y 30 meses",
    x = NULL,
    y = expression(beta[w]),
    linetype = "Ventana"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    carpeta_salida,
    "15_rolling_36_vs_30_meses.png"
  ),
  plot = p_ventanas,
  width = 10,
  height = 6,
  dpi = 300
)


# 15.5 Trayectoria de rho -------------------------------------------------------

p_rho <- ggplot(
  rolling_36,
  aes(
    x = fecha_fin,
    y = rho
  )
) +
  geom_line(linewidth = 0.9) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    title =
      "Persistencia de la brecha de expectativas",
    subtitle =
      "Coeficiente rho en ventanas móviles de 36 meses",
    x = NULL,
    y = expression(rho[w])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    carpeta_salida,
    "16_trayectoria_rho_rolling_36m.png"
  ),
  plot = p_rho,
  width = 10,
  height = 5.8,
  dpi = 300
)


# 15.6 Trayectoria de gamma cambiario ------------------------------------------

p_gamma <- ggplot(
  rolling_36,
  aes(
    x = fecha_fin,
    y = gamma_tc
  )
) +
  geom_line(linewidth = 0.9) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title =
      "Sensibilidad rolling al tipo de cambio",
    subtitle =
      "Coeficiente de la depreciación cambiaria rezagada; ventana de 36 meses",
    x = NULL,
    y = expression(gamma[w])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    carpeta_salida,
    "17_trayectoria_gamma_tc_rolling_36m.png"
  ),
  plot = p_gamma,
  width = 10,
  height = 5.8,
  dpi = 300
)


# 15.7 Promedios anuales estandarizados ----------------------------------------

datos_anuales <- comparacion_anual |>
  mutate(
    z_rolling_36 =
      estandarizar(beta_rolling_36),

    z_tvp_base =
      estandarizar(beta_tvp_base)
  ) |>
  select(
    anio,
    `Rolling 36m` = z_rolling_36,
    `TVP-AR(1)` = z_tvp_base
  ) |>
  pivot_longer(
    cols = -anio,
    names_to = "metodo",
    values_to = "valor"
  )

p_anual <- ggplot(
  datos_anuales,
  aes(
    x = anio,
    y = valor,
    linetype = metodo
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point() +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    title =
      "Promedios anuales de la sensibilidad estimada",
    subtitle =
      "Series estandarizadas en la muestra común",
    x = NULL,
    y = "Promedio anual estandarizado",
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    carpeta_salida,
    "18_comparacion_promedios_anuales.png"
  ),
  plot = p_anual,
  width = 10,
  height = 6,
  dpi = 300
)


# 15.8 TVP base y robusto frente a rolling -------------------------------------

if (
  tvp_pulso_disponible &&
    "z_beta_tvp_pulso" %in% names(comparacion)
) {

  datos_robustez <- comparacion |>
    select(
      fecha,
      `Rolling 36m` =
        z_beta_rolling_36,
      `TVP base` =
        z_beta_tvp_base,
      `TVP con pulso` =
        z_beta_tvp_pulso
    ) |>
    pivot_longer(
      cols = -fecha,
      names_to = "metodo",
      values_to = "valor"
    )

  p_robustez <- ggplot(
    datos_robustez,
    aes(
      x = fecha,
      y = valor,
      linetype = metodo
    )
  ) +
    geom_line(linewidth = 0.85) +
    geom_hline(
      yintercept = 0,
      linetype = "dotted"
    ) +
    labs(
      title =
        "Rolling revisado y robustez de la trayectoria TVP",
      subtitle =
        "Series estandarizadas en la muestra común",
      x = NULL,
      y = "Sensibilidad estandarizada",
      linetype = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title.position = "plot",
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  ggsave(
    file.path(
      carpeta_salida,
      "19_rolling36_tvp_base_y_pulso.png"
    ),
    plot = p_robustez,
    width = 10,
    height = 6,
    dpi = 300
  )
}


# 16. RESUMEN EJECUTIVO --------------------------------------------------------

fila_principal <- tabla_correlaciones |>
  filter(
    serie_1 ==
      "Rolling 36m: beta corto plazo",
    serie_2 ==
      "TVP-AR(1) base"
  )

fila_largo_plazo <- tabla_correlaciones |>
  filter(
    serie_1 ==
      "Rolling 36m: beta largo plazo",
    serie_2 ==
      "TVP-AR(1) base"
  )

fila_ventanas <- tabla_correlaciones |>
  filter(
    serie_1 == "Rolling 36m",
    serie_2 == "Rolling 30m"
  )

resumen_roll <- tabla_resumen_tendencias |>
  filter(
    metodo ==
      "Rolling 36m: beta corto plazo"
  )

resumen_roll_lp <- tabla_resumen_tendencias |>
  filter(
    metodo ==
      "Rolling 36m: beta largo plazo"
  )

resumen_tvp <- tabla_resumen_tendencias |>
  filter(
    metodo == "TVP-AR(1) base"
  )


formato_num <- function(x, digitos = 4) {

  if (
    length(x) == 0 ||
      !is.finite(x[[1]])
  ) {
    return("NA")
  }

  formatC(
    x[[1]],
    digits = digitos,
    format = "f"
  )
}


lineas_resumen <- c(
  "",
  "======================================================================",
  "COMPARACIÓN: ROLLING REVISADO VS. TVP-AR(1)",
  "======================================================================",
  paste0("Base: ", archivo_bd),
  paste0("Variable cambiaria utilizada: ", col_tc),
  paste0("TVP base: ", archivo_tvp_base),
  paste0(
    "TVP con pulso disponible: ",
    tvp_pulso_disponible
  ),
  paste0(
    "Muestra común: ",
    format(min(comparacion$fecha), "%Y-%m"),
    " a ",
    format(max(comparacion$fecha), "%Y-%m")
  ),
  paste0(
    "Observaciones comunes: ",
    nrow(comparacion)
  ),
  "",
  "ESPECIFICACIÓN ROLLING",
  paste0(
    "Ventana principal/robustez: ",
    ventana_principal,
    " / ",
    ventana_robustez,
    " meses"
  ),
  paste0(
    "Inferencia: Newey-West HAC, lag = ",
    hac_lag
  ),
  paste(
    "Controles: brecha de expectativa rezagada",
    "y depreciación cambiaria rezagada."
  ),
  "",
  "COMPARACIÓN PRINCIPAL: BETA ROLLING 36M VS. TVP BASE",
  paste0(
    "Pearson en niveles: ",
    formato_num(
      fila_principal$
        correlacion_pearson_nivel
    )
  ),
  paste0(
    "Spearman en niveles: ",
    formato_num(
      fila_principal$
        correlacion_spearman_nivel
    )
  ),
  paste0(
    "Pearson en primeras diferencias: ",
    formato_num(
      fila_principal$
        correlacion_pearson_primeras_diferencias
    )
  ),
  paste0(
    "Spearman en primeras diferencias: ",
    formato_num(
      fila_principal$
        correlacion_spearman_primeras_diferencias
    )
  ),
  paste0(
    "Coincidencia mensual de dirección: ",
    formato_num(
      resumen_concordancia$
        porcentaje_misma_direccion_corto_plazo,
      2
    ),
    "%"
  ),
  "",
  "COMPARACIÓN COMPLEMENTARIA: EFECTO ROLLING DE LARGO PLAZO",
  paste0(
    "Pearson en niveles: ",
    formato_num(
      fila_largo_plazo$
        correlacion_pearson_nivel
    )
  ),
  paste0(
    "Spearman en niveles: ",
    formato_num(
      fila_largo_plazo$
        correlacion_spearman_nivel
    )
  ),
  "",
  "TENDENCIAS EN LA MUESTRA COMÚN",
  paste0(
    "Rolling corto plazo: ",
    resumen_roll$direccion_tendencia,
    " | reducción = ",
    formato_num(
      resumen_roll$reduccion_porcentual,
      2
    ),
    "%"
  ),
  paste0(
    "Rolling largo plazo: ",
    resumen_roll_lp$direccion_tendencia,
    " | reducción = ",
    formato_num(
      resumen_roll_lp$reduccion_porcentual,
      2
    ),
    "%"
  ),
  paste0(
    "TVP base: ",
    resumen_tvp$direccion_tendencia,
    " | reducción = ",
    formato_num(
      resumen_tvp$reduccion_porcentual,
      2
    ),
    "%"
  ),
  "",
  "ROBUSTEZ DE LA VENTANA ROLLING",
  paste0(
    "Correlación rolling 36m vs. 30m: ",
    formato_num(
      fila_ventanas$
        correlacion_pearson_nivel
    )
  ),
  "",
  "NOTA METODOLÓGICA",
  paste(
    "Los coeficientes rolling y TVP no son equivalentes en unidades.",
    "La evidencia comparada se basa en trayectoria, tendencia y comovimiento."
  ),
  "",
  paste0(
    "Resultados guardados en: ",
    carpeta_salida
  ),
  "======================================================================"
)

writeLines(
  lineas_resumen,
  con = file.path(
    carpeta_salida,
    "00_resumen_ejecutivo.txt"
  )
)

cat(
  paste(
    lineas_resumen,
    collapse = "\n"
  ),
  "\n"
)
