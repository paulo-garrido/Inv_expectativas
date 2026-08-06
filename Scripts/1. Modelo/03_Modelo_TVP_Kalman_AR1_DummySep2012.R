# ==============================================================================
# OBJETIVO 2 - TERCERA ESPECIFICACIÓN / PRUEBA DE ROBUSTEZ
# Modelo TVP-Kalman AR(1) con intervención de pulso en septiembre de 2012
#
# MODELO BASE:
#   Delta exp_inf_24m_t = alpha
#                         + rho * Delta exp_inf_24m_(t-1)
#                         + beta_t * shock_inf_t
#                         + epsilon_t
#
# MODELO CON INTERVENCIÓN:
#   Delta exp_inf_24m_t = alpha
#                         + rho * Delta exp_inf_24m_(t-1)
#                         + gamma * D_sep2012,t
#                         + beta_t * shock_inf_t
#                         + epsilon_t
#
# Transición:
#   beta_t = beta_(t-1) + eta_beta,t
#
# La dummy D_sep2012,t toma valor 1 únicamente en septiembre de 2012.
# No representa un cambio metodológico; controla una revisión extraordinaria.
#
# Objetivos:
#   1. Verificar si la trayectoria descendente de beta_t es robusta.
#   2. Evaluar cuánto mejora el ajuste al controlar septiembre de 2012.
#   3. Revisar si mejoran normalidad y autocorrelación residual.
#
# Archivos de entrada admitidos:
#   bd_sample_modelo.csv
#   bd_sample_modelo(1).csv
# ==============================================================================


# 0. PAQUETES ------------------------------------------------------------------

paquetes <- c(
  "dplyr",
  "readr",
  "lubridate",
  "ggplot2",
  "tidyr",
  "tibble",
  "broom",
  "KFAS"
)

paquetes_faltantes <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(paquetes_faltantes) > 0) {
  install.packages(paquetes_faltantes, dependencies = TRUE)
}

invisible(lapply(paquetes, library, character.only = TRUE))


# 1. CONFIGURACIÓN -------------------------------------------------------------

archivos_candidatos <- c(
  "bd_sample_modelo.csv",
  "bd_sample_modelo(1).csv"
)

archivo_entrada <- archivos_candidatos[
  file.exists(archivos_candidatos)
][1]

if (is.na(archivo_entrada)) {
  stop(
    paste0(
      "No se encontró ninguno de estos archivos:\n",
      paste0("- ", archivos_candidatos, collapse = "\n"), "\n\n",
      "Directorio actual: ", getwd()
    )
  )
}

carpeta_salida <- file.path(
  getwd(),
  "resultados_objetivo2_tvp_ar1_dummy_sep2012"
)

fecha_inicio <- as.Date("2011-01-01")
fecha_intervencion <- as.Date("2012-09-01")

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


# 2. IMPORTACIÓN Y VALIDACIÓN --------------------------------------------------

bd_original <- readr::read_csv(
  archivo_entrada,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", ".", "null")
)

variables_requeridas <- c(
  "fecha",
  "exp_inf_12m",
  "exp_inf_24m",
  "infl_gt"
)

variables_faltantes <- setdiff(
  variables_requeridas,
  names(bd_original)
)

if (length(variables_faltantes) > 0) {
  stop(
    paste(
      "Faltan las siguientes variables requeridas:",
      paste(variables_faltantes, collapse = ", ")
    )
  )
}


# 3. CONSTRUCCIÓN DE VARIABLES -------------------------------------------------

bd_modelo <- bd_original %>%
  mutate(
    fecha = lubridate::ymd(fecha)
  ) %>%
  arrange(fecha) %>%
  mutate(
    d_exp24 = exp_inf_24m - lag(exp_inf_24m),
    d_exp24_lag1 = lag(d_exp24),
    exp12_lag1 = lag(exp_inf_12m),
    shock_inf = infl_gt - exp12_lag1,

    # Dummy de pulso: 1 solo en septiembre de 2012.
    dummy_sep2012 = as.numeric(fecha == fecha_intervencion)
  ) %>%
  filter(fecha >= fecha_inicio) %>%
  filter(
    is.finite(d_exp24),
    is.finite(d_exp24_lag1),
    is.finite(shock_inf)
  )

if (nrow(bd_modelo) < 36) {
  stop("La muestra efectiva tiene menos de 36 observaciones.")
}

if (sum(bd_modelo$dummy_sep2012) != 1) {
  stop(
    "La dummy de septiembre de 2012 no contiene exactamente una observación. ",
    "Revise las fechas de la base."
  )
}


# 4. FUNCIONES AUXILIARES ------------------------------------------------------

# 4.1 Construcción del modelo de espacio de estados ----------------------------

construir_modelo_tvp <- function(
    y,
    y_lag1,
    x,
    valores_iniciales_estado,
    dummy = NULL
) {

  n <- length(y)
  usar_dummy <- !is.null(dummy)

  if (usar_dummy) {

    # Estados:
    # 1 alpha constante
    # 2 rho constante
    # 3 gamma constante (pulso septiembre 2012)
    # 4 beta_t variante
    m <- 4

    Z_t <- array(
      NA_real_,
      dim = c(1, m, n)
    )

    Z_t[1, 1, ] <- 1
    Z_t[1, 2, ] <- y_lag1
    Z_t[1, 3, ] <- dummy
    Z_t[1, 4, ] <- x

    R_beta <- matrix(
      c(0, 0, 0, 1),
      nrow = m,
      ncol = 1
    )

  } else {

    # Estados:
    # 1 alpha constante
    # 2 rho constante
    # 3 beta_t variante
    m <- 3

    Z_t <- array(
      NA_real_,
      dim = c(1, m, n)
    )

    Z_t[1, 1, ] <- 1
    Z_t[1, 2, ] <- y_lag1
    Z_t[1, 3, ] <- x

    R_beta <- matrix(
      c(0, 0, 1),
      nrow = m,
      ncol = 1
    )
  }

  modelo <- KFAS::SSModel(
    y ~ -1 +
      SSMcustom(
        Z = Z_t,
        T = diag(m),
        R = R_beta,
        Q = matrix(NA_real_, nrow = 1, ncol = 1),
        n = n,
        a1 = valores_iniciales_estado,
        P1 = matrix(0, nrow = m, ncol = m),
        P1inf = diag(m)
      ),
    H = matrix(NA_real_, nrow = 1, ncol = 1)
  )

  modelo
}


# 4.2 Actualización de hiperparámetros -----------------------------------------

actualizar_varianzas <- function(par, modelo) {

  # par[1] = log(H)
  # par[2] = log(Q_beta)
  modelo$H[1, 1, 1] <- exp(par[1])
  modelo$Q[1, 1, 1] <- exp(par[2])

  modelo
}


# 4.3 Estimación con múltiples valores iniciales -------------------------------

estimar_modelo <- function(
    modelo_sin_estimar,
    var_residuo_ols,
    var_y,
    etiqueta
) {

  lista_inicios <- list(
    log(c(var_residuo_ols, var_y * 1e-2)),
    log(c(var_residuo_ols, var_y * 1e-3)),
    log(c(var_residuo_ols, var_y * 1e-4)),
    log(c(var_residuo_ols * 0.5, var_y * 1e-3)),
    log(c(var_residuo_ols * 2.0, var_y * 1e-3)),
    log(c(var_residuo_ols, 1e-6))
  )

  ajustes <- lapply(
    seq_along(lista_inicios),
    function(i) {

      intento <- tryCatch(
        KFAS::fitSSM(
          model = modelo_sin_estimar,
          inits = lista_inicios[[i]],
          updatefn = actualizar_varianzas,
          method = "L-BFGS-B",
          lower = rep(log(1e-10), 2),
          upper = rep(log(100), 2),
          control = list(
            maxit = 5000,
            factr = 1e7
          )
        ),
        error = function(e) NULL
      )

      if (is.null(intento)) {
        return(NULL)
      }

      log_verosimilitud <- tryCatch(
        as.numeric(logLik(intento$model)),
        error = function(e) NA_real_
      )

      list(
        ajuste = intento,
        inicio = i,
        convergencia = intento$optim.out$convergence,
        mensaje = intento$optim.out$message,
        logLik = log_verosimilitud
      )
    }
  )

  ajustes <- Filter(Negate(is.null), ajustes)

  if (length(ajustes) == 0) {
    stop("No fue posible estimar el modelo: ", etiqueta)
  }

  tabla_intentos <- tibble::tibble(
    modelo = etiqueta,
    intento = vapply(ajustes, function(z) z$inicio, numeric(1)),
    convergencia = vapply(
      ajustes,
      function(z) z$convergencia,
      numeric(1)
    ),
    logLik = vapply(ajustes, function(z) z$logLik, numeric(1)),
    mensaje = vapply(
      ajustes,
      function(z) {
        if (is.null(z$mensaje)) "" else as.character(z$mensaje)
      },
      character(1)
    )
  ) %>%
    arrange(desc(logLik))

  ajustes_convergentes <- ajustes[
    vapply(
      ajustes,
      function(z) z$convergencia == 0 && is.finite(z$logLik),
      logical(1)
    )
  ]

  if (length(ajustes_convergentes) == 0) {
    warning(
      "Ningún intento convergió normalmente para: ",
      etiqueta,
      ". Se utilizará la mayor log-verosimilitud disponible."
    )
    elegibles <- ajustes
  } else {
    elegibles <- ajustes_convergentes
  }

  indice_mejor <- which.max(
    vapply(elegibles, function(z) z$logLik, numeric(1))
  )

  mejor <- elegibles[[indice_mejor]]

  list(
    mejor = mejor,
    modelo = mejor$ajuste$model,
    tabla_intentos = tabla_intentos
  )
}


# 4.4 Extracción de resultados -------------------------------------------------

extraer_resultados <- function(
    modelo_estimado,
    bd,
    usar_dummy = FALSE
) {

  kfs <- KFAS::KFS(
    modelo_estimado,
    filtering = c("state", "mean"),
    smoothing = c("state", "mean", "disturbance")
  )

  if (usar_dummy) {

    alpha_t <- as.numeric(kfs$alphahat[, 1])
    rho_t   <- as.numeric(kfs$alphahat[, 2])
    gamma_t <- as.numeric(kfs$alphahat[, 3])
    beta_t  <- as.numeric(kfs$alphahat[, 4])

    se_alpha_t <- sqrt(pmax(as.numeric(kfs$V[1, 1, ]), 0))
    se_rho_t   <- sqrt(pmax(as.numeric(kfs$V[2, 2, ]), 0))
    se_gamma_t <- sqrt(pmax(as.numeric(kfs$V[3, 3, ]), 0))
    se_beta_t  <- sqrt(pmax(as.numeric(kfs$V[4, 4, ]), 0))

    ajustado <- alpha_t +
      rho_t * bd$d_exp24_lag1 +
      gamma_t * bd$dummy_sep2012 +
      beta_t * bd$shock_inf

  } else {

    alpha_t <- as.numeric(kfs$alphahat[, 1])
    rho_t   <- as.numeric(kfs$alphahat[, 2])
    gamma_t <- rep(NA_real_, nrow(bd))
    beta_t  <- as.numeric(kfs$alphahat[, 3])

    se_alpha_t <- sqrt(pmax(as.numeric(kfs$V[1, 1, ]), 0))
    se_rho_t   <- sqrt(pmax(as.numeric(kfs$V[2, 2, ]), 0))
    se_gamma_t <- rep(NA_real_, nrow(bd))
    se_beta_t  <- sqrt(pmax(as.numeric(kfs$V[3, 3, ]), 0))

    ajustado <- alpha_t +
      rho_t * bd$d_exp24_lag1 +
      beta_t * bd$shock_inf
  }

  residuos_estandarizados <- tryCatch(
    as.numeric(stats::rstandard(kfs, type = "recursive")),
    error = function(e) rep(NA_real_, nrow(bd))
  )

  tabla <- bd %>%
    transmute(
      fecha,
      exp_inf_12m,
      exp_inf_24m,
      infl_gt,
      d_exp24,
      d_exp24_lag1,
      exp12_lag1,
      shock_inf,
      dummy_sep2012,

      alpha = alpha_t,
      se_alpha = se_alpha_t,

      rho = rho_t,
      se_rho = se_rho_t,

      gamma_sep2012 = gamma_t,
      se_gamma_sep2012 = se_gamma_t,

      beta_t = beta_t,
      se_beta_t = se_beta_t,
      beta_li_95 = beta_t - z_critico * se_beta_t,
      beta_ls_95 = beta_t + z_critico * se_beta_t,

      indice_anclaje = -beta_t,

      ajustado_tvp = ajustado,
      residuo_tvp = d_exp24 - ajustado,
      residuo_estandarizado = residuos_estandarizados
    )

  list(
    kfs = kfs,
    tabla = tabla
  )
}


# 4.5 Diagnósticos -------------------------------------------------------------

calcular_diagnosticos <- function(
    resultados,
    nombre_modelo
) {

  residuos <- resultados$residuo_estandarizado
  residuos <- residuos[is.finite(residuos)]

  if (length(residuos) < 30) {
    return(
      tibble::tibble(
        modelo = nombre_modelo,
        diagnostico = "No disponible",
        rezago = NA_real_,
        estadistico = NA_real_,
        grados_libertad = NA_real_,
        p_valor = NA_real_,
        conclusion_5pct = "No disponible"
      )
    )
  }

  # fitdf = 1 descuenta el componente AR(1).
  rezagos_lb <- c(6, 12, 18, 24)
  rezagos_lb <- rezagos_lb[rezagos_lb < length(residuos)]

  tabla_lb <- lapply(
    rezagos_lb,
    function(L) {

      prueba <- Box.test(
        residuos,
        lag = L,
        type = "Ljung-Box",
        fitdf = 1
      )

      tibble::tibble(
        modelo = nombre_modelo,
        diagnostico = "Ljung-Box",
        rezago = L,
        estadistico = unname(prueba$statistic),
        grados_libertad = unname(prueba$parameter),
        p_valor = prueba$p.value,
        conclusion_5pct = ifelse(
          prueba$p.value < 0.05,
          "Se rechaza ausencia de autocorrelación",
          "No se rechaza ausencia de autocorrelación"
        )
      )
    }
  ) %>%
    dplyr::bind_rows()

  media_r <- mean(residuos)
  desv_r  <- sd(residuos)

  if (!is.finite(desv_r) || desv_r <= 0) {

    jb <- tibble::tibble(
      modelo = nombre_modelo,
      diagnostico = "Jarque-Bera aproximado",
      rezago = NA_real_,
      estadistico = NA_real_,
      grados_libertad = 2,
      p_valor = NA_real_,
      conclusion_5pct = "No disponible"
    )

  } else {

    asimetria <- mean(((residuos - media_r) / desv_r)^3)
    curtosis  <- mean(((residuos - media_r) / desv_r)^4)

    estadistico_jb <- length(residuos) / 6 *
      (
        asimetria^2 +
          ((curtosis - 3)^2) / 4
      )

    p_jb <- pchisq(
      estadistico_jb,
      df = 2,
      lower.tail = FALSE
    )

    jb <- tibble::tibble(
      modelo = nombre_modelo,
      diagnostico = "Jarque-Bera aproximado",
      rezago = NA_real_,
      estadistico = estadistico_jb,
      grados_libertad = 2,
      p_valor = p_jb,
      conclusion_5pct = ifelse(
        p_jb < 0.05,
        "Se rechaza normalidad",
        "No se rechaza normalidad"
      )
    )
  }

  bind_rows(tabla_lb, jb)
}


# 5. PREPARACIÓN DE MODELOS ----------------------------------------------------

y      <- bd_modelo$d_exp24
y_lag1 <- bd_modelo$d_exp24_lag1
x      <- bd_modelo$shock_inf
d      <- bd_modelo$dummy_sep2012

# Modelos OLS usados únicamente para inicialización.
ols_base <- lm(
  d_exp24 ~ d_exp24_lag1 + shock_inf,
  data = bd_modelo
)

ols_dummy <- lm(
  d_exp24 ~ d_exp24_lag1 + dummy_sep2012 + shock_inf,
  data = bd_modelo
)

coef_base <- coef(ols_base)
coef_dummy <- coef(ols_dummy)

if (any(!is.finite(coef_base)) || any(!is.finite(coef_dummy))) {
  stop("Los modelos OLS de inicialización produjeron coeficientes no finitos.")
}

estado_inicial_base <- c(
  alpha = unname(coef_base["(Intercept)"]),
  rho   = unname(coef_base["d_exp24_lag1"]),
  beta  = unname(coef_base["shock_inf"])
)

estado_inicial_dummy <- c(
  alpha = unname(coef_dummy["(Intercept)"]),
  rho   = unname(coef_dummy["d_exp24_lag1"]),
  gamma = unname(coef_dummy["dummy_sep2012"]),
  beta  = unname(coef_dummy["shock_inf"])
)

modelo_base_sin_estimar <- construir_modelo_tvp(
  y = y,
  y_lag1 = y_lag1,
  x = x,
  valores_iniciales_estado = estado_inicial_base,
  dummy = NULL
)

modelo_dummy_sin_estimar <- construir_modelo_tvp(
  y = y,
  y_lag1 = y_lag1,
  x = x,
  valores_iniciales_estado = estado_inicial_dummy,
  dummy = d
)

var_residuo_base <- max(var(residuals(ols_base), na.rm = TRUE), 1e-8)
var_residuo_dummy <- max(var(residuals(ols_dummy), na.rm = TRUE), 1e-8)
var_y <- max(var(y, na.rm = TRUE), 1e-8)


# 6. ESTIMACIÓN ---------------------------------------------------------------

ajuste_base <- estimar_modelo(
  modelo_sin_estimar = modelo_base_sin_estimar,
  var_residuo_ols = var_residuo_base,
  var_y = var_y,
  etiqueta = "TVP-AR(1) base"
)

ajuste_dummy <- estimar_modelo(
  modelo_sin_estimar = modelo_dummy_sin_estimar,
  var_residuo_ols = var_residuo_dummy,
  var_y = var_y,
  etiqueta = "TVP-AR(1) + pulso septiembre 2012"
)

readr::write_csv(
  ajuste_base$tabla_intentos,
  file.path(carpeta_salida, "01_intentos_modelo_base.csv")
)

readr::write_csv(
  ajuste_dummy$tabla_intentos,
  file.path(carpeta_salida, "02_intentos_modelo_dummy.csv")
)


# 7. FILTRO, SUAVIZADOR Y RESULTADOS ------------------------------------------

salida_base <- extraer_resultados(
  modelo_estimado = ajuste_base$modelo,
  bd = bd_modelo,
  usar_dummy = FALSE
)

salida_dummy <- extraer_resultados(
  modelo_estimado = ajuste_dummy$modelo,
  bd = bd_modelo,
  usar_dummy = TRUE
)

resultados_base <- salida_base$tabla
resultados_dummy <- salida_dummy$tabla

readr::write_csv(
  resultados_base,
  file.path(carpeta_salida, "03_resultados_modelo_base.csv")
)

readr::write_csv(
  resultados_dummy,
  file.path(carpeta_salida, "04_resultados_modelo_dummy.csv")
)


# 8. COEFICIENTES FIJOS DEL MODELO CON DUMMY ----------------------------------

alpha_dummy <- mean(resultados_dummy$alpha, na.rm = TRUE)
rho_dummy   <- mean(resultados_dummy$rho, na.rm = TRUE)
gamma_dummy <- mean(resultados_dummy$gamma_sep2012, na.rm = TRUE)

se_alpha_dummy <- mean(resultados_dummy$se_alpha, na.rm = TRUE)
se_rho_dummy   <- mean(resultados_dummy$se_rho, na.rm = TRUE)
se_gamma_dummy <- mean(resultados_dummy$se_gamma_sep2012, na.rm = TRUE)

tabla_coeficientes_dummy <- tibble::tibble(
  parametro = c("alpha", "rho", "gamma_sep2012"),
  estimacion = c(alpha_dummy, rho_dummy, gamma_dummy),
  error_estandar = c(
    se_alpha_dummy,
    se_rho_dummy,
    se_gamma_dummy
  ),
  limite_inferior_95 = c(
    alpha_dummy - z_critico * se_alpha_dummy,
    rho_dummy - z_critico * se_rho_dummy,
    gamma_dummy - z_critico * se_gamma_dummy
  ),
  limite_superior_95 = c(
    alpha_dummy + z_critico * se_alpha_dummy,
    rho_dummy + z_critico * se_rho_dummy,
    gamma_dummy + z_critico * se_gamma_dummy
  )
) %>%
  mutate(
    estadistico_z = estimacion / error_estandar,
    p_valor = 2 * pnorm(abs(estadistico_z), lower.tail = FALSE),
    significativo_5pct = p_valor < 0.05
  )

readr::write_csv(
  tabla_coeficientes_dummy,
  file.path(carpeta_salida, "05_coeficientes_fijos_modelo_dummy.csv")
)


# 9. VARIANZAS Y COMPARACIÓN DE AJUSTE ----------------------------------------

H_base <- as.numeric(ajuste_base$modelo$H[1, 1, 1])
Q_base <- as.numeric(ajuste_base$modelo$Q[1, 1, 1])

H_dummy <- as.numeric(ajuste_dummy$modelo$H[1, 1, 1])
Q_dummy <- as.numeric(ajuste_dummy$modelo$Q[1, 1, 1])

loglik_base <- as.numeric(logLik(ajuste_base$modelo))
loglik_dummy <- as.numeric(logLik(ajuste_dummy$modelo))

# Conteo homogéneo aproximado:
# Base: alpha, rho, beta inicial, H, Q_beta = 5
# Dummy: alpha, rho, gamma, beta inicial, H, Q_beta = 6
k_base <- 5
k_dummy <- 6
n_obs <- nrow(bd_modelo)

LR_dummy <- 2 * (loglik_dummy - loglik_base)
p_LR_dummy <- pchisq(
  LR_dummy,
  df = 1,
  lower.tail = FALSE
)

tabla_comparacion_modelos <- tibble::tibble(
  modelo = c(
    "TVP-AR(1) base",
    "TVP-AR(1) + pulso septiembre 2012"
  ),
  log_verosimilitud = c(loglik_base, loglik_dummy),
  parametros_aproximados = c(k_base, k_dummy),
  AIC_aproximado = c(
    -2 * loglik_base + 2 * k_base,
    -2 * loglik_dummy + 2 * k_dummy
  ),
  BIC_aproximado = c(
    -2 * loglik_base + log(n_obs) * k_base,
    -2 * loglik_dummy + log(n_obs) * k_dummy
  ),
  H = c(H_base, H_dummy),
  Q_beta = c(Q_base, Q_dummy),
  Q_beta_sobre_H = c(Q_base / H_base, Q_dummy / H_dummy),
  beta_inicial = c(
    first(resultados_base$beta_t),
    first(resultados_dummy$beta_t)
  ),
  beta_final = c(
    last(resultados_base$beta_t),
    last(resultados_dummy$beta_t)
  )
)

readr::write_csv(
  tabla_comparacion_modelos,
  file.path(carpeta_salida, "06_comparacion_modelos_kfas.csv")
)

tabla_prueba_intervencion <- tibble::tibble(
  contraste = "Pulso septiembre 2012 frente al modelo base",
  estadistico_LR = LR_dummy,
  grados_libertad = 1,
  p_valor = p_LR_dummy,
  conclusion_5pct = ifelse(
    p_LR_dummy < 0.05,
    "La intervención mejora significativamente el ajuste",
    "No se identifica una mejora significativa del ajuste"
  )
)

readr::write_csv(
  tabla_prueba_intervencion,
  file.path(carpeta_salida, "07_prueba_intervencion_sep2012.csv")
)


# 10. COMPARACIÓN DE TRAYECTORIAS DE BETA -------------------------------------

comparacion_beta <- resultados_base %>%
  select(
    fecha,
    beta_base = beta_t,
    beta_base_li_95 = beta_li_95,
    beta_base_ls_95 = beta_ls_95
  ) %>%
  left_join(
    resultados_dummy %>%
      select(
        fecha,
        beta_dummy = beta_t,
        beta_dummy_li_95 = beta_li_95,
        beta_dummy_ls_95 = beta_ls_95
      ),
    by = "fecha"
  ) %>%
  mutate(
    diferencia_dummy_menos_base = beta_dummy - beta_base,
    diferencia_absoluta = abs(diferencia_dummy_menos_base)
  )

resumen_comparacion_beta <- comparacion_beta %>%
  summarise(
    correlacion_trayectorias = cor(
      beta_base,
      beta_dummy,
      use = "complete.obs"
    ),
    diferencia_media = mean(
      diferencia_dummy_menos_base,
      na.rm = TRUE
    ),
    diferencia_absoluta_media = mean(
      diferencia_absoluta,
      na.rm = TRUE
    ),
    diferencia_absoluta_maxima = max(
      diferencia_absoluta,
      na.rm = TRUE
    ),
    fecha_diferencia_maxima = fecha[
      which.max(diferencia_absoluta)
    ],
    RMSE_trayectorias = sqrt(
      mean(
        (beta_dummy - beta_base)^2,
        na.rm = TRUE
      )
    )
  )

readr::write_csv(
  comparacion_beta,
  file.path(carpeta_salida, "08_comparacion_mensual_beta.csv")
)

readr::write_csv(
  resumen_comparacion_beta,
  file.path(carpeta_salida, "09_resumen_comparacion_beta.csv")
)


# 11. DIAGNÓSTICOS COMPARADOS --------------------------------------------------

diagnosticos_base <- calcular_diagnosticos(
  resultados_base,
  "TVP-AR(1) base"
)

diagnosticos_dummy <- calcular_diagnosticos(
  resultados_dummy,
  "TVP-AR(1) + pulso septiembre 2012"
)

tabla_diagnosticos <- bind_rows(
  diagnosticos_base,
  diagnosticos_dummy
)

readr::write_csv(
  tabla_diagnosticos,
  file.path(carpeta_salida, "10_diagnosticos_comparados.csv")
)


# 12. OBSERVACIONES EXTREMAS DESPUÉS DE LA INTERVENCIÓN -----------------------

extremos_dummy <- resultados_dummy %>%
  filter(is.finite(residuo_estandarizado)) %>%
  arrange(desc(abs(residuo_estandarizado))) %>%
  select(
    fecha,
    exp_inf_12m,
    exp_inf_24m,
    infl_gt,
    d_exp24,
    d_exp24_lag1,
    shock_inf,
    dummy_sep2012,
    ajustado_tvp,
    residuo_tvp,
    residuo_estandarizado
  ) %>%
  slice_head(n = 15)

readr::write_csv(
  extremos_dummy,
  file.path(carpeta_salida, "11_extremos_residuales_modelo_dummy.csv")
)


# 13. GRÁFICAS -----------------------------------------------------------------

# 13.1 Trayectorias base y con intervención ------------------------------------

grafica_comparacion_beta <- comparacion_beta %>%
  select(
    fecha,
    `Modelo base` = beta_base,
    `Con pulso septiembre 2012` = beta_dummy
  ) %>%
  pivot_longer(
    cols = -fecha,
    names_to = "modelo",
    values_to = "beta"
  ) %>%
  ggplot(
    aes(
      x = fecha,
      y = beta,
      linetype = modelo
    )
  ) +
  geom_line(linewidth = 0.85) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    title = "Comparación de la trayectoria de la sensibilidad",
    subtitle = "Modelo TVP-AR(1) base y modelo con intervención en septiembre de 2012",
    x = NULL,
    y = expression(beta[t]),
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    carpeta_salida,
    "12_comparacion_trayectorias_beta.png"
  ),
  plot = grafica_comparacion_beta,
  width = 10,
  height = 6,
  dpi = 300
)


# 13.2 Trayectoria del modelo con intervención y bandas ------------------------

grafica_beta_dummy <- ggplot(
  resultados_dummy,
  aes(x = fecha, y = beta_t)
) +
  geom_ribbon(
    aes(ymin = beta_li_95, ymax = beta_ls_95),
    alpha = 0.20
  ) +
  geom_line(linewidth = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(
    xintercept = as.numeric(fecha_intervencion),
    linetype = "dotted"
  ) +
  labs(
    title = "Sensibilidad variable de las expectativas a 24 meses",
    subtitle = paste0(
      "TVP-AR(1) con pulso en septiembre de 2012; bandas del ",
      nivel_confianza * 100,
      "%"
    ),
    x = NULL,
    y = expression(beta[t])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(
    carpeta_salida,
    "13_trayectoria_beta_modelo_dummy.png"
  ),
  plot = grafica_beta_dummy,
  width = 10,
  height = 6,
  dpi = 300
)


# 13.3 Residuos estandarizados -------------------------------------------------

datos_residuos_dummy <- resultados_dummy %>%
  filter(is.finite(residuo_estandarizado))

grafica_residuos_dummy <- ggplot(
  datos_residuos_dummy,
  aes(x = fecha, y = residuo_estandarizado)
) +
  geom_line(linewidth = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted") +
  geom_vline(
    xintercept = as.numeric(fecha_intervencion),
    linetype = "dotted"
  ) +
  labs(
    title = "Residuos recursivos estandarizados",
    subtitle = "Modelo TVP-AR(1) con pulso en septiembre de 2012",
    x = NULL,
    y = "Residuo estandarizado"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(
    carpeta_salida,
    "14_residuos_estandarizados_modelo_dummy.png"
  ),
  plot = grafica_residuos_dummy,
  width = 10,
  height = 5.5,
  dpi = 300
)


# 13.4 ACF y PACF --------------------------------------------------------------

residuos_dummy_validos <- datos_residuos_dummy$residuo_estandarizado

png(
  filename = file.path(
    carpeta_salida,
    "15_acf_residuos_modelo_dummy.png"
  ),
  width = 1800,
  height = 1100,
  res = 180
)
acf(
  residuos_dummy_validos,
  lag.max = 24,
  main = "ACF de residuos: modelo con pulso septiembre 2012"
)
dev.off()

png(
  filename = file.path(
    carpeta_salida,
    "16_pacf_residuos_modelo_dummy.png"
  ),
  width = 1800,
  height = 1100,
  res = 180
)
pacf(
  residuos_dummy_validos,
  lag.max = 24,
  main = "PACF de residuos: modelo con pulso septiembre 2012"
)
dev.off()


# 13.5 Ajuste observado frente a estimado --------------------------------------

grafica_ajuste <- resultados_dummy %>%
  select(
    fecha,
    observado = d_exp24,
    ajustado = ajustado_tvp
  ) %>%
  pivot_longer(
    cols = c(observado, ajustado),
    names_to = "serie",
    values_to = "valor"
  ) %>%
  ggplot(
    aes(
      x = fecha,
      y = valor,
      linetype = serie
    )
  ) +
  geom_line(linewidth = 0.7) +
  geom_vline(
    xintercept = as.numeric(fecha_intervencion),
    linetype = "dotted"
  ) +
  labs(
    title = "Cambio observado y ajustado de las expectativas a 24 meses",
    subtitle = "Modelo TVP-AR(1) con pulso en septiembre de 2012",
    x = NULL,
    y = expression(Delta * e[t]^24),
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title.position = "plot",
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    carpeta_salida,
    "17_ajuste_modelo_dummy.png"
  ),
  plot = grafica_ajuste,
  width = 10,
  height = 6,
  dpi = 300
)


# 14. RESUMEN FINAL EN CONSOLA -------------------------------------------------

cat("\n")
cat("====================================================================\n")
cat("TVP-KALMAN AR(1): ROBUSTEZ CON PULSO EN SEPTIEMBRE DE 2012\n")
cat("====================================================================\n")
cat("Archivo de entrada: ", archivo_entrada, "\n", sep = "")
cat(
  "Muestra efectiva: ",
  format(min(bd_modelo$fecha), "%Y-%m"),
  " a ",
  format(max(bd_modelo$fecha), "%Y-%m"),
  "\n",
  sep = ""
)
cat("Número de observaciones: ", nrow(bd_modelo), "\n", sep = "")

cat("\nMODELO BASE\n")
cat(
  "Convergencia: ",
  ajuste_base$mejor$convergencia,
  "\n",
  sep = ""
)
cat(
  "Log-verosimilitud: ",
  round(loglik_base, 5),
  "\n",
  sep = ""
)
cat("H: ", format(H_base, scientific = TRUE), "\n", sep = "")
cat("Q_beta: ", format(Q_base, scientific = TRUE), "\n", sep = "")
cat(
  "Beta inicial/final: ",
  round(first(resultados_base$beta_t), 5),
  " / ",
  round(last(resultados_base$beta_t), 5),
  "\n",
  sep = ""
)

cat("\nMODELO CON PULSO SEPTIEMBRE 2012\n")
cat(
  "Convergencia: ",
  ajuste_dummy$mejor$convergencia,
  "\n",
  sep = ""
)
cat(
  "Log-verosimilitud: ",
  round(loglik_dummy, 5),
  "\n",
  sep = ""
)
cat(
  "Gamma septiembre 2012: ",
  round(gamma_dummy, 6),
  "\n",
  sep = ""
)
cat(
  "p-valor gamma: ",
  format(
    tabla_coeficientes_dummy$p_valor[
      tabla_coeficientes_dummy$parametro == "gamma_sep2012"
    ],
    scientific = TRUE
  ),
  "\n",
  sep = ""
)
cat("H: ", format(H_dummy, scientific = TRUE), "\n", sep = "")
cat("Q_beta: ", format(Q_dummy, scientific = TRUE), "\n", sep = "")
cat(
  "Beta inicial/final: ",
  round(first(resultados_dummy$beta_t), 5),
  " / ",
  round(last(resultados_dummy$beta_t), 5),
  "\n",
  sep = ""
)

cat("\nCOMPARACIÓN\n")
cat(
  "LR de la intervención: ",
  round(LR_dummy, 5),
  "\n",
  sep = ""
)
cat(
  "p-valor LR: ",
  format(p_LR_dummy, scientific = TRUE),
  "\n",
  sep = ""
)
cat(
  "Correlación entre trayectorias beta: ",
  round(
    resumen_comparacion_beta$correlacion_trayectorias,
    6
  ),
  "\n",
  sep = ""
)
cat(
  "Máxima diferencia absoluta entre beta_t: ",
  round(
    resumen_comparacion_beta$diferencia_absoluta_maxima,
    6
  ),
  "\n",
  sep = ""
)

cat("\nDiagnósticos comparados:\n")
print(tabla_diagnosticos)

cat("\nResultados guardados en:\n", carpeta_salida, "\n", sep = "")
cat("====================================================================\n")


