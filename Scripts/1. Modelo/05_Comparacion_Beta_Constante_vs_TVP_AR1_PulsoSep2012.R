# ==============================================================================
# OBJETIVO 2 - SCRIPT 05
# SENSIBILIDAD CONSTANTE VS. VARIANTE CON PULSO EN SEPTIEMBRE DE 2012
#
# MODELO RESTRINGIDO: SENSIBILIDAD CONSTANTE
#   Delta e24_t = alpha
#                 + rho * Delta e24_(t-1)
#                 + gamma * D_sep2012,t
#                 + beta * shock_inf_t
#                 + epsilon_t
#   Q_beta = 0
#
# MODELO NO RESTRINGIDO: SENSIBILIDAD VARIANTE
#   Delta e24_t = alpha
#                 + rho * Delta e24_(t-1)
#                 + gamma * D_sep2012,t
#                 + beta_t * shock_inf_t
#                 + epsilon_t
#
#   beta_t = beta_(t-1) + eta_beta,t
#   Q_beta > 0
#
# Propósito:
#   Verificar si la evidencia a favor de beta_t variante permanece después
#   de controlar el valor extraordinario de septiembre de 2012.
#
# IMPORTANTE:
# - La dummy es un pulso aditivo, no un cambio metodológico.
# - Ambos modelos usan la misma muestra, variables e inicialización difusa.
# - La única diferencia estructural es Q_beta = 0 frente a Q_beta > 0.
# - Como H0: Q_beta = 0 está en la frontera del espacio paramétrico, se
#   reportan el p-valor convencional y una aproximación de mezcla 50:50.
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
      paste0("- ", archivos_candidatos, collapse = "\n"),
      "\n\nDirectorio actual: ",
      getwd()
    )
  )
}

carpeta_salida <- file.path(
  getwd(),
  "resultados_objetivo2_constante_vs_tvp_ar1_pulso_sep2012"
)

fecha_inicio <- as.Date("2011-01-01")
fecha_pulso <- as.Date("2012-09-01")

nivel_confianza <- 0.95
z_critico <- qnorm(1 - (1 - nivel_confianza) / 2)

limite_inferior_varianza <- 1e-12

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
    dummy_sep2012 = as.numeric(fecha == fecha_pulso)
  ) %>%
  filter(fecha >= fecha_inicio) %>%
  filter(
    is.finite(d_exp24),
    is.finite(d_exp24_lag1),
    is.finite(shock_inf),
    is.finite(dummy_sep2012)
  )

if (nrow(bd_modelo) < 36) {
  stop("La muestra efectiva tiene menos de 36 observaciones.")
}

if (anyDuplicated(bd_modelo$fecha) > 0) {
  stop("Existen fechas duplicadas en la muestra efectiva.")
}

if (sum(bd_modelo$dummy_sep2012) != 1) {
  stop(
    "La dummy de septiembre de 2012 debe contener exactamente un valor igual a 1."
  )
}


# 4. FUNCIONES DE CONSTRUCCIÓN -------------------------------------------------

# 4.1 Matriz de observación común ----------------------------------------------

crear_Z <- function(y_lag1, dummy, shock) {

  n <- length(shock)

  Z_t <- array(
    NA_real_,
    dim = c(1, 4, n)
  )

  # Estados:
  # 1 alpha
  # 2 rho
  # 3 gamma_sep2012
  # 4 beta o beta_t
  Z_t[1, 1, ] <- 1
  Z_t[1, 2, ] <- y_lag1
  Z_t[1, 3, ] <- dummy
  Z_t[1, 4, ] <- shock

  Z_t
}


# 4.2 Modelo con sensibilidad constante ---------------------------------------

construir_modelo_constante <- function(
    y,
    y_lag1,
    dummy,
    shock,
    estado_inicial
) {

  n <- length(y)
  Z_t <- crear_Z(y_lag1, dummy, shock)

  # Todos los estados son constantes:
  # alpha, rho, gamma y beta.
  # R = 0 y Q = 0 impiden evolución estocástica.
  KFAS::SSModel(
    y ~ -1 +
      SSMcustom(
        Z = Z_t,
        T = diag(4),
        R = matrix(0, nrow = 4, ncol = 1),
        Q = matrix(0, nrow = 1, ncol = 1),
        n = n,
        a1 = estado_inicial,
        P1 = matrix(0, nrow = 4, ncol = 4),
        P1inf = diag(4)
      ),
    H = matrix(NA_real_, nrow = 1, ncol = 1)
  )
}


# 4.3 Modelo con sensibilidad variante ----------------------------------------

construir_modelo_tvp <- function(
    y,
    y_lag1,
    dummy,
    shock,
    estado_inicial
) {

  n <- length(y)
  Z_t <- crear_Z(y_lag1, dummy, shock)

  # alpha, rho y gamma son constantes.
  # beta_t evoluciona como paseo aleatorio.
  KFAS::SSModel(
    y ~ -1 +
      SSMcustom(
        Z = Z_t,
        T = diag(4),
        R = matrix(
          c(0, 0, 0, 1),
          nrow = 4,
          ncol = 1
        ),
        Q = matrix(NA_real_, nrow = 1, ncol = 1),
        n = n,
        a1 = estado_inicial,
        P1 = matrix(0, nrow = 4, ncol = 4),
        P1inf = diag(4)
      ),
    H = matrix(NA_real_, nrow = 1, ncol = 1)
  )
}


# 4.4 Funciones de actualización ----------------------------------------------

actualizar_H <- function(par, modelo) {

  modelo$H[1, 1, 1] <- exp(par[1])

  modelo
}


actualizar_H_Q <- function(par, modelo) {

  modelo$H[1, 1, 1] <- exp(par[1])
  modelo$Q[1, 1, 1] <- exp(par[2])

  modelo
}


# 5. FUNCIONES DE ESTIMACIÓN ---------------------------------------------------

estimar_constante <- function(
    modelo_sin_estimar,
    var_residuo
) {

  inicios_H <- log(
    pmax(
      c(
        var_residuo * 0.25,
        var_residuo * 0.50,
        var_residuo,
        var_residuo * 2,
        var_residuo * 4
      ),
      limite_inferior_varianza
    )
  )

  ajustes <- lapply(
    seq_along(inicios_H),
    function(i) {

      intento <- tryCatch(
        KFAS::fitSSM(
          model = modelo_sin_estimar,
          inits = inicios_H[i],
          updatefn = actualizar_H,
          method = "L-BFGS-B",
          lower = log(limite_inferior_varianza),
          upper = log(100),
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

      ll <- tryCatch(
        as.numeric(logLik(intento$model)),
        error = function(e) NA_real_
      )

      list(
        ajuste = intento,
        intento = i,
        convergencia = intento$optim.out$convergence,
        mensaje = ifelse(
          is.null(intento$optim.out$message),
          "",
          as.character(intento$optim.out$message)
        ),
        logLik = ll
      )
    }
  )

  ajustes <- Filter(Negate(is.null), ajustes)

  if (length(ajustes) == 0) {
    stop("No fue posible estimar el modelo constante con pulso.")
  }

  tabla <- tibble::tibble(
    modelo = "Constante + pulso sep-2012",
    intento = vapply(ajustes, function(z) z$intento, numeric(1)),
    convergencia = vapply(
      ajustes,
      function(z) z$convergencia,
      numeric(1)
    ),
    log_verosimilitud = vapply(
      ajustes,
      function(z) z$logLik,
      numeric(1)
    ),
    mensaje = vapply(
      ajustes,
      function(z) z$mensaje,
      character(1)
    )
  ) %>%
    arrange(desc(log_verosimilitud))

  elegibles <- ajustes[
    vapply(
      ajustes,
      function(z) z$convergencia == 0 && is.finite(z$logLik),
      logical(1)
    )
  ]

  if (length(elegibles) == 0) {
    warning(
      "Ningún intento del modelo constante convergió normalmente. ",
      "Se utilizará la mayor log-verosimilitud disponible."
    )
    elegibles <- ajustes
  }

  mejor <- elegibles[[
    which.max(
      vapply(elegibles, function(z) z$logLik, numeric(1))
    )
  ]]

  list(
    mejor = mejor,
    modelo = mejor$ajuste$model,
    tabla_intentos = tabla
  )
}


estimar_tvp <- function(
    modelo_sin_estimar,
    var_residuo,
    var_y
) {

  lista_inicios <- list(
    log(c(var_residuo,       var_y * 1e-2)),
    log(c(var_residuo,       var_y * 1e-3)),
    log(c(var_residuo,       var_y * 1e-4)),
    log(c(var_residuo,       var_y * 1e-5)),
    log(c(var_residuo * 0.5, var_y * 1e-3)),
    log(c(var_residuo * 2.0, var_y * 1e-3)),
    log(c(var_residuo,       1e-8)),
    log(c(var_residuo,       1e-10))
  )

  ajustes <- lapply(
    seq_along(lista_inicios),
    function(i) {

      intento <- tryCatch(
        KFAS::fitSSM(
          model = modelo_sin_estimar,
          inits = lista_inicios[[i]],
          updatefn = actualizar_H_Q,
          method = "L-BFGS-B",
          lower = rep(log(limite_inferior_varianza), 2),
          upper = rep(log(100), 2),
          control = list(
            maxit = 10000,
            factr = 1e7
          )
        ),
        error = function(e) NULL
      )

      if (is.null(intento)) {
        return(NULL)
      }

      ll <- tryCatch(
        as.numeric(logLik(intento$model)),
        error = function(e) NA_real_
      )

      list(
        ajuste = intento,
        intento = i,
        convergencia = intento$optim.out$convergence,
        mensaje = ifelse(
          is.null(intento$optim.out$message),
          "",
          as.character(intento$optim.out$message)
        ),
        logLik = ll
      )
    }
  )

  ajustes <- Filter(Negate(is.null), ajustes)

  if (length(ajustes) == 0) {
    stop("No fue posible estimar el modelo TVP con pulso.")
  }

  tabla <- tibble::tibble(
    modelo = "TVP + pulso sep-2012",
    intento = vapply(ajustes, function(z) z$intento, numeric(1)),
    convergencia = vapply(
      ajustes,
      function(z) z$convergencia,
      numeric(1)
    ),
    log_verosimilitud = vapply(
      ajustes,
      function(z) z$logLik,
      numeric(1)
    ),
    mensaje = vapply(
      ajustes,
      function(z) z$mensaje,
      character(1)
    )
  ) %>%
    arrange(desc(log_verosimilitud))

  elegibles <- ajustes[
    vapply(
      ajustes,
      function(z) z$convergencia == 0 && is.finite(z$logLik),
      logical(1)
    )
  ]

  if (length(elegibles) == 0) {
    warning(
      "Ningún intento del modelo TVP convergió normalmente. ",
      "Se utilizará la mayor log-verosimilitud disponible."
    )
    elegibles <- ajustes
  }

  mejor <- elegibles[[
    which.max(
      vapply(elegibles, function(z) z$logLik, numeric(1))
    )
  ]]

  list(
    mejor = mejor,
    modelo = mejor$ajuste$model,
    tabla_intentos = tabla
  )
}


# 6. FUNCIONES DE EXTRACCIÓN ---------------------------------------------------

extraer_estado_constante <- function(kfs, indice_estado) {

  estimaciones <- as.numeric(kfs$alphahat[, indice_estado])
  varianzas <- as.numeric(kfs$V[indice_estado, indice_estado, ])

  estimacion <- mean(
    estimaciones[is.finite(estimaciones)],
    na.rm = TRUE
  )

  error_estandar <- mean(
    sqrt(pmax(varianzas[is.finite(varianzas)], 0)),
    na.rm = TRUE
  )

  c(
    estimacion = estimacion,
    error_estandar = error_estandar
  )
}


extraer_residuos_estandarizados <- function(kfs, n_esperado) {

  residuos <- tryCatch(
    as.numeric(
      stats::rstandard(
        kfs,
        type = "recursive"
      )
    ),
    error = function(e) rep(NA_real_, n_esperado)
  )

  if (length(residuos) != n_esperado) {
    residuos <- rep(NA_real_, n_esperado)
  }

  residuos
}


extraer_innovaciones <- function(kfs, n_esperado) {

  innovaciones <- tryCatch(
    as.numeric(kfs$v),
    error = function(e) numeric(0)
  )

  if (length(innovaciones) != n_esperado) {
    innovaciones <- tryCatch(
      as.numeric(
        residuals(
          kfs,
          type = "recursive",
          standardize = FALSE
        )
      ),
      error = function(e) numeric(0)
    )
  }

  if (length(innovaciones) != n_esperado) {
    innovaciones <- rep(NA_real_, n_esperado)
  }

  innovaciones
}


extraer_varianza_innovaciones <- function(kfs, n_esperado) {

  F_t <- tryCatch(
    as.numeric(kfs$F),
    error = function(e) numeric(0)
  )

  if (length(F_t) != n_esperado) {
    F_t <- rep(NA_real_, n_esperado)
  }

  F_t
}


calcular_metricas_prediccion <- function(
    y,
    innovaciones,
    F_t,
    nombre_modelo
) {

  prediccion <- y - innovaciones

  validos <- is.finite(y) &
    is.finite(innovaciones) &
    is.finite(prediccion)

  if (sum(validos) == 0) {
    return(
      list(
        tabla = tibble::tibble(
          modelo = nombre_modelo,
          observaciones_validas = 0,
          RMSE_un_paso = NA_real_,
          MAE_un_paso = NA_real_,
          MedAE_un_paso = NA_real_,
          correlacion_observado_predicho = NA_real_,
          perdida_logaritmica_media = NA_real_
        ),
        predicciones = tibble::tibble(
          observado = y,
          predicho_un_paso = prediccion,
          error_prediccion = innovaciones,
          varianza_prediccion = F_t
        )
      )
    )
  }

  validos_log <- validos &
    is.finite(F_t) &
    F_t > 0

  perdida_log <- if (sum(validos_log) > 0) {
    mean(
      0.5 * (
        log(2 * pi) +
          log(F_t[validos_log]) +
          innovaciones[validos_log]^2 / F_t[validos_log]
      )
    )
  } else {
    NA_real_
  }

  tabla <- tibble::tibble(
    modelo = nombre_modelo,
    observaciones_validas = sum(validos),
    RMSE_un_paso = sqrt(
      mean(innovaciones[validos]^2)
    ),
    MAE_un_paso = mean(
      abs(innovaciones[validos])
    ),
    MedAE_un_paso = median(
      abs(innovaciones[validos])
    ),
    correlacion_observado_predicho = suppressWarnings(
      cor(
        y[validos],
        prediccion[validos]
      )
    ),
    perdida_logaritmica_media = perdida_log
  )

  predicciones <- tibble::tibble(
    observado = y,
    predicho_un_paso = prediccion,
    error_prediccion = innovaciones,
    varianza_prediccion = F_t
  )

  list(
    tabla = tabla,
    predicciones = predicciones
  )
}


calcular_diagnosticos <- function(
    residuos,
    nombre_modelo
) {

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
    bind_rows()

  media_r <- mean(residuos)
  desv_r <- sd(residuos)

  if (!is.finite(desv_r) || desv_r <= 0) {

    tabla_jb <- tibble::tibble(
      modelo = nombre_modelo,
      diagnostico = "Jarque-Bera aproximado",
      rezago = NA_real_,
      estadistico = NA_real_,
      grados_libertad = 2,
      p_valor = NA_real_,
      conclusion_5pct = "No disponible"
    )

  } else {

    asimetria <- mean(
      ((residuos - media_r) / desv_r)^3
    )

    curtosis <- mean(
      ((residuos - media_r) / desv_r)^4
    )

    JB <- length(residuos) / 6 *
      (
        asimetria^2 +
          ((curtosis - 3)^2) / 4
      )

    p_JB <- pchisq(
      JB,
      df = 2,
      lower.tail = FALSE
    )

    tabla_jb <- tibble::tibble(
      modelo = nombre_modelo,
      diagnostico = "Jarque-Bera aproximado",
      rezago = NA_real_,
      estadistico = JB,
      grados_libertad = 2,
      p_valor = p_JB,
      conclusion_5pct = ifelse(
        p_JB < 0.05,
        "Se rechaza normalidad",
        "No se rechaza normalidad"
      )
    )
  }

  bind_rows(
    tabla_lb,
    tabla_jb
  )
}


# 7. INICIALIZACIÓN COMÚN ------------------------------------------------------

modelo_ols_inicial <- lm(
  d_exp24 ~
    d_exp24_lag1 +
    dummy_sep2012 +
    shock_inf,
  data = bd_modelo
)

coef_ols <- coef(modelo_ols_inicial)

if (any(!is.finite(coef_ols))) {
  stop("El modelo OLS de inicialización produjo coeficientes no finitos.")
}

estado_inicial <- c(
  alpha = unname(coef_ols["(Intercept)"]),
  rho = unname(coef_ols["d_exp24_lag1"]),
  gamma = unname(coef_ols["dummy_sep2012"]),
  beta = unname(coef_ols["shock_inf"])
)

var_residuo_ols <- max(
  var(residuals(modelo_ols_inicial), na.rm = TRUE),
  1e-8
)

var_y <- max(
  var(bd_modelo$d_exp24, na.rm = TRUE),
  1e-8
)


# 8. CONSTRUCCIÓN Y ESTIMACIÓN -------------------------------------------------

modelo_constante_sin_estimar <- construir_modelo_constante(
  y = bd_modelo$d_exp24,
  y_lag1 = bd_modelo$d_exp24_lag1,
  dummy = bd_modelo$dummy_sep2012,
  shock = bd_modelo$shock_inf,
  estado_inicial = estado_inicial
)

modelo_tvp_sin_estimar <- construir_modelo_tvp(
  y = bd_modelo$d_exp24,
  y_lag1 = bd_modelo$d_exp24_lag1,
  dummy = bd_modelo$dummy_sep2012,
  shock = bd_modelo$shock_inf,
  estado_inicial = estado_inicial
)

ajuste_constante <- estimar_constante(
  modelo_sin_estimar = modelo_constante_sin_estimar,
  var_residuo = var_residuo_ols
)

ajuste_tvp <- estimar_tvp(
  modelo_sin_estimar = modelo_tvp_sin_estimar,
  var_residuo = var_residuo_ols,
  var_y = var_y
)

tabla_intentos <- bind_rows(
  ajuste_constante$tabla_intentos,
  ajuste_tvp$tabla_intentos
)

readr::write_csv(
  tabla_intentos,
  file.path(carpeta_salida, "01_intentos_estimacion.csv")
)


# 9. FILTRO Y SUAVIZADOR -------------------------------------------------------

kfs_constante <- KFAS::KFS(
  ajuste_constante$modelo,
  filtering = c("state", "mean"),
  smoothing = c("state", "mean", "disturbance")
)

kfs_tvp <- KFAS::KFS(
  ajuste_tvp$modelo,
  filtering = c("state", "mean"),
  smoothing = c("state", "mean", "disturbance")
)


# 10. COEFICIENTES DEL MODELO CONSTANTE ---------------------------------------

alpha_constante <- extraer_estado_constante(kfs_constante, 1)
rho_constante <- extraer_estado_constante(kfs_constante, 2)
gamma_constante <- extraer_estado_constante(kfs_constante, 3)
beta_constante <- extraer_estado_constante(kfs_constante, 4)

tabla_coeficientes_constante <- tibble::tibble(
  parametro = c(
    "alpha",
    "rho",
    "gamma_sep2012",
    "beta"
  ),
  estimacion = c(
    alpha_constante["estimacion"],
    rho_constante["estimacion"],
    gamma_constante["estimacion"],
    beta_constante["estimacion"]
  ),
  error_estandar = c(
    alpha_constante["error_estandar"],
    rho_constante["error_estandar"],
    gamma_constante["error_estandar"],
    beta_constante["error_estandar"]
  )
) %>%
  mutate(
    limite_inferior_95 =
      estimacion - z_critico * error_estandar,
    limite_superior_95 =
      estimacion + z_critico * error_estandar,
    estadistico_z =
      estimacion / error_estandar,
    p_valor =
      2 * pnorm(abs(estadistico_z), lower.tail = FALSE),
    significativo_5pct =
      p_valor < 0.05
  )

readr::write_csv(
  tabla_coeficientes_constante,
  file.path(
    carpeta_salida,
    "02_coeficientes_modelo_constante_pulso.csv"
  )
)


# 11. COEFICIENTES FIJOS DEL MODELO TVP ---------------------------------------

alpha_tvp <- extraer_estado_constante(kfs_tvp, 1)
rho_tvp <- extraer_estado_constante(kfs_tvp, 2)
gamma_tvp <- extraer_estado_constante(kfs_tvp, 3)

tabla_coeficientes_fijos_tvp <- tibble::tibble(
  parametro = c(
    "alpha",
    "rho",
    "gamma_sep2012"
  ),
  estimacion = c(
    alpha_tvp["estimacion"],
    rho_tvp["estimacion"],
    gamma_tvp["estimacion"]
  ),
  error_estandar = c(
    alpha_tvp["error_estandar"],
    rho_tvp["error_estandar"],
    gamma_tvp["error_estandar"]
  )
) %>%
  mutate(
    limite_inferior_95 =
      estimacion - z_critico * error_estandar,
    limite_superior_95 =
      estimacion + z_critico * error_estandar,
    estadistico_z =
      estimacion / error_estandar,
    p_valor =
      2 * pnorm(abs(estadistico_z), lower.tail = FALSE),
    significativo_5pct =
      p_valor < 0.05
  )

readr::write_csv(
  tabla_coeficientes_fijos_tvp,
  file.path(
    carpeta_salida,
    "03_coeficientes_fijos_modelo_tvp_pulso.csv"
  )
)


# 12. TRAYECTORIA DE BETA_T ----------------------------------------------------

beta_t <- as.numeric(kfs_tvp$alphahat[, 4])

se_beta_t <- sqrt(
  pmax(
    as.numeric(kfs_tvp$V[4, 4, ]),
    0
  )
)

residuos_constante <- extraer_residuos_estandarizados(
  kfs_constante,
  nrow(bd_modelo)
)

residuos_tvp <- extraer_residuos_estandarizados(
  kfs_tvp,
  nrow(bd_modelo)
)

beta_constante_estimado <- unname(
  beta_constante["estimacion"]
)

beta_constante_li <- tabla_coeficientes_constante %>%
  filter(parametro == "beta") %>%
  pull(limite_inferior_95)

beta_constante_ls <- tabla_coeficientes_constante %>%
  filter(parametro == "beta") %>%
  pull(limite_superior_95)

resultados_beta_tvp <- bd_modelo %>%
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
    beta_t = beta_t,
    se_beta_t = se_beta_t,
    beta_li_95 = beta_t - z_critico * se_beta_t,
    beta_ls_95 = beta_t + z_critico * se_beta_t,
    indice_anclaje = -beta_t,
    beta_constante = beta_constante_estimado,
    beta_constante_li_95 = beta_constante_li,
    beta_constante_ls_95 = beta_constante_ls,
    diferencia_tvp_menos_constante =
      beta_t - beta_constante_estimado,
    residuo_estandarizado_constante =
      residuos_constante,
    residuo_estandarizado_tvp =
      residuos_tvp
  )

readr::write_csv(
  resultados_beta_tvp,
  file.path(
    carpeta_salida,
    "04_trayectoria_beta_tvp_pulso.csv"
  )
)


# 13. VARIANZAS ESTIMADAS ------------------------------------------------------

H_constante <- as.numeric(
  ajuste_constante$modelo$H[1, 1, 1]
)

H_tvp <- as.numeric(
  ajuste_tvp$modelo$H[1, 1, 1]
)

Q_beta_tvp <- as.numeric(
  ajuste_tvp$modelo$Q[1, 1, 1]
)

tabla_varianzas <- tibble::tibble(
  modelo = c(
    "Constante + pulso sep-2012",
    "TVP + pulso sep-2012"
  ),
  H = c(
    H_constante,
    H_tvp
  ),
  Q_beta = c(
    0,
    Q_beta_tvp
  ),
  Q_beta_sobre_H = c(
    0,
    Q_beta_tvp / H_tvp
  ),
  desviacion_estandar_observacion =
    sqrt(c(H_constante, H_tvp)),
  desviacion_estandar_innovacion_beta = c(
    0,
    sqrt(Q_beta_tvp)
  ),
  razon_Q_sobre_limite_inferior = c(
    NA_real_,
    Q_beta_tvp / limite_inferior_varianza
  ),
  Q_beta_en_limite_numerico = c(
    FALSE,
    Q_beta_tvp <= limite_inferior_varianza * 10
  )
)

readr::write_csv(
  tabla_varianzas,
  file.path(
    carpeta_salida,
    "05_varianzas_estimadas.csv"
  )
)


# 14. PREDICCIÓN RECURSIVA DE UN PASO -----------------------------------------

innovaciones_constante <- extraer_innovaciones(
  kfs_constante,
  nrow(bd_modelo)
)

innovaciones_tvp <- extraer_innovaciones(
  kfs_tvp,
  nrow(bd_modelo)
)

F_constante <- extraer_varianza_innovaciones(
  kfs_constante,
  nrow(bd_modelo)
)

F_tvp <- extraer_varianza_innovaciones(
  kfs_tvp,
  nrow(bd_modelo)
)

metricas_constante <- calcular_metricas_prediccion(
  y = bd_modelo$d_exp24,
  innovaciones = innovaciones_constante,
  F_t = F_constante,
  nombre_modelo = "Constante + pulso sep-2012"
)

metricas_tvp <- calcular_metricas_prediccion(
  y = bd_modelo$d_exp24,
  innovaciones = innovaciones_tvp,
  F_t = F_tvp,
  nombre_modelo = "TVP + pulso sep-2012"
)

tabla_metricas_prediccion <- bind_rows(
  metricas_constante$tabla,
  metricas_tvp$tabla
)

predicciones_comparadas <- bd_modelo %>%
  select(
    fecha,
    observado = d_exp24
  ) %>%
  bind_cols(
    metricas_constante$predicciones %>%
      select(
        predicho_constante = predicho_un_paso,
        error_constante = error_prediccion,
        varianza_prediccion_constante =
          varianza_prediccion
      ),
    metricas_tvp$predicciones %>%
      select(
        predicho_tvp = predicho_un_paso,
        error_tvp = error_prediccion,
        varianza_prediccion_tvp =
          varianza_prediccion
      )
  )

readr::write_csv(
  tabla_metricas_prediccion,
  file.path(
    carpeta_salida,
    "06_metricas_prediccion_un_paso.csv"
  )
)

readr::write_csv(
  predicciones_comparadas,
  file.path(
    carpeta_salida,
    "07_predicciones_un_paso.csv"
  )
)


# 15. COMPARACIÓN DE MODELOS ---------------------------------------------------

loglik_constante <- as.numeric(
  logLik(ajuste_constante$modelo)
)

loglik_tvp <- as.numeric(
  logLik(ajuste_tvp$modelo)
)

# Parámetros estructurales:
# Constante + pulso: alpha, rho, gamma, beta, H = 5
# TVP + pulso: alpha, rho, gamma, beta inicial, H, Q_beta = 6
k_constante <- 5
k_tvp <- 6
n_obs <- nrow(bd_modelo)

AIC_constante <-
  -2 * loglik_constante + 2 * k_constante

AIC_tvp <-
  -2 * loglik_tvp + 2 * k_tvp

BIC_constante <-
  -2 * loglik_constante + log(n_obs) * k_constante

BIC_tvp <-
  -2 * loglik_tvp + log(n_obs) * k_tvp

AICc_constante <-
  AIC_constante +
  (
    2 * k_constante * (k_constante + 1)
  ) /
  (
    n_obs - k_constante - 1
  )

AICc_tvp <-
  AIC_tvp +
  (
    2 * k_tvp * (k_tvp + 1)
  ) /
  (
    n_obs - k_tvp - 1
  )

tabla_comparacion_modelos <- tibble::tibble(
  modelo = c(
    "Constante + pulso sep-2012",
    "TVP + pulso sep-2012"
  ),
  log_verosimilitud = c(
    loglik_constante,
    loglik_tvp
  ),
  parametros_contados = c(
    k_constante,
    k_tvp
  ),
  AIC = c(
    AIC_constante,
    AIC_tvp
  ),
  AICc = c(
    AICc_constante,
    AICc_tvp
  ),
  BIC = c(
    BIC_constante,
    BIC_tvp
  ),
  H = c(
    H_constante,
    H_tvp
  ),
  Q_beta = c(
    0,
    Q_beta_tvp
  ),
  RMSE_un_paso = c(
    metricas_constante$tabla$RMSE_un_paso,
    metricas_tvp$tabla$RMSE_un_paso
  ),
  MAE_un_paso = c(
    metricas_constante$tabla$MAE_un_paso,
    metricas_tvp$tabla$MAE_un_paso
  ),
  perdida_logaritmica_media = c(
    metricas_constante$tabla$perdida_logaritmica_media,
    metricas_tvp$tabla$perdida_logaritmica_media
  )
) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    peso_Akaike =
      exp(-0.5 * delta_AIC) /
      sum(exp(-0.5 * delta_AIC)),
    delta_BIC = BIC - min(BIC),
    peso_BIC_aproximado =
      exp(-0.5 * delta_BIC) /
      sum(exp(-0.5 * delta_BIC))
  )

readr::write_csv(
  tabla_comparacion_modelos,
  file.path(
    carpeta_salida,
    "08_comparacion_modelos_constante_vs_tvp_pulso.csv"
  )
)


# 16. PRUEBA DE RAZÓN DE VEROSIMILITUD ----------------------------------------

estadistico_LR <- max(
  0,
  2 * (
    loglik_tvp -
      loglik_constante
  )
)

p_LR_convencional <- pchisq(
  estadistico_LR,
  df = 1,
  lower.tail = FALSE
)

# Aproximación de mezcla 50:50:
# 0.5 * chi-cuadrado(0) + 0.5 * chi-cuadrado(1)
p_LR_frontera <- ifelse(
  estadistico_LR > 0,
  0.5 * p_LR_convencional,
  1
)

tabla_prueba_LR <- tibble::tibble(
  hipotesis_nula =
    "H0: Q_beta = 0; sensibilidad constante con pulso",
  hipotesis_alternativa =
    "H1: Q_beta > 0; sensibilidad variante con pulso",
  estadistico_LR = estadistico_LR,
  grados_libertad_referencia = 1,
  p_valor_convencional_chi2_1 =
    p_LR_convencional,
  p_valor_mezcla_frontera_50_50 =
    p_LR_frontera,
  conclusion_frontera_5pct = ifelse(
    p_LR_frontera < 0.05,
    "Se rechaza H0; evidencia a favor de sensibilidad variante",
    "No se rechaza H0; evidencia insuficiente para variar beta"
  ),
  advertencia = paste(
    "Q_beta está en la frontera bajo H0.",
    "El p-valor de mezcla es una aproximación;",
    "debe evaluarse junto con AIC, BIC y predicción."
  )
)

readr::write_csv(
  tabla_prueba_LR,
  file.path(
    carpeta_salida,
    "09_prueba_LR_variacion_beta_pulso.csv"
  )
)


# 17. DIAGNÓSTICOS RESIDUALES --------------------------------------------------

diagnosticos_constante <- calcular_diagnosticos(
  residuos_constante,
  "Constante + pulso sep-2012"
)

diagnosticos_tvp <- calcular_diagnosticos(
  residuos_tvp,
  "TVP + pulso sep-2012"
)

tabla_diagnosticos <- bind_rows(
  diagnosticos_constante,
  diagnosticos_tvp
)

readr::write_csv(
  tabla_diagnosticos,
  file.path(
    carpeta_salida,
    "10_diagnosticos_residuales_comparados.csv"
  )
)


# 18. RESUMEN DE SENSIBILIDAD --------------------------------------------------

porcentaje_beta_positivo <- mean(
  resultados_beta_tvp$beta_t > 0,
  na.rm = TRUE
) * 100

porcentaje_beta_significativo_positivo <- mean(
  resultados_beta_tvp$beta_li_95 > 0,
  na.rm = TRUE
) * 100

porcentaje_intervalo_incluye_cero <- mean(
  resultados_beta_tvp$beta_li_95 <= 0 &
    resultados_beta_tvp$beta_ls_95 >= 0,
  na.rm = TRUE
) * 100

indice_maximo_beta <- which.max(
  resultados_beta_tvp$beta_t
)

indice_minimo_beta <- which.min(
  resultados_beta_tvp$beta_t
)

tabla_resumen_beta <- tibble::tibble(
  beta_constante =
    beta_constante_estimado,
  beta_constante_p_valor =
    tabla_coeficientes_constante %>%
      filter(parametro == "beta") %>%
      pull(p_valor),
  gamma_constante =
    unname(gamma_constante["estimacion"]),
  gamma_constante_p_valor =
    tabla_coeficientes_constante %>%
      filter(parametro == "gamma_sep2012") %>%
      pull(p_valor),
  gamma_tvp =
    unname(gamma_tvp["estimacion"]),
  gamma_tvp_p_valor =
    tabla_coeficientes_fijos_tvp %>%
      filter(parametro == "gamma_sep2012") %>%
      pull(p_valor),
  beta_tvp_inicial =
    first(resultados_beta_tvp$beta_t),
  beta_tvp_final =
    last(resultados_beta_tvp$beta_t),
  cambio_total_beta_tvp =
    last(resultados_beta_tvp$beta_t) -
      first(resultados_beta_tvp$beta_t),
  reduccion_porcentual_desde_inicio = ifelse(
    abs(first(resultados_beta_tvp$beta_t)) > 0,
    100 * (
      1 -
        last(resultados_beta_tvp$beta_t) /
          first(resultados_beta_tvp$beta_t)
    ),
    NA_real_
  ),
  beta_tvp_maximo =
    resultados_beta_tvp$beta_t[indice_maximo_beta],
  fecha_beta_tvp_maximo =
    resultados_beta_tvp$fecha[indice_maximo_beta],
  beta_tvp_minimo =
    resultados_beta_tvp$beta_t[indice_minimo_beta],
  fecha_beta_tvp_minimo =
    resultados_beta_tvp$fecha[indice_minimo_beta],
  porcentaje_beta_positivo =
    porcentaje_beta_positivo,
  porcentaje_beta_significativo_positivo =
    porcentaje_beta_significativo_positivo,
  porcentaje_intervalo_incluye_cero =
    porcentaje_intervalo_incluye_cero,
  Q_beta =
    Q_beta_tvp,
  Q_beta_sobre_H =
    Q_beta_tvp / H_tvp
)

readr::write_csv(
  tabla_resumen_beta,
  file.path(
    carpeta_salida,
    "11_resumen_sensibilidad_constante_y_tvp_pulso.csv"
  )
)


# 19. GRÁFICAS -----------------------------------------------------------------

# 19.1 Beta constante frente a beta_t ------------------------------------------

grafica_beta_comparada <- ggplot(
  resultados_beta_tvp,
  aes(x = fecha)
) +
  geom_ribbon(
    aes(
      ymin = beta_li_95,
      ymax = beta_ls_95
    ),
    alpha = 0.18
  ) +
  geom_line(
    aes(
      y = beta_t,
      linetype = "Beta variante"
    ),
    linewidth = 0.9
  ) +
  geom_hline(
    aes(
      yintercept = beta_constante,
      linetype = "Beta constante"
    ),
    linewidth = 0.85
  ) +
  geom_hline(
    yintercept = beta_constante_li,
    linetype = "dotted",
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept = beta_constante_ls,
    linetype = "dotted",
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "longdash",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = as.numeric(fecha_pulso),
    linetype = "dotted",
    linewidth = 0.55
  ) +
  labs(
    title =
      "Sensibilidad constante y variante con pulso en septiembre de 2012",
    subtitle =
      "Comparación homogénea bajo una representación de espacio de estados",
    x = NULL,
    y = expression(beta),
    linetype = NULL,
    caption = paste(
      "La banda sombreada corresponde al intervalo del 95% de beta_t.",
      "Las líneas horizontales punteadas delimitan el intervalo del beta constante."
    )
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
    "12_beta_constante_vs_tvp_pulso.png"
  ),
  plot = grafica_beta_comparada,
  width = 10,
  height = 6,
  dpi = 300
)


# 19.2 Trayectoria TVP con pulso -----------------------------------------------

grafica_beta_tvp <- ggplot(
  resultados_beta_tvp,
  aes(x = fecha, y = beta_t)
) +
  geom_ribbon(
    aes(
      ymin = beta_li_95,
      ymax = beta_ls_95
    ),
    alpha = 0.18
  ) +
  geom_line(linewidth = 0.9) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = as.numeric(fecha_pulso),
    linetype = "dotted",
    linewidth = 0.55
  ) +
  labs(
    title =
      "Trayectoria de la sensibilidad variante con intervención",
    subtitle =
      paste0(
        "Pulso en septiembre de 2012; bandas del ",
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
    "13_trayectoria_beta_tvp_pulso.png"
  ),
  plot = grafica_beta_tvp,
  width = 10,
  height = 6,
  dpi = 300
)


# 19.3 Residuos comparados ------------------------------------------------------

datos_residuos_largos <- resultados_beta_tvp %>%
  select(
    fecha,
    `Constante + pulso` =
      residuo_estandarizado_constante,
    `TVP + pulso` =
      residuo_estandarizado_tvp
  ) %>%
  pivot_longer(
    cols = -fecha,
    names_to = "modelo",
    values_to = "residuo"
  ) %>%
  filter(is.finite(residuo))

grafica_residuos <- ggplot(
  datos_residuos_largos,
  aes(
    x = fecha,
    y = residuo,
    linetype = modelo
  )
) +
  geom_line(linewidth = 0.65) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = c(-2, 2),
    linetype = "dotted"
  ) +
  geom_vline(
    xintercept = as.numeric(fecha_pulso),
    linetype = "dotted",
    linewidth = 0.55
  ) +
  labs(
    title =
      "Residuos recursivos estandarizados con intervención",
    x = NULL,
    y = "Residuo estandarizado",
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
    "14_residuos_estandarizados_comparados.png"
  ),
  plot = grafica_residuos,
  width = 10,
  height = 5.8,
  dpi = 300
)


# 19.4 Predicciones de un paso -------------------------------------------------

datos_prediccion_largos <- predicciones_comparadas %>%
  select(
    fecha,
    Observado = observado,
    `Predicho: beta constante` =
      predicho_constante,
    `Predicho: beta variante` =
      predicho_tvp
  ) %>%
  pivot_longer(
    cols = -fecha,
    names_to = "serie",
    values_to = "valor"
  ) %>%
  filter(is.finite(valor))

grafica_prediccion <- ggplot(
  datos_prediccion_largos,
  aes(
    x = fecha,
    y = valor,
    linetype = serie
  )
) +
  geom_line(linewidth = 0.65) +
  geom_vline(
    xintercept = as.numeric(fecha_pulso),
    linetype = "dotted",
    linewidth = 0.55
  ) +
  labs(
    title =
      "Cambio observado y predicciones de un paso con intervención",
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
    "15_predicciones_un_paso_comparadas.png"
  ),
  plot = grafica_prediccion,
  width = 10,
  height = 6,
  dpi = 300
)


# 19.5 ACF y PACF --------------------------------------------------------------

res_const_validos <- residuos_constante[
  is.finite(residuos_constante)
]

res_tvp_validos <- residuos_tvp[
  is.finite(residuos_tvp)
]

png(
  filename = file.path(
    carpeta_salida,
    "16_acf_residuos_constante_pulso.png"
  ),
  width = 1800,
  height = 1100,
  res = 180
)
acf(
  res_const_validos,
  lag.max = 24,
  main =
    "ACF de residuos: sensibilidad constante con pulso"
)
dev.off()

png(
  filename = file.path(
    carpeta_salida,
    "17_acf_residuos_tvp_pulso.png"
  ),
  width = 1800,
  height = 1100,
  res = 180
)
acf(
  res_tvp_validos,
  lag.max = 24,
  main =
    "ACF de residuos: sensibilidad variante con pulso"
)
dev.off()

png(
  filename = file.path(
    carpeta_salida,
    "18_pacf_residuos_constante_pulso.png"
  ),
  width = 1800,
  height = 1100,
  res = 180
)
pacf(
  res_const_validos,
  lag.max = 24,
  main =
    "PACF de residuos: sensibilidad constante con pulso"
)
dev.off()

png(
  filename = file.path(
    carpeta_salida,
    "19_pacf_residuos_tvp_pulso.png"
  ),
  width = 1800,
  height = 1100,
  res = 180
)
pacf(
  res_tvp_validos,
  lag.max = 24,
  main =
    "PACF de residuos: sensibilidad variante con pulso"
)
dev.off()


# 20. GUARDAR MODELOS ----------------------------------------------------------

saveRDS(
  ajuste_constante$modelo,
  file.path(
    carpeta_salida,
    "20_modelo_constante_pulso_kfas.rds"
  )
)

saveRDS(
  ajuste_tvp$modelo,
  file.path(
    carpeta_salida,
    "21_modelo_tvp_pulso_kfas.rds"
  )
)


# 21. RESUMEN EJECUTIVO --------------------------------------------------------

mejor_AIC <- tabla_comparacion_modelos$modelo[
  which.min(tabla_comparacion_modelos$AIC)
]

mejor_BIC <- tabla_comparacion_modelos$modelo[
  which.min(tabla_comparacion_modelos$BIC)
]

formato_num <- function(x, digitos = 5) {
  formatC(
    x,
    digits = digitos,
    format = "f"
  )
}

formato_cientifico <- function(x, digitos = 4) {
  formatC(
    x,
    digits = digitos,
    format = "e"
  )
}

lineas_resumen <- c(
  "",
  "======================================================================",
  "COMPARACIÓN CON PULSO: SENSIBILIDAD CONSTANTE VS. VARIANTE",
  "======================================================================",
  paste0("Archivo: ", archivo_entrada),
  paste0(
    "Muestra: ",
    format(min(bd_modelo$fecha), "%Y-%m"),
    " a ",
    format(max(bd_modelo$fecha), "%Y-%m")
  ),
  paste0("Observaciones: ", n_obs),
  paste0(
    "Intervención: pulso en ",
    format(fecha_pulso, "%Y-%m")
  ),
  "",
  "MODELO CONSTANTE + PULSO",
  paste0(
    "Convergencia: ",
    ajuste_constante$mejor$convergencia
  ),
  paste0(
    "Log-verosimilitud: ",
    formato_num(loglik_constante)
  ),
  paste0(
    "alpha: ",
    formato_num(alpha_constante["estimacion"])
  ),
  paste0(
    "rho: ",
    formato_num(rho_constante["estimacion"]),
    " | p = ",
    formato_cientifico(
      tabla_coeficientes_constante$p_valor[
        tabla_coeficientes_constante$parametro == "rho"
      ]
    )
  ),
  paste0(
    "gamma sep-2012: ",
    formato_num(gamma_constante["estimacion"]),
    " | p = ",
    formato_cientifico(
      tabla_coeficientes_constante$p_valor[
        tabla_coeficientes_constante$parametro ==
          "gamma_sep2012"
      ]
    )
  ),
  paste0(
    "beta constante: ",
    formato_num(beta_constante["estimacion"]),
    " | p = ",
    formato_cientifico(
      tabla_coeficientes_constante$p_valor[
        tabla_coeficientes_constante$parametro == "beta"
      ]
    )
  ),
  paste0(
    "H: ",
    formato_cientifico(H_constante)
  ),
  "",
  "MODELO TVP + PULSO",
  paste0(
    "Convergencia: ",
    ajuste_tvp$mejor$convergencia
  ),
  paste0(
    "Log-verosimilitud: ",
    formato_num(loglik_tvp)
  ),
  paste0(
    "alpha: ",
    formato_num(alpha_tvp["estimacion"])
  ),
  paste0(
    "rho: ",
    formato_num(rho_tvp["estimacion"]),
    " | p = ",
    formato_cientifico(
      tabla_coeficientes_fijos_tvp$p_valor[
        tabla_coeficientes_fijos_tvp$parametro == "rho"
      ]
    )
  ),
  paste0(
    "gamma sep-2012: ",
    formato_num(gamma_tvp["estimacion"]),
    " | p = ",
    formato_cientifico(
      tabla_coeficientes_fijos_tvp$p_valor[
        tabla_coeficientes_fijos_tvp$parametro ==
          "gamma_sep2012"
      ]
    )
  ),
  paste0(
    "H: ",
    formato_cientifico(H_tvp)
  ),
  paste0(
    "Q_beta: ",
    formato_cientifico(Q_beta_tvp)
  ),
  paste0(
    "Q_beta/H: ",
    formato_cientifico(
      Q_beta_tvp / H_tvp
    )
  ),
  paste0(
    "Beta inicial/final: ",
    formato_num(first(resultados_beta_tvp$beta_t)),
    " / ",
    formato_num(last(resultados_beta_tvp$beta_t))
  ),
  "",
  "COMPARACIÓN FORMAL",
  paste0(
    "LR: ",
    formato_num(estadistico_LR)
  ),
  paste0(
    "p convencional chi2(1): ",
    formato_cientifico(p_LR_convencional)
  ),
  paste0(
    "p mezcla de frontera 50:50: ",
    formato_cientifico(p_LR_frontera)
  ),
  paste0(
    "Mejor AIC: ",
    mejor_AIC
  ),
  paste0(
    "Mejor BIC: ",
    mejor_BIC
  ),
  paste0(
    "RMSE constante/TVP: ",
    formato_num(
      metricas_constante$tabla$RMSE_un_paso
    ),
    " / ",
    formato_num(
      metricas_tvp$tabla$RMSE_un_paso
    )
  ),
  paste0(
    "MAE constante/TVP: ",
    formato_num(
      metricas_constante$tabla$MAE_un_paso
    ),
    " / ",
    formato_num(
      metricas_tvp$tabla$MAE_un_paso
    )
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
