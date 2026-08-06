# ==============================================================================
# OBJETIVO 2 - COMPARACIÓN FORMAL DE SENSIBILIDAD CONSTANTE Y VARIANTE
# Modelos AR(1) estimados bajo una representación homogénea en KFAS
#
# MODELO RESTRINGIDO: SENSIBILIDAD CONSTANTE
#   Delta e24_t = alpha + rho * Delta e24_(t-1) + beta * shock_inf_t + epsilon_t
#   Q_beta = 0
#
# MODELO NO RESTRINGIDO: SENSIBILIDAD VARIANTE EN EL TIEMPO
#   Delta e24_t = alpha + rho * Delta e24_(t-1) + beta_t * shock_inf_t + epsilon_t
#   beta_t = beta_(t-1) + eta_beta,t
#   Q_beta > 0
#
# IMPORTANTE:
# - Esta comparación NO incluye la intervención de septiembre de 2012.
# - Ambos modelos usan la misma muestra y la misma inicialización difusa.
# - La única diferencia estructural es si beta permanece constante o evoluciona.
# - Como H0: Q_beta = 0 está en la frontera del espacio paramétrico, se reporta:
#     a) p-valor convencional chi-cuadrado(1);
#     b) p-valor aproximado con mezcla 50:50.
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
  "resultados_objetivo2_constante_vs_tvp_ar1"
)

fecha_inicio <- as.Date("2011-01-01")
nivel_confianza <- 0.95
z_critico <- qnorm(1 - (1 - nivel_confianza) / 2)

# Límite inferior usado para Q_beta en el optimizador.
# Sirve para verificar si la estimación colapsa numéricamente hacia cero.
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
    shock_inf = infl_gt - exp12_lag1
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

if (anyDuplicated(bd_modelo$fecha) > 0) {
  stop("Existen fechas duplicadas en la muestra efectiva.")
}


# 4. FUNCIONES DE CONSTRUCCIÓN -------------------------------------------------

# 4.1 Matriz de observación común ----------------------------------------------

crear_Z <- function(y_lag1, x) {

  n <- length(x)

  Z_t <- array(
    NA_real_,
    dim = c(1, 3, n)
  )

  Z_t[1, 1, ] <- 1
  Z_t[1, 2, ] <- y_lag1
  Z_t[1, 3, ] <- x

  Z_t
}


# 4.2 Modelo con sensibilidad constante ---------------------------------------

construir_modelo_constante <- function(
    y,
    y_lag1,
    x,
    estado_inicial
) {

  n <- length(y)
  Z_t <- crear_Z(y_lag1, x)

  # Estados:
  # 1 alpha constante
  # 2 rho constante
  # 3 beta constante
  #
  # Se incluye una perturbación ficticia con carga cero para mantener
  # una estructura compatible con SSMcustom. Al ser R = 0 y Q = 0,
  # ningún estado evoluciona estocásticamente.
  KFAS::SSModel(
    y ~ -1 +
      SSMcustom(
        Z = Z_t,
        T = diag(3),
        R = matrix(0, nrow = 3, ncol = 1),
        Q = matrix(0, nrow = 1, ncol = 1),
        n = n,
        a1 = estado_inicial,
        P1 = matrix(0, nrow = 3, ncol = 3),
        P1inf = diag(3)
      ),
    H = matrix(NA_real_, nrow = 1, ncol = 1)
  )
}


# 4.3 Modelo con sensibilidad variante ----------------------------------------

construir_modelo_tvp <- function(
    y,
    y_lag1,
    x,
    estado_inicial
) {

  n <- length(y)
  Z_t <- crear_Z(y_lag1, x)

  # Estados:
  # 1 alpha constante
  # 2 rho constante
  # 3 beta_t variante como paseo aleatorio
  KFAS::SSModel(
    y ~ -1 +
      SSMcustom(
        Z = Z_t,
        T = diag(3),
        R = matrix(
          c(0, 0, 1),
          nrow = 3,
          ncol = 1
        ),
        Q = matrix(NA_real_, nrow = 1, ncol = 1),
        n = n,
        a1 = estado_inicial,
        P1 = matrix(0, nrow = 3, ncol = 3),
        P1inf = diag(3)
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
    stop("No fue posible estimar el modelo de sensibilidad constante.")
  }

  tabla <- tibble::tibble(
    modelo = "Sensibilidad constante",
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
    stop("No fue posible estimar el modelo TVP.")
  }

  tabla <- tibble::tibble(
    modelo = "Sensibilidad variante",
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
  d_exp24 ~ d_exp24_lag1 + shock_inf,
  data = bd_modelo
)

coef_ols <- coef(modelo_ols_inicial)

if (any(!is.finite(coef_ols))) {
  stop("El modelo OLS de inicialización produjo coeficientes no finitos.")
}

estado_inicial <- c(
  alpha = unname(coef_ols["(Intercept)"]),
  rho = unname(coef_ols["d_exp24_lag1"]),
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
  x = bd_modelo$shock_inf,
  estado_inicial = estado_inicial
)

modelo_tvp_sin_estimar <- construir_modelo_tvp(
  y = bd_modelo$d_exp24,
  y_lag1 = bd_modelo$d_exp24_lag1,
  x = bd_modelo$shock_inf,
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


# 10. PARÁMETROS DEL MODELO CONSTANTE -----------------------------------------

alpha_constante <- extraer_estado_constante(
  kfs_constante,
  1
)

rho_constante <- extraer_estado_constante(
  kfs_constante,
  2
)

beta_constante <- extraer_estado_constante(
  kfs_constante,
  3
)

tabla_coeficientes_constante <- tibble::tibble(
  parametro = c("alpha", "rho", "beta"),
  estimacion = c(
    alpha_constante["estimacion"],
    rho_constante["estimacion"],
    beta_constante["estimacion"]
  ),
  error_estandar = c(
    alpha_constante["error_estandar"],
    rho_constante["error_estandar"],
    beta_constante["error_estandar"]
  )
) %>%
  mutate(
    limite_inferior_95 = estimacion -
      z_critico * error_estandar,
    limite_superior_95 = estimacion +
      z_critico * error_estandar,
    estadistico_z = estimacion / error_estandar,
    p_valor = 2 * pnorm(
      abs(estadistico_z),
      lower.tail = FALSE
    ),
    significativo_5pct = p_valor < 0.05
  )

readr::write_csv(
  tabla_coeficientes_constante,
  file.path(
    carpeta_salida,
    "02_coeficientes_modelo_constante.csv"
  )
)


# 11. PARÁMETROS DEL MODELO TVP ------------------------------------------------

alpha_tvp <- extraer_estado_constante(
  kfs_tvp,
  1
)

rho_tvp <- extraer_estado_constante(
  kfs_tvp,
  2
)

tabla_coeficientes_fijos_tvp <- tibble::tibble(
  parametro = c("alpha", "rho"),
  estimacion = c(
    alpha_tvp["estimacion"],
    rho_tvp["estimacion"]
  ),
  error_estandar = c(
    alpha_tvp["error_estandar"],
    rho_tvp["error_estandar"]
  )
) %>%
  mutate(
    limite_inferior_95 = estimacion -
      z_critico * error_estandar,
    limite_superior_95 = estimacion +
      z_critico * error_estandar,
    estadistico_z = estimacion / error_estandar,
    p_valor = 2 * pnorm(
      abs(estadistico_z),
      lower.tail = FALSE
    ),
    significativo_5pct = p_valor < 0.05
  )

readr::write_csv(
  tabla_coeficientes_fijos_tvp,
  file.path(
    carpeta_salida,
    "03_coeficientes_fijos_modelo_tvp.csv"
  )
)


# 12. TRAYECTORIA DE BETA_T ----------------------------------------------------

beta_t <- as.numeric(kfs_tvp$alphahat[, 3])

se_beta_t <- sqrt(
  pmax(
    as.numeric(kfs_tvp$V[3, 3, ]),
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
    beta_t = beta_t,
    se_beta_t = se_beta_t,
    beta_li_95 = beta_t - z_critico * se_beta_t,
    beta_ls_95 = beta_t + z_critico * se_beta_t,
    indice_anclaje = -beta_t,
    beta_constante = unname(
      beta_constante["estimacion"]
    ),
    beta_constante_li_95 = tabla_coeficientes_constante %>%
      filter(parametro == "beta") %>%
      pull(limite_inferior_95),
    beta_constante_ls_95 = tabla_coeficientes_constante %>%
      filter(parametro == "beta") %>%
      pull(limite_superior_95),
    diferencia_tvp_menos_constante =
      beta_t - unname(beta_constante["estimacion"]),
    residuo_estandarizado_constante = residuos_constante,
    residuo_estandarizado_tvp = residuos_tvp
  )

readr::write_csv(
  resultados_beta_tvp,
  file.path(
    carpeta_salida,
    "04_trayectoria_beta_tvp.csv"
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
    "Sensibilidad constante",
    "Sensibilidad variante"
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
  desviacion_estandar_observacion = sqrt(
    c(H_constante, H_tvp)
  ),
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
  nombre_modelo = "Sensibilidad constante"
)

metricas_tvp <- calcular_metricas_prediccion(
  y = bd_modelo$d_exp24,
  innovaciones = innovaciones_tvp,
  F_t = F_tvp,
  nombre_modelo = "Sensibilidad variante"
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


# 15. COMPARACIÓN DE VEROSIMILITUD E INFORMACIÓN ------------------------------

loglik_constante <- as.numeric(
  logLik(ajuste_constante$modelo)
)

loglik_tvp <- as.numeric(
  logLik(ajuste_tvp$modelo)
)

# Conteo de parámetros estructurales:
# Constante: alpha, rho, beta, H = 4
# TVP: alpha, rho, beta inicial, H, Q_beta = 5
k_constante <- 4
k_tvp <- 5
n_obs <- nrow(bd_modelo)

AIC_constante <- -2 * loglik_constante +
  2 * k_constante

AIC_tvp <- -2 * loglik_tvp +
  2 * k_tvp

BIC_constante <- -2 * loglik_constante +
  log(n_obs) * k_constante

BIC_tvp <- -2 * loglik_tvp +
  log(n_obs) * k_tvp

AICc_constante <- AIC_constante +
  (
    2 * k_constante * (k_constante + 1)
  ) /
  (
    n_obs - k_constante - 1
  )

AICc_tvp <- AIC_tvp +
  (
    2 * k_tvp * (k_tvp + 1)
  ) /
  (
    n_obs - k_tvp - 1
  )

tabla_comparacion_modelos <- tibble::tibble(
  modelo = c(
    "Sensibilidad constante",
    "Sensibilidad variante"
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
    peso_Akaike = exp(-0.5 * delta_AIC) /
      sum(exp(-0.5 * delta_AIC)),
    delta_BIC = BIC - min(BIC),
    peso_BIC_aproximado = exp(-0.5 * delta_BIC) /
      sum(exp(-0.5 * delta_BIC))
  )

readr::write_csv(
  tabla_comparacion_modelos,
  file.path(
    carpeta_salida,
    "08_comparacion_modelos_constante_vs_tvp.csv"
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

# Referencia convencional, no ajustada por frontera.
p_LR_convencional <- pchisq(
  estadistico_LR,
  df = 1,
  lower.tail = FALSE
)

# Aproximación 50:50 de chi-cuadrado(0) y chi-cuadrado(1).
# Para LR > 0, la cola equivale a 0.5 * P[chi2(1) >= LR].
p_LR_frontera <- ifelse(
  estadistico_LR > 0,
  0.5 * p_LR_convencional,
  1
)

tabla_prueba_LR <- tibble::tibble(
  hipotesis_nula = "H0: Q_beta = 0; sensibilidad constante",
  hipotesis_alternativa =
    "H1: Q_beta > 0; sensibilidad variante",
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
  advertencia =
    paste(
      "Q_beta está en la frontera bajo H0.",
      "El p-valor de mezcla es una aproximación;",
      "AIC, BIC y desempeño predictivo deben evaluarse conjuntamente."
    )
)

readr::write_csv(
  tabla_prueba_LR,
  file.path(
    carpeta_salida,
    "09_prueba_LR_variacion_beta.csv"
  )
)


# 17. DIAGNÓSTICOS RESIDUALES --------------------------------------------------

diagnosticos_constante <- calcular_diagnosticos(
  residuos_constante,
  "Sensibilidad constante"
)

diagnosticos_tvp <- calcular_diagnosticos(
  residuos_tvp,
  "Sensibilidad variante"
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


# 18. RESUMEN DE BETA ----------------------------------------------------------

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
  beta_constante = unname(
    beta_constante["estimacion"]
  ),
  beta_constante_p_valor =
    tabla_coeficientes_constante %>%
      filter(parametro == "beta") %>%
      pull(p_valor),
  beta_tvp_inicial = first(
    resultados_beta_tvp$beta_t
  ),
  beta_tvp_final = last(
    resultados_beta_tvp$beta_t
  ),
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
  Q_beta = Q_beta_tvp,
  Q_beta_sobre_H = Q_beta_tvp / H_tvp
)

readr::write_csv(
  tabla_resumen_beta,
  file.path(
    carpeta_salida,
    "11_resumen_sensibilidad_constante_y_tvp.csv"
  )
)


# 19. GRÁFICAS -----------------------------------------------------------------

# 19.1 Beta constante frente a trayectoria TVP --------------------------------

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
    yintercept =
      unique(
        resultados_beta_tvp$beta_constante_li_95
      ),
    linetype = "dotted",
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept =
      unique(
        resultados_beta_tvp$beta_constante_ls_95
      ),
    linetype = "dotted",
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "longdash",
    linewidth = 0.5
  ) +
  labs(
    title =
      "Sensibilidad constante y sensibilidad variante en el tiempo",
    subtitle =
      "Modelo AR(1) estimado bajo una representación homogénea en KFAS",
    x = NULL,
    y = expression(beta),
    linetype = NULL,
    caption =
      "La banda sombreada corresponde al intervalo del 95% de beta_t; las líneas punteadas, al intervalo del beta constante."
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
    "12_beta_constante_vs_tvp.png"
  ),
  plot = grafica_beta_comparada,
  width = 10,
  height = 6,
  dpi = 300
)


# 19.2 Trayectoria TVP ----------------------------------------------------------

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
  labs(
    title =
      "Trayectoria de la sensibilidad variante de las expectativas",
    subtitle =
      paste0(
        "Bandas del ",
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
    "13_trayectoria_beta_tvp.png"
  ),
  plot = grafica_beta_tvp,
  width = 10,
  height = 6,
  dpi = 300
)


# 19.3 Residuos estandarizados comparados -------------------------------------

datos_residuos_largos <- resultados_beta_tvp %>%
  select(
    fecha,
    `Sensibilidad constante` =
      residuo_estandarizado_constante,
    `Sensibilidad variante` =
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
  labs(
    title =
      "Residuos recursivos estandarizados",
    subtitle =
      "Comparación entre sensibilidad constante y variante",
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
  labs(
    title =
      "Cambio observado y predicciones recursivas de un paso",
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
    "16_acf_residuos_constante.png"
  ),
  width = 1800,
  height = 1100,
  res = 180
)
acf(
  res_const_validos,
  lag.max = 24,
  main =
    "ACF de residuos: sensibilidad constante"
)
dev.off()

png(
  filename = file.path(
    carpeta_salida,
    "17_acf_residuos_tvp.png"
  ),
  width = 1800,
  height = 1100,
  res = 180
)
acf(
  res_tvp_validos,
  lag.max = 24,
  main =
    "ACF de residuos: sensibilidad variante"
)
dev.off()

png(
  filename = file.path(
    carpeta_salida,
    "18_pacf_residuos_constante.png"
  ),
  width = 1800,
  height = 1100,
  res = 180
)
pacf(
  res_const_validos,
  lag.max = 24,
  main =
    "PACF de residuos: sensibilidad constante"
)
dev.off()

png(
  filename = file.path(
    carpeta_salida,
    "19_pacf_residuos_tvp.png"
  ),
  width = 1800,
  height = 1100,
  res = 180
)
pacf(
  res_tvp_validos,
  lag.max = 24,
  main =
    "PACF de residuos: sensibilidad variante"
)
dev.off()


# 20. GUARDAR OBJETOS PARA REPRODUCIBILIDAD -----------------------------------

saveRDS(
  ajuste_constante$modelo,
  file.path(
    carpeta_salida,
    "20_modelo_constante_kfas.rds"
  )
)

saveRDS(
  ajuste_tvp$modelo,
  file.path(
    carpeta_salida,
    "21_modelo_tvp_kfas.rds"
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
  "COMPARACIÓN: SENSIBILIDAD CONSTANTE VS. VARIANTE EN EL TIEMPO",
  "======================================================================",
  paste0("Archivo: ", archivo_entrada),
  paste0(
    "Muestra: ",
    format(min(bd_modelo$fecha), "%Y-%m"),
    " a ",
    format(max(bd_modelo$fecha), "%Y-%m")
  ),
  paste0("Observaciones: ", n_obs),
  "",
  "MODELO CON SENSIBILIDAD CONSTANTE",
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
  "MODELO CON SENSIBILIDAD VARIANTE",
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

modelo_constante <- readRDS(
  "resultados_objetivo2_constante_vs_tvp_ar1/20_modelo_constante_kfas.rds"
)

modelo_tvp <- readRDS(
  "resultados_objetivo2_constante_vs_tvp_ar1/21_modelo_tvp_kfas.rds"
)
logLik(modelo_constante)
logLik(modelo_tvp)

resultado_tvp <- KFAS::KFS(
  modelo_tvp,
  filtering = c("state", "mean"),
  smoothing = c("state", "mean", "disturbance")
)
