# ==============================================================================
# OBJETIVO 2 - SCRIPT 10
# DEMERTZIS: PRUEBA DEL TVP RESTRINGIDO PARA CREDIBILIDAD EN GUATEMALA
# ==============================================================================
#
# PROPÓSITO
#
# Evaluar si la evidencia a favor del TVP-VAR de Demertzis proviene realmente
# de la ecuación de expectativas —que determina la proxy de credibilidad— o de
# la flexibilización de la ecuación de inflación.
#
# SISTEMA VAR(1)
#
#   Inflación:
#     pi_t = a0_t + a_t pi_(t-1) + b_t e_(t-1) + eps_pi,t
#
#   Expectativas a 24 meses:
#     e_t  = c0_t + c_t pi_(t-1) + d_t e_(t-1) + eps_e,t
#
# PROXY DE CREDIBILIDAD
#
#     lambda_t = 1 - c_t / (1 - d_t)
#
# ANCLA INFLACIONARIA IMPLÍCITA
#
#     pi_estrella_t = c0_t / (1 - d_t - c_t)
#
# MODELOS COMPARADOS
#
# M0. Todos los coeficientes constantes.
# M1. TVP solo en la ecuación de inflación: a0_t, a_t y b_t.
# M2. TVP completa en expectativas: c0_t, c_t y d_t;
#     ecuación de inflación constante.
# M3. CANDIDATO RESTRINGIDO:
#     c0_t y c_t variantes; d constante; inflación constante.
# M4. TVP-VAR irrestricto: los seis coeficientes variantes.
#
# DECISIÓN EMPÍRICA
#
# 1. M3 vs M0 determina si existe variación temporal en los parámetros que
#    construyen la proxy, sin permitir TVP en inflación.
# 2. M2 vs M3 determina si d_t necesita variar.
# 3. M4 vs M2 determina si agregar TVP en inflación aporta información una vez
#    que la ecuación de expectativas ya es TVP.
# 4. M1 vs M0 muestra cuánto de la mejora del TVP irrestricto puede provenir
#    exclusivamente de la ecuación de inflación.
#
# ADAPTACIÓN A GUATEMALA
#
# - Inflación observada: infl_gt.
# - Expectativa de mayor plazo disponible: exp_inf_24m.
# - Frecuencia mensual.
# - Muestra hasta 2026-06.
# - La meta explícita no entra en la estimación; se usa para contrastar el ancla
#   implícita estimada.
#
# PRECAUCIÓN
#
# La expectativa a 24 meses no equivale a una expectativa de largo plazo de
# seis a diez años. La interpretación es credibilidad/anclaje al horizonte de
# 24 meses disponible para Guatemala.
# ==============================================================================


# 0. LIMPIEZA Y PAQUETES -------------------------------------------------------

rm(list = ls())
cat("\014")
graphics.off()

paquetes <- c(
  "dplyr",
  "tidyr",
  "readr",
  "lubridate",
  "ggplot2",
  "tibble",
  "KFAS"
)

faltantes <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(faltantes) > 0) {
  install.packages(faltantes, dependencies = TRUE)
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
      "No se encontró la base de datos.\n",
      "Rutas revisadas:\n",
      paste0("- ", archivos_candidatos, collapse = "\n"),
      "\nDirectorio actual: ", getwd()
    )
  )
}

carpeta_salida <- file.path(
  getwd(),
  "resultados_objetivo2_demertzis_tvp_restringido"
)

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

columnas_requeridas <- c(
  "fecha",
  "infl_gt",
  "exp_inf_24m"
)

fecha_fin_muestra <- as.Date("2026-06-01")

nivel_confianza <- 0.95
z_critico <- stats::qnorm(1 - (1 - nivel_confianza) / 2)

limite_varianza <- 1e-12
tolerancia_denominador <- 0.05
rezagos_diagnostico <- c(1, 2, 5, 12, 18, 24)

semilla <- 12345

# Bootstrap paramétrico opcional para M3 vs M0.
# Puede tardar varios minutos. Primero conviene ejecutar con FALSE.
ejecutar_bootstrap_LR <- FALSE
replicas_bootstrap_LR <- 199


# 2. FUNCIONES GENERALES -------------------------------------------------------

verificar_columnas <- function(datos, columnas) {

  faltan <- setdiff(columnas, names(datos))

  if (length(faltan) > 0) {
    stop(
      "Faltan las siguientes variables requeridas: ",
      paste(faltan, collapse = ", ")
    )
  }
}


formato_num <- function(x, digitos = 5) {

  if (length(x) == 0 || !is.finite(x[[1]])) {
    return("NA")
  }

  formatC(
    x[[1]],
    digits = digitos,
    format = "f"
  )
}


formato_num_vector <- function(x, digitos = 5) {

  vapply(
    x,
    function(z) {
      if (!is.finite(z)) {
        "NA"
      } else {
        formatC(z, digits = digitos, format = "f")
      }
    },
    character(1)
  )
}


formato_cientifico_vector <- function(x, digitos = 4) {

  vapply(
    x,
    function(z) {
      if (!is.finite(z)) {
        "NA"
      } else {
        formatC(z, digits = digitos, format = "e")
      }
    },
    character(1)
  )
}


calcular_AICc <- function(AIC, k, n) {

  if (n <= k + 1) {
    return(NA_real_)
  }

  AIC + (2 * k * (k + 1)) / (n - k - 1)
}


jarque_bera_manual <- function(residuos) {

  residuos <- residuos[is.finite(residuos)]
  n <- length(residuos)

  if (n < 8 || !is.finite(stats::sd(residuos)) || stats::sd(residuos) == 0) {
    return(
      tibble::tibble(
        estadistico = NA_real_,
        p_valor = NA_real_
      )
    )
  }

  z <- (residuos - mean(residuos)) / stats::sd(residuos)
  asimetria <- mean(z^3)
  curtosis <- mean(z^4)

  JB <- n / 6 * (
    asimetria^2 +
      (curtosis - 3)^2 / 4
  )

  tibble::tibble(
    estadistico = JB,
    p_valor = stats::pchisq(JB, df = 2, lower.tail = FALSE)
  )
}


correlacion_segura <- function(x, y) {

  validos <- is.finite(x) & is.finite(y)

  if (sum(validos) < 4) {
    return(NA_real_)
  }

  x <- x[validos]
  y <- y[validos]

  if (
    !is.finite(stats::sd(x)) ||
      !is.finite(stats::sd(y)) ||
      stats::sd(x) == 0 ||
      stats::sd(y) == 0
  ) {
    return(NA_real_)
  }

  stats::cor(x, y)
}


# 3. CARGA Y PREPARACIÓN DE LOS DATOS -----------------------------------------

bd_original <- readr::read_csv(
  archivo_entrada,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", ".", "null")
)

verificar_columnas(
  bd_original,
  columnas_requeridas
)

bd <- bd_original |>
  dplyr::mutate(
    fecha = lubridate::ymd(fecha)
  ) |>
  dplyr::arrange(fecha) |>
  dplyr::filter(
    fecha <= fecha_fin_muestra
  )

if (anyNA(bd$fecha)) {
  stop("Existen fechas que no pudieron convertirse al formato Date.")
}

if (anyDuplicated(bd$fecha) > 0) {
  stop("Existen fechas duplicadas.")
}

if (!"meta_inflacion" %in% names(bd)) {

  bd <- bd |>
    dplyr::mutate(
      meta_inflacion = dplyr::case_when(
        fecha < lubridate::ymd("2012-01-01") ~ 5.0,
        fecha < lubridate::ymd("2013-01-01") ~ 4.5,
        TRUE ~ 4.0
      )
    )
}

datos_niveles <- bd |>
  dplyr::transmute(
    fecha,
    inflacion = infl_gt,
    expectativa = exp_inf_24m,
    meta_inflacion
  ) |>
  dplyr::filter(
    is.finite(inflacion),
    is.finite(expectativa),
    is.finite(meta_inflacion)
  )

if (nrow(datos_niveles) < 80) {
  stop("La muestra efectiva tiene menos de 80 observaciones.")
}

saltos_meses <- diff(
  lubridate::year(datos_niveles$fecha) * 12 +
    lubridate::month(datos_niveles$fecha)
)

if (any(saltos_meses != 1)) {
  stop(
    "La muestra presenta saltos mensuales. ",
    "El modelo requiere una secuencia temporal regular."
  )
}

# Datos utilizados por el VAR(1): se pierde la primera observación por el rezago.
datos_modelo <- datos_niveles |>
  dplyr::mutate(
    inflacion_lag1 = dplyr::lag(inflacion, 1),
    expectativa_lag1 = dplyr::lag(expectativa, 1)
  ) |>
  dplyr::filter(
    is.finite(inflacion_lag1),
    is.finite(expectativa_lag1)
  )

Y_modelo <- as.matrix(
  datos_modelo |>
    dplyr::select(
      inflacion,
      expectativa
    )
)

n_obs <- nrow(datos_modelo)

readr::write_csv(
  datos_modelo,
  file.path(
    carpeta_salida,
    "01_datos_utilizados_demertzis_restringido.csv"
  )
)


# 4. VAR(1) CONSTANTE PARA INICIALIZACIÓN -------------------------------------

inicializar_desde_OLS <- function(datos) {

  eq_inflacion <- stats::lm(
    inflacion ~ inflacion_lag1 + expectativa_lag1,
    data = datos,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  eq_expectativa <- stats::lm(
    expectativa ~ inflacion_lag1 + expectativa_lag1,
    data = datos,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  coef_inflacion <- stats::coef(eq_inflacion)
  coef_expectativa <- stats::coef(eq_expectativa)

  estado_inicial <- c(
    a0 = unname(coef_inflacion["(Intercept)"]),
    a = unname(coef_inflacion["inflacion_lag1"]),
    b = unname(coef_inflacion["expectativa_lag1"]),
    c0 = unname(coef_expectativa["(Intercept)"]),
    c = unname(coef_expectativa["inflacion_lag1"]),
    d = unname(coef_expectativa["expectativa_lag1"])
  )

  if (any(!is.finite(estado_inicial))) {
    stop(
      "La inicialización OLS produjo coeficientes no finitos. ",
      "Revise los datos y la colinealidad."
    )
  }

  H_inicial <- pmax(
    c(
      stats::var(stats::residuals(eq_inflacion), na.rm = TRUE),
      stats::var(stats::residuals(eq_expectativa), na.rm = TRUE)
    ),
    1e-8
  )

  list(
    estado_inicial = estado_inicial,
    H_inicial = H_inicial,
    eq_inflacion = eq_inflacion,
    eq_expectativa = eq_expectativa
  )
}

inicializacion <- inicializar_desde_OLS(datos_modelo)

estado_inicial <- inicializacion$estado_inicial
H_inicial <- inicializacion$H_inicial

coeficientes_OLS <- dplyr::bind_rows(
  tibble::tibble(
    ecuacion = "Inflación",
    parametro = names(stats::coef(inicializacion$eq_inflacion)),
    estimacion = as.numeric(stats::coef(inicializacion$eq_inflacion))
  ),
  tibble::tibble(
    ecuacion = "Expectativa 24m",
    parametro = names(stats::coef(inicializacion$eq_expectativa)),
    estimacion = as.numeric(stats::coef(inicializacion$eq_expectativa))
  )
)

readr::write_csv(
  coeficientes_OLS,
  file.path(
    carpeta_salida,
    "02_coeficientes_VAR1_OLS_inicializacion.csv"
  )
)


# 5. ESPECIFICACIONES DE LOS MODELOS ------------------------------------------

# Índices de los estados:
# 1 a0, 2 a, 3 b, 4 c0, 5 c, 6 d.
especificaciones <- list(
  M0 = list(
    nombre = "M0: VAR(1) constante SSM",
    q_activos = integer(0),
    descripcion = "Todos los coeficientes constantes"
  ),
  M1 = list(
    nombre = "M1: TVP solo inflación",
    q_activos = c(1L, 2L, 3L),
    descripcion = "a0_t, a_t y b_t variantes; expectativas constantes"
  ),
  M2 = list(
    nombre = "M2: TVP completa en expectativas",
    q_activos = c(4L, 5L, 6L),
    descripcion = "c0_t, c_t y d_t variantes; inflación constante"
  ),
  M3 = list(
    nombre = "M3: Demertzis restringido candidato",
    q_activos = c(4L, 5L),
    descripcion = "c0_t y c_t variantes; d e inflación constantes"
  ),
  M4 = list(
    nombre = "M4: TVP-VAR irrestricto",
    q_activos = 1:6,
    descripcion = "Los seis coeficientes variantes"
  )
)


# 6. CONSTRUCCIÓN GENÉRICA DEL MODELO SSM ------------------------------------

crear_modelo_demertzis <- function(
    datos,
    Y,
    estado_inicial,
    q_activos
) {

  n <- nrow(datos)
  m <- 6L

  Z_t <- array(
    0,
    dim = c(2, m, n)
  )

  # Ecuación de inflación.
  Z_t[1, 1, ] <- 1
  Z_t[1, 2, ] <- datos$inflacion_lag1
  Z_t[1, 3, ] <- datos$expectativa_lag1

  # Ecuación de expectativas.
  Z_t[2, 4, ] <- 1
  Z_t[2, 5, ] <- datos$inflacion_lag1
  Z_t[2, 6, ] <- datos$expectativa_lag1

  Q_mat <- matrix(
    0,
    nrow = m,
    ncol = m
  )

  if (length(q_activos) > 0) {
    for (j in q_activos) {
      Q_mat[j, j] <- NA_real_
    }
  }

  KFAS::SSModel(
    Y ~ -1 +
      SSMcustom(
        Z = Z_t,
        T = diag(m),
        R = diag(m),
        Q = Q_mat,
        a1 = estado_inicial,
        P1 = matrix(0, nrow = m, ncol = m),
        P1inf = diag(m),
        n = n
      ),
    H = diag(rep(NA_real_, 2))
  )
}


crear_updatefn <- function(q_activos) {

  force(q_activos)

  function(par, modelo) {

    modelo$H[1, 1, 1] <- exp(par[1])
    modelo$H[2, 2, 1] <- exp(par[2])
    modelo$H[1, 2, 1] <- 0
    modelo$H[2, 1, 1] <- 0

    if (length(q_activos) > 0) {

      # Mantener en cero todas las varianzas restringidas.
      for (j in 1:6) {
        modelo$Q[j, j, 1] <- 0
      }

      valores_Q <- exp(par[2 + seq_along(q_activos)])

      for (j in seq_along(q_activos)) {
        pos <- q_activos[j]
        modelo$Q[pos, pos, 1] <- valores_Q[j]
      }
    }

    modelo
  }
}


# 7. MÚLTIPLES VALORES INICIALES ----------------------------------------------

crear_inicios <- function(H_inicial, numero_Q) {

  if (numero_Q == 0) {

    return(
      lapply(
        c(0.5, 1, 2),
        function(f) {
          log(pmax(H_inicial * f, limite_varianza))
        }
      )
    )
  }

  escalas_Q <- c(
    1e-10,
    1e-8,
    1e-7,
    1e-6,
    1e-5,
    1e-4,
    1e-3,
    1e-2
  )

  inicios <- lapply(
    escalas_Q,
    function(q) {
      log(
        c(
          H_inicial,
          rep(q, numero_Q)
        )
      )
    }
  )

  # Cambios en las varianzas de observación.
  inicios <- c(
    inicios,
    list(
      log(c(H_inicial * 0.5, rep(1e-5, numero_Q))),
      log(c(H_inicial * 2.0, rep(1e-5, numero_Q)))
    )
  )

  # Inicios asimétricos para reducir el riesgo de soluciones locales.
  if (numero_Q >= 2) {

    for (j in seq_len(numero_Q)) {
      q_asimetrico <- rep(1e-6, numero_Q)
      q_asimetrico[j] <- 1e-3

      inicios <- c(
        inicios,
        list(log(c(H_inicial, q_asimetrico)))
      )
    }
  }

  inicios
}


estimar_especificacion <- function(
    codigo,
    especificacion,
    datos,
    Y,
    estado_inicial,
    H_inicial,
    maxit = 30000
) {

  q_activos <- especificacion$q_activos

  modelo_sin_estimar <- crear_modelo_demertzis(
    datos = datos,
    Y = Y,
    estado_inicial = estado_inicial,
    q_activos = q_activos
  )

  updatefn <- crear_updatefn(q_activos)
  inicios <- crear_inicios(H_inicial, length(q_activos))

  resultados <- lapply(
    seq_along(inicios),
    function(i) {

      mensaje_error <- NA_character_

      ajuste <- tryCatch(
        KFAS::fitSSM(
          model = modelo_sin_estimar,
          inits = inicios[[i]],
          updatefn = updatefn,
          method = "L-BFGS-B",
          lower = rep(
            log(limite_varianza),
            2 + length(q_activos)
          ),
          upper = rep(
            log(100),
            2 + length(q_activos)
          ),
          control = list(
            maxit = maxit,
            factr = 1e7
          )
        ),
        error = function(e) {
          mensaje_error <<- conditionMessage(e)
          NULL
        }
      )

      if (is.null(ajuste)) {
        return(
          list(
            exito = FALSE,
            intento = i,
            convergencia = NA_integer_,
            logLik = NA_real_,
            mensaje = mensaje_error,
            ajuste = NULL
          )
        )
      }

      ll <- tryCatch(
        as.numeric(stats::logLik(ajuste$model)),
        error = function(e) NA_real_
      )

      list(
        exito = is.finite(ll),
        intento = i,
        convergencia = ajuste$optim.out$convergence,
        logLik = ll,
        mensaje = ifelse(
          is.null(ajuste$optim.out$message),
          "",
          as.character(ajuste$optim.out$message)
        ),
        ajuste = ajuste
      )
    }
  )

  tabla_intentos <- dplyr::bind_rows(
    lapply(
      resultados,
      function(x) {
        tibble::tibble(
          codigo_modelo = codigo,
          modelo = especificacion$nombre,
          intento = x$intento,
          exito = x$exito,
          convergencia = x$convergencia,
          log_verosimilitud = x$logLik,
          mensaje = x$mensaje
        )
      }
    )
  )

  exitosos <- resultados[
    vapply(
      resultados,
      function(x) {
        isTRUE(x$exito) && is.finite(x$logLik)
      },
      logical(1)
    )
  ]

  if (length(exitosos) == 0) {

    readr::write_csv(
      tabla_intentos,
      file.path(
        carpeta_salida,
        paste0("ERROR_intentos_", codigo, ".csv")
      )
    )

    stop(
      "No fue posible estimar ", especificacion$nombre,
      ". Revise el archivo de intentos guardado."
    )
  }

  convergentes <- exitosos[
    vapply(
      exitosos,
      function(x) x$convergencia == 0,
      logical(1)
    )
  ]

  elegibles <- if (length(convergentes) > 0) {
    convergentes
  } else {
    warning(
      "Ningún intento convergió normalmente para ",
      especificacion$nombre,
      ". Se utilizará el de mayor log-verosimilitud."
    )
    exitosos
  }

  mejor <- elegibles[[
    which.max(
      vapply(elegibles, function(x) x$logLik, numeric(1))
    )
  ]]

  list(
    codigo = codigo,
    especificacion = especificacion,
    mejor = mejor,
    modelo = mejor$ajuste$model,
    tabla_intentos = tabla_intentos
  )
}


# 8. ESTIMACIÓN DE LOS CINCO MODELOS ------------------------------------------

ajustes <- list()

for (codigo in names(especificaciones)) {

  cat("\nEstimando ", especificaciones[[codigo]]$nombre, "...\n", sep = "")

  ajustes[[codigo]] <- estimar_especificacion(
    codigo = codigo,
    especificacion = especificaciones[[codigo]],
    datos = datos_modelo,
    Y = Y_modelo,
    estado_inicial = estado_inicial,
    H_inicial = H_inicial
  )
}

tabla_intentos <- dplyr::bind_rows(
  lapply(ajustes, function(x) x$tabla_intentos)
)

readr::write_csv(
  tabla_intentos,
  file.path(
    carpeta_salida,
    "03_intentos_estimacion_cinco_modelos.csv"
  )
)


# 9. FILTRO Y SUAVIZADOR DE KALMAN --------------------------------------------

resultados_kfs <- lapply(
  ajustes,
  function(x) {
    KFAS::KFS(
      x$modelo,
      filtering = c("state", "mean"),
      smoothing = c("state", "mean", "disturbance")
    )
  }
)


# 10. HIPERPARÁMETROS Y CRITERIOS DE INFORMACIÓN ------------------------------

nombres_Q <- c(
  "Q_a0",
  "Q_a",
  "Q_b",
  "Q_c0",
  "Q_c",
  "Q_d"
)

extraer_hiperparametros <- function(codigo, ajuste) {

  H <- diag(ajuste$modelo$H[, , 1])
  Q <- diag(ajuste$modelo$Q[, , 1])

  fila <- tibble::tibble(
    codigo_modelo = codigo,
    modelo = ajuste$especificacion$nombre,
    descripcion = ajuste$especificacion$descripcion,
    H_inflacion = H[1],
    H_expectativa = H[2],
    convergencia = ajuste$mejor$convergencia,
    log_verosimilitud = ajuste$mejor$logLik
  )

  for (j in seq_along(nombres_Q)) {
    fila[[nombres_Q[j]]] <- Q[j]
  }

  fila
}

tabla_hiperparametros <- dplyr::bind_rows(
  lapply(
    names(ajustes),
    function(codigo) {
      extraer_hiperparametros(codigo, ajustes[[codigo]])
    }
  )
)

calcular_fila_comparacion <- function(codigo, ajuste) {

  ll <- as.numeric(stats::logLik(ajuste$modelo))
  numero_Q <- length(ajuste$especificacion$q_activos)

  # Se contabilizan seis coeficientes del VAR, dos varianzas H y las Q activas.
  k <- 8 + numero_Q

  AIC <- -2 * ll + 2 * k
  BIC <- -2 * ll + log(n_obs) * k
  AICc <- calcular_AICc(AIC, k, n_obs)

  tibble::tibble(
    codigo_modelo = codigo,
    modelo = ajuste$especificacion$nombre,
    descripcion = ajuste$especificacion$descripcion,
    observaciones = n_obs,
    parametros_contados = k,
    log_verosimilitud = ll,
    AIC = AIC,
    AICc = AICc,
    BIC = BIC
  )
}

tabla_comparacion <- dplyr::bind_rows(
  lapply(
    names(ajustes),
    function(codigo) {
      calcular_fila_comparacion(codigo, ajustes[[codigo]])
    }
  )
) |>
  dplyr::mutate(
    delta_AIC = AIC - min(AIC),
    peso_Akaike = exp(-0.5 * delta_AIC) /
      sum(exp(-0.5 * delta_AIC)),
    delta_BIC = BIC - min(BIC),
    peso_BIC_aproximado = exp(-0.5 * delta_BIC) /
      sum(exp(-0.5 * delta_BIC))
  ) |>
  dplyr::arrange(AIC)

readr::write_csv(
  tabla_hiperparametros,
  file.path(
    carpeta_salida,
    "04_hiperparametros_cinco_modelos.csv"
  )
)

readr::write_csv(
  tabla_comparacion,
  file.path(
    carpeta_salida,
    "05_comparacion_cinco_modelos.csv"
  )
)


# 11. PRUEBAS LR ENTRE MODELOS ANIDADOS ---------------------------------------

obtener_ll <- function(codigo) {
  as.numeric(stats::logLik(ajustes[[codigo]]$modelo))
}

crear_prueba_LR <- function(
    restringido,
    irrestricto,
    restricciones,
    descripcion
) {

  ll_r <- obtener_ll(restringido)
  ll_u <- obtener_ll(irrestricto)

  LR <- max(0, 2 * (ll_u - ll_r))

  p_chi <- stats::pchisq(
    LR,
    df = restricciones,
    lower.tail = FALSE
  )

  p_mixto_una_varianza <- if (restricciones == 1) {
    0.5 * stats::pchisq(LR, df = 1, lower.tail = FALSE)
  } else {
    NA_real_
  }

  tibble::tibble(
    modelo_restringido = restringido,
    modelo_irrestricto = irrestricto,
    descripcion = descripcion,
    logLik_restringido = ll_r,
    logLik_irrestricto = ll_u,
    estadistico_LR = LR,
    numero_varianzas_en_frontera = restricciones,
    p_valor_chi_cuadrado_referencial = p_chi,
    p_valor_mixto_una_varianza = p_mixto_una_varianza,
    advertencia = paste(
      "Las varianzas Q son cero bajo H0 y se encuentran en la frontera.",
      "Los p-valores son referenciales; deben evaluarse junto con AIC, BIC,",
      "magnitud de Q y, para M3 vs M0, el bootstrap paramétrico opcional."
    )
  )
}

tabla_pruebas_LR <- dplyr::bind_rows(
  crear_prueba_LR(
    "M0", "M1", 3,
    "¿La TVP exclusivamente en inflación mejora al modelo constante?"
  ),
  crear_prueba_LR(
    "M0", "M3", 2,
    "Prueba central: ¿c0_t y c_t variantes mejoran al modelo constante?"
  ),
  crear_prueba_LR(
    "M3", "M2", 1,
    "¿Es necesario permitir que d_t varíe adicionalmente?"
  ),
  crear_prueba_LR(
    "M2", "M4", 3,
    "¿La TVP en inflación añade información después de hacer TVP la ecuación de expectativas?"
  ),
  crear_prueba_LR(
    "M3", "M4", 4,
    "¿El TVP irrestricto mejora al candidato restringido?"
  ),
  crear_prueba_LR(
    "M0", "M4", 6,
    "Prueba global referencial: TVP-VAR irrestricto frente a constante"
  )
)

readr::write_csv(
  tabla_pruebas_LR,
  file.path(
    carpeta_salida,
    "06_pruebas_LR_modelos_anidados.csv"
  )
)


# 12. EXTRACCIÓN DE ESTADOS Y PROXY -------------------------------------------

nombres_estados <- c(
  "a0_t",
  "a_t",
  "b_t",
  "c0_t",
  "c_t",
  "d_t"
)

calcular_proxy_modelo <- function(codigo, kfs) {

  estados <- as.data.frame(kfs$alphahat)
  names(estados) <- nombres_estados

  trayectoria <- dplyr::bind_cols(
    datos_modelo |>
      dplyr::select(
        fecha,
        inflacion,
        expectativa,
        meta_inflacion
      ),
    estados
  )

  lista_proxy <- lapply(
    seq_len(nrow(trayectoria)),
    function(t) {

      c0 <- trayectoria$c0_t[t]
      c <- trayectoria$c_t[t]
      d <- trayectoria$d_t[t]

      den_d <- 1 - d
      den_pi <- 1 - d - c

      if (
        !is.finite(den_d) ||
          !is.finite(den_pi) ||
          abs(den_d) < tolerancia_denominador ||
          abs(den_pi) < tolerancia_denominador
      ) {
        return(
          tibble::tibble(
            lambda_t = NA_real_,
            lambda_se = NA_real_,
            pi_estrella_t = NA_real_,
            pi_estrella_se = NA_real_,
            denominador_1_menos_d = den_d,
            denominador_ancla = den_pi
          )
        )
      }

      V_sub <- kfs$V[
        c(4, 5, 6),
        c(4, 5, 6),
        t,
        drop = FALSE
      ]

      V_sub <- matrix(V_sub, nrow = 3, ncol = 3)

      lambda <- 1 - c / den_d
      pi_estrella <- c0 / den_pi

      grad_lambda <- matrix(
        c(
          0,
          -1 / den_d,
          -c / den_d^2
        ),
        ncol = 1
      )

      grad_pi <- matrix(
        c(
          1 / den_pi,
          c0 / den_pi^2,
          c0 / den_pi^2
        ),
        ncol = 1
      )

      var_lambda <- tryCatch(
        as.numeric(t(grad_lambda) %*% V_sub %*% grad_lambda),
        error = function(e) NA_real_
      )

      var_pi <- tryCatch(
        as.numeric(t(grad_pi) %*% V_sub %*% grad_pi),
        error = function(e) NA_real_
      )

      tibble::tibble(
        lambda_t = lambda,
        lambda_se = ifelse(
          is.finite(var_lambda) && var_lambda >= 0,
          sqrt(var_lambda),
          NA_real_
        ),
        pi_estrella_t = pi_estrella,
        pi_estrella_se = ifelse(
          is.finite(var_pi) && var_pi >= 0,
          sqrt(var_pi),
          NA_real_
        ),
        denominador_1_menos_d = den_d,
        denominador_ancla = den_pi
      )
    }
  )

  proxy <- dplyr::bind_rows(lista_proxy)

  dplyr::bind_cols(trayectoria, proxy) |>
    dplyr::mutate(
      codigo_modelo = codigo,
      lambda_li_95 = lambda_t - z_critico * lambda_se,
      lambda_ls_95 = lambda_t + z_critico * lambda_se,
      pi_estrella_li_95 = pi_estrella_t - z_critico * pi_estrella_se,
      pi_estrella_ls_95 = pi_estrella_t + z_critico * pi_estrella_se,
      lambda_fuera_0_1 = lambda_t < 0 | lambda_t > 1,
      lambda_significativamente_menor_1 = lambda_ls_95 < 1,
      lambda_no_distinta_de_1 = lambda_li_95 <= 1 & lambda_ls_95 >= 1,
      brecha_ancla_meta = pi_estrella_t - meta_inflacion,
      ancla_dentro_banda_meta =
        pi_estrella_t >= meta_inflacion - 1 &
        pi_estrella_t <= meta_inflacion + 1
    )
}

trayectorias <- lapply(
  names(resultados_kfs),
  function(codigo) {
    calcular_proxy_modelo(codigo, resultados_kfs[[codigo]])
  }
)

names(trayectorias) <- names(resultados_kfs)

trayectoria_candidata <- trayectorias$M3

readr::write_csv(
  trayectoria_candidata,
  file.path(
    carpeta_salida,
    "07_trayectoria_candidato_M3_lambda_pi_estrella.csv"
  )
)

trayectorias_comparables <- dplyr::bind_rows(
  lapply(
    c("M0", "M2", "M3", "M4"),
    function(codigo) {
      trayectorias[[codigo]] |>
        dplyr::select(
          codigo_modelo,
          fecha,
          lambda_t,
          lambda_li_95,
          lambda_ls_95,
          pi_estrella_t,
          c0_t,
          c_t,
          d_t
        )
    }
  )
)

readr::write_csv(
  trayectorias_comparables,
  file.path(
    carpeta_salida,
    "08_trayectorias_proxy_modelos_comparables.csv"
  )
)


# 13. RESUMEN DEL CANDIDATO M3 -------------------------------------------------

indice_min_lambda <- which.min(trayectoria_candidata$lambda_t)
indice_max_lambda <- which.max(trayectoria_candidata$lambda_t)

resumen_candidato <- tibble::tibble(
  fecha_inicial = dplyr::first(trayectoria_candidata$fecha),
  fecha_final = dplyr::last(trayectoria_candidata$fecha),
  observaciones = nrow(trayectoria_candidata),
  lambda_promedio = mean(trayectoria_candidata$lambda_t, na.rm = TRUE),
  lambda_inicial = dplyr::first(trayectoria_candidata$lambda_t),
  lambda_final = dplyr::last(trayectoria_candidata$lambda_t),
  lambda_minimo = trayectoria_candidata$lambda_t[indice_min_lambda],
  fecha_lambda_minimo = trayectoria_candidata$fecha[indice_min_lambda],
  lambda_maximo = trayectoria_candidata$lambda_t[indice_max_lambda],
  fecha_lambda_maximo = trayectoria_candidata$fecha[indice_max_lambda],
  porcentaje_lambda_fuera_0_1 = 100 * mean(
    trayectoria_candidata$lambda_fuera_0_1,
    na.rm = TRUE
  ),
  porcentaje_lambda_significativamente_menor_1 = 100 * mean(
    trayectoria_candidata$lambda_significativamente_menor_1,
    na.rm = TRUE
  ),
  pi_estrella_promedio = mean(
    trayectoria_candidata$pi_estrella_t,
    na.rm = TRUE
  ),
  pi_estrella_inicial = dplyr::first(
    trayectoria_candidata$pi_estrella_t
  ),
  pi_estrella_final = dplyr::last(
    trayectoria_candidata$pi_estrella_t
  ),
  MAE_pi_estrella_meta = mean(
    abs(trayectoria_candidata$brecha_ancla_meta),
    na.rm = TRUE
  ),
  RMSE_pi_estrella_meta = sqrt(
    mean(
      trayectoria_candidata$brecha_ancla_meta^2,
      na.rm = TRUE
    )
  ),
  porcentaje_ancla_dentro_banda_meta = 100 * mean(
    trayectoria_candidata$ancla_dentro_banda_meta,
    na.rm = TRUE
  ),
  c0_promedio = mean(trayectoria_candidata$c0_t, na.rm = TRUE),
  c_promedio = mean(trayectoria_candidata$c_t, na.rm = TRUE),
  d_estimado = mean(trayectoria_candidata$d_t, na.rm = TRUE)
)

resumen_subperiodos <- trayectoria_candidata |>
  dplyr::mutate(
    subperiodo = dplyr::case_when(
      fecha <= lubridate::ymd("2014-12-01") ~ "2010-2014",
      fecha <= lubridate::ymd("2019-12-01") ~ "2015-2019",
      fecha <= lubridate::ymd("2022-12-01") ~ "2020-2022",
      TRUE ~ "2023-2026"
    )
  ) |>
  dplyr::group_by(subperiodo) |>
  dplyr::summarise(
    fecha_inicial = min(fecha),
    fecha_final = max(fecha),
    observaciones = dplyr::n(),
    lambda_promedio = mean(lambda_t, na.rm = TRUE),
    c0_promedio = mean(c0_t, na.rm = TRUE),
    c_promedio = mean(c_t, na.rm = TRUE),
    d_promedio = mean(d_t, na.rm = TRUE),
    pi_estrella_promedio = mean(pi_estrella_t, na.rm = TRUE),
    meta_promedio = mean(meta_inflacion, na.rm = TRUE),
    brecha_ancla_meta_promedio = mean(brecha_ancla_meta, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  resumen_candidato,
  file.path(
    carpeta_salida,
    "09_resumen_candidato_M3.csv"
  )
)

readr::write_csv(
  resumen_subperiodos,
  file.path(
    carpeta_salida,
    "10_resumen_candidato_por_subperiodos.csv"
  )
)


# 14. DIAGNÓSTICOS RESIDUALES --------------------------------------------------

extraer_residuos_estandarizados <- function(kfs, n_esperado) {

  residuos <- tryCatch(
    as.matrix(
      stats::rstandard(
        kfs,
        type = "recursive"
      )
    ),
    error = function(e) {
      matrix(
        NA_real_,
        nrow = n_esperado,
        ncol = 2
      )
    }
  )

  if (nrow(residuos) != n_esperado || ncol(residuos) != 2) {
    residuos <- matrix(
      NA_real_,
      nrow = n_esperado,
      ncol = 2
    )
  }

  colnames(residuos) <- c("inflacion", "expectativa")
  residuos
}

residuos_modelos <- lapply(
  resultados_kfs,
  extraer_residuos_estandarizados,
  n_esperado = n_obs
)

crear_diagnosticos_modelo <- function(codigo, residuos) {

  filas <- list()

  for (ecuacion in c("inflacion", "expectativa")) {

    r <- residuos[, ecuacion]
    r <- r[is.finite(r)]

    for (L in rezagos_diagnostico) {

      if (length(r) > L + 5) {

        prueba <- stats::Box.test(
          r,
          lag = L,
          type = "Ljung-Box",
          fitdf = 0
        )

        filas[[length(filas) + 1]] <- tibble::tibble(
          codigo_modelo = codigo,
          ecuacion = ecuacion,
          diagnostico = "Ljung-Box",
          rezago = L,
          estadistico = unname(prueba$statistic),
          p_valor = prueba$p.value
        )
      }
    }

    jb <- jarque_bera_manual(r)

    filas[[length(filas) + 1]] <- tibble::tibble(
      codigo_modelo = codigo,
      ecuacion = ecuacion,
      diagnostico = "Jarque-Bera",
      rezago = NA_real_,
      estadistico = jb$estadistico,
      p_valor = jb$p_valor
    )
  }

  rho <- correlacion_segura(
    residuos[, "inflacion"],
    residuos[, "expectativa"]
  )

  validos_sistema <- is.finite(residuos[, "inflacion"]) &
    is.finite(residuos[, "expectativa"])

  n_sistema <- sum(validos_sistema)

  fisher_z <- if (
    is.finite(rho) &&
      abs(rho) < 1 &&
      n_sistema > 3
  ) {
    atanh(rho) * sqrt(n_sistema - 3)
  } else {
    NA_real_
  }

  p_fisher <- if (is.finite(fisher_z)) {
    2 * stats::pnorm(abs(fisher_z), lower.tail = FALSE)
  } else {
    NA_real_
  }

  filas[[length(filas) + 1]] <- tibble::tibble(
    codigo_modelo = codigo,
    ecuacion = "sistema",
    diagnostico = "Correlación contemporánea Fisher",
    rezago = NA_real_,
    estadistico = fisher_z,
    p_valor = p_fisher
  )

  dplyr::bind_rows(filas) |>
    dplyr::mutate(
      conclusion_5pct = dplyr::case_when(
        !is.finite(p_valor) ~ "No disponible",
        p_valor < 0.05 ~ "Se rechaza la hipótesis nula",
        TRUE ~ "No se rechaza la hipótesis nula"
      )
    )
}

tabla_diagnosticos <- dplyr::bind_rows(
  lapply(
    names(residuos_modelos),
    function(codigo) {
      crear_diagnosticos_modelo(
        codigo,
        residuos_modelos[[codigo]]
      )
    }
  )
)

readr::write_csv(
  tabla_diagnosticos,
  file.path(
    carpeta_salida,
    "11_diagnosticos_residuales_cinco_modelos.csv"
  )
)

# Residuos extremos del candidato.
residuos_M3 <- residuos_modelos$M3

residuos_candidato <- tibble::tibble(
  fecha = datos_modelo$fecha,
  residuo_inflacion = residuos_M3[, "inflacion"],
  residuo_expectativa = residuos_M3[, "expectativa"]
)

extremos_candidato <- dplyr::bind_rows(
  residuos_candidato |>
    dplyr::filter(is.finite(residuo_inflacion)) |>
    dplyr::transmute(
      fecha,
      ecuacion = "inflacion",
      residuo_estandarizado = residuo_inflacion
    ),
  residuos_candidato |>
    dplyr::filter(is.finite(residuo_expectativa)) |>
    dplyr::transmute(
      fecha,
      ecuacion = "expectativa",
      residuo_estandarizado = residuo_expectativa
    )
) |>
  dplyr::arrange(dplyr::desc(abs(residuo_estandarizado))) |>
  dplyr::slice_head(n = 20)

readr::write_csv(
  extremos_candidato,
  file.path(
    carpeta_salida,
    "12_residuos_extremos_candidato_M3.csv"
  )
)


# 15. ESTABILIDAD DEL VAR CONSTANTE -------------------------------------------

estados_M0 <- as.data.frame(resultados_kfs$M0$alphahat)
names(estados_M0) <- nombres_estados

coef_M0 <- colMeans(estados_M0, na.rm = TRUE)

matriz_VAR_M0 <- matrix(
  c(
    coef_M0["a_t"], coef_M0["b_t"],
    coef_M0["c_t"], coef_M0["d_t"]
  ),
  nrow = 2,
  byrow = TRUE
)

raices_M0 <- eigen(matriz_VAR_M0, only.values = TRUE)$values
radio_espectral_M0 <- max(Mod(raices_M0))

estabilidad_M0 <- tibble::tibble(
  radio_espectral = radio_espectral_M0,
  VAR_estable = radio_espectral_M0 < 1,
  raiz_1_real = Re(raices_M0[1]),
  raiz_1_imaginaria = Im(raices_M0[1]),
  raiz_2_real = Re(raices_M0[2]),
  raiz_2_imaginaria = Im(raices_M0[2])
)

readr::write_csv(
  estabilidad_M0,
  file.path(
    carpeta_salida,
    "13_estabilidad_VAR1_constante.csv"
  )
)


# 16. BOOTSTRAP PARAMÉTRICO OPCIONAL: M3 VS M0 --------------------------------

resultado_bootstrap <- tibble::tibble(
  replicas_solicitadas = replicas_bootstrap_LR,
  replicas_validas = 0,
  LR_observado = tabla_pruebas_LR |>
    dplyr::filter(
      modelo_restringido == "M0",
      modelo_irrestricto == "M3"
    ) |>
    dplyr::pull(estadistico_LR),
  p_valor_bootstrap = NA_real_,
  ejecutado = ejecutar_bootstrap_LR,
  nota = "Bootstrap no ejecutado. Cambie ejecutar_bootstrap_LR a TRUE."
)

if (ejecutar_bootstrap_LR) {

  if (radio_espectral_M0 >= 1) {
    warning(
      "No se ejecutará el bootstrap: el VAR(1) constante estimado no es estable."
    )
  } else {

    set.seed(semilla)

    H_M0 <- diag(ajustes$M0$modelo$H[, , 1])

    a0 <- coef_M0["a0_t"]
    a <- coef_M0["a_t"]
    b <- coef_M0["b_t"]
    c0 <- coef_M0["c0_t"]
    c_coef <- coef_M0["c_t"]
    d_coef <- coef_M0["d_t"]

    valor_inicial <- c(
      inflacion = datos_niveles$inflacion[1],
      expectativa = datos_niveles$expectativa[1]
    )

    simular_VAR1_M0 <- function() {

      pi_prev <- valor_inicial["inflacion"]
      e_prev <- valor_inicial["expectativa"]

      pi_sim <- numeric(n_obs)
      e_sim <- numeric(n_obs)
      pi_lag <- numeric(n_obs)
      e_lag <- numeric(n_obs)

      for (t in seq_len(n_obs)) {

        pi_lag[t] <- pi_prev
        e_lag[t] <- e_prev

        innov_pi <- stats::rnorm(1, 0, sqrt(H_M0[1]))
        innov_e <- stats::rnorm(1, 0, sqrt(H_M0[2]))

        pi_actual <- a0 + a * pi_prev + b * e_prev + innov_pi
        e_actual <- c0 + c_coef * pi_prev + d_coef * e_prev + innov_e

        pi_sim[t] <- pi_actual
        e_sim[t] <- e_actual

        pi_prev <- pi_actual
        e_prev <- e_actual
      }

      tibble::tibble(
        fecha = datos_modelo$fecha,
        inflacion = pi_sim,
        expectativa = e_sim,
        meta_inflacion = datos_modelo$meta_inflacion,
        inflacion_lag1 = pi_lag,
        expectativa_lag1 = e_lag
      )
    }

    LR_boot <- rep(NA_real_, replicas_bootstrap_LR)

    for (r in seq_len(replicas_bootstrap_LR)) {

      datos_sim <- simular_VAR1_M0()
      Y_sim <- as.matrix(datos_sim[, c("inflacion", "expectativa")])

      init_sim <- tryCatch(
        inicializar_desde_OLS(datos_sim),
        error = function(e) NULL
      )

      if (is.null(init_sim)) {
        next
      }

      ajuste_sim_M0 <- tryCatch(
        estimar_especificacion(
          codigo = "M0_boot",
          especificacion = especificaciones$M0,
          datos = datos_sim,
          Y = Y_sim,
          estado_inicial = init_sim$estado_inicial,
          H_inicial = init_sim$H_inicial,
          maxit = 10000
        ),
        error = function(e) NULL
      )

      ajuste_sim_M3 <- tryCatch(
        estimar_especificacion(
          codigo = "M3_boot",
          especificacion = especificaciones$M3,
          datos = datos_sim,
          Y = Y_sim,
          estado_inicial = init_sim$estado_inicial,
          H_inicial = init_sim$H_inicial,
          maxit = 15000
        ),
        error = function(e) NULL
      )

      if (is.null(ajuste_sim_M0) || is.null(ajuste_sim_M3)) {
        next
      }

      ll0 <- as.numeric(stats::logLik(ajuste_sim_M0$modelo))
      ll3 <- as.numeric(stats::logLik(ajuste_sim_M3$modelo))

      LR_boot[r] <- max(0, 2 * (ll3 - ll0))

      if (r %% 10 == 0) {
        cat("Bootstrap LR: réplica ", r, " de ", replicas_bootstrap_LR, "\n", sep = "")
      }
    }

    LR_validos <- LR_boot[is.finite(LR_boot)]
    LR_observado <- resultado_bootstrap$LR_observado[1]

    p_boot <- if (length(LR_validos) > 0) {
      (1 + sum(LR_validos >= LR_observado)) /
        (1 + length(LR_validos))
    } else {
      NA_real_
    }

    resultado_bootstrap <- tibble::tibble(
      replicas_solicitadas = replicas_bootstrap_LR,
      replicas_validas = length(LR_validos),
      LR_observado = LR_observado,
      p_valor_bootstrap = p_boot,
      ejecutado = TRUE,
      nota = ifelse(
        length(LR_validos) >= 0.8 * replicas_bootstrap_LR,
        "Número de réplicas válidas adecuado.",
        "Menos de 80% de réplicas válidas; interpretar con cautela."
      )
    )

    readr::write_csv(
      tibble::tibble(LR_bootstrap = LR_validos),
      file.path(
        carpeta_salida,
        "14_distribucion_bootstrap_LR_M3_vs_M0.csv"
      )
    )
  }
}

readr::write_csv(
  resultado_bootstrap,
  file.path(
    carpeta_salida,
    "15_resultado_bootstrap_LR_M3_vs_M0.csv"
  )
)


# 17. GRÁFICAS -----------------------------------------------------------------

# 17.1 Proxy lambda del candidato restringido.
p_lambda_M3 <- ggplot2::ggplot(
  trayectoria_candidata,
  ggplot2::aes(x = fecha, y = lambda_t)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = lambda_li_95,
      ymax = lambda_ls_95
    ),
    alpha = 0.18
  ) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
  ggplot2::labs(
    title = "Proxy de credibilidad: Demertzis restringido",
    subtitle = expression(lambda[t] == 1 - c[t] / (1 - d)),
    x = NULL,
    y = expression(lambda[t]),
    caption = paste(
      "c0_t y c_t son variantes; d y los parámetros de inflación son constantes.",
      "Bandas del 95%."
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = file.path(
    carpeta_salida,
    "16_lambda_t_candidato_M3.png"
  ),
  plot = p_lambda_M3,
  width = 10,
  height = 6,
  dpi = 300
)

# 17.2 Comparación de lambda entre modelos de expectativas.
datos_lambda_comparacion <- trayectorias_comparables |>
  dplyr::filter(codigo_modelo %in% c("M0", "M2", "M3", "M4"))

p_lambda_comparacion <- ggplot2::ggplot(
  datos_lambda_comparacion,
  ggplot2::aes(
    x = fecha,
    y = lambda_t,
    linetype = codigo_modelo
  )
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dotted") +
  ggplot2::labs(
    title = "Comparación de la proxy de credibilidad",
    subtitle = "M0 constante, M2 expectativas TVP, M3 restringido y M4 irrestricto",
    x = NULL,
    y = expression(lambda[t]),
    linetype = "Modelo"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = file.path(
    carpeta_salida,
    "17_comparacion_lambda_modelos.png"
  ),
  plot = p_lambda_comparacion,
  width = 10,
  height = 6,
  dpi = 300
)

# 17.3 Ancla implícita del candidato frente a meta y expectativas.
datos_ancla <- trayectoria_candidata |>
  dplyr::select(
    fecha,
    `Ancla implícita M3` = pi_estrella_t,
    `Expectativa 24 meses` = expectativa,
    Meta = meta_inflacion
  ) |>
  tidyr::pivot_longer(
    cols = -fecha,
    names_to = "serie",
    values_to = "valor"
  )

p_ancla <- ggplot2::ggplot(
  datos_ancla,
  ggplot2::aes(
    x = fecha,
    y = valor,
    linetype = serie
  )
) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::labs(
    title = "Ancla inflacionaria implícita del modelo restringido",
    subtitle = expression(pi[t]^"*" == c[0*t] / (1 - d - c[t])),
    x = NULL,
    y = "Porcentaje",
    linetype = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = file.path(
    carpeta_salida,
    "18_pi_estrella_M3_expectativa_meta.png"
  ),
  plot = p_ancla,
  width = 10,
  height = 6,
  dpi = 300
)

# 17.4 Parámetros de la ecuación de expectativas del candidato.
datos_parametros_M3 <- trayectoria_candidata |>
  dplyr::select(
    fecha,
    `c0_t: intercepto` = c0_t,
    `c_t: sensibilidad a inflación` = c_t,
    `d: persistencia fija` = d_t
  ) |>
  tidyr::pivot_longer(
    cols = -fecha,
    names_to = "parametro",
    values_to = "estimacion"
  )

p_parametros_M3 <- ggplot2::ggplot(
  datos_parametros_M3,
  ggplot2::aes(
    x = fecha,
    y = estimacion,
    linetype = parametro
  )
) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::labs(
    title = "Ecuación de expectativas del Demertzis restringido",
    x = NULL,
    y = "Coeficiente",
    linetype = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = file.path(
    carpeta_salida,
    "19_parametros_expectativas_M3.png"
  ),
  plot = p_parametros_M3,
  width = 10,
  height = 6,
  dpi = 300
)

# 17.5 Criterios de información.
datos_ic <- tabla_comparacion |>
  dplyr::select(
    codigo_modelo,
    AIC,
    BIC
  ) |>
  tidyr::pivot_longer(
    cols = c(AIC, BIC),
    names_to = "criterio",
    values_to = "valor"
  )

p_ic <- ggplot2::ggplot(
  datos_ic,
  ggplot2::aes(
    x = codigo_modelo,
    y = valor,
    linetype = criterio,
    group = criterio
  )
) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::labs(
    title = "Comparación de AIC y BIC entre especificaciones",
    subtitle = "Valores menores indican mejor balance entre ajuste y parsimonia",
    x = "Modelo",
    y = "Valor del criterio",
    linetype = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = file.path(
    carpeta_salida,
    "20_comparacion_AIC_BIC.png"
  ),
  plot = p_ic,
  width = 9,
  height = 5.5,
  dpi = 300
)

# 17.6 ACF de residuos recursivos del candidato.
for (ecuacion in c("inflacion", "expectativa")) {

  residuos_validos <- residuos_M3[, ecuacion]
  residuos_validos <- residuos_validos[is.finite(residuos_validos)]

  grDevices::png(
    filename = file.path(
      carpeta_salida,
      paste0("21_ACF_residuos_M3_", ecuacion, ".png")
    ),
    width = 1800,
    height = 1100,
    res = 180
  )

  stats::acf(
    residuos_validos,
    lag.max = 24,
    main = paste("ACF residuos recursivos M3:", ecuacion)
  )

  grDevices::dev.off()
}


# 18. GUARDAR OBJETOS ----------------------------------------------------------

for (codigo in names(ajustes)) {
  saveRDS(
    ajustes[[codigo]]$modelo,
    file.path(
      carpeta_salida,
      paste0("modelo_", codigo, ".rds")
    )
  )
}


# 19. REGLAS AUTOMÁTICAS DE DECISIÓN ------------------------------------------

fila_M0 <- tabla_comparacion |>
  dplyr::filter(codigo_modelo == "M0")

fila_M1 <- tabla_comparacion |>
  dplyr::filter(codigo_modelo == "M1")

fila_M2 <- tabla_comparacion |>
  dplyr::filter(codigo_modelo == "M2")

fila_M3 <- tabla_comparacion |>
  dplyr::filter(codigo_modelo == "M3")

fila_M4 <- tabla_comparacion |>
  dplyr::filter(codigo_modelo == "M4")

Q_M3 <- tabla_hiperparametros |>
  dplyr::filter(codigo_modelo == "M3")

Q_M2 <- tabla_hiperparametros |>
  dplyr::filter(codigo_modelo == "M2")

p_LR_M3_M0 <- tabla_pruebas_LR |>
  dplyr::filter(
    modelo_restringido == "M0",
    modelo_irrestricto == "M3"
  ) |>
  dplyr::pull(p_valor_chi_cuadrado_referencial)

p_LR_M2_M3_mixto <- tabla_pruebas_LR |>
  dplyr::filter(
    modelo_restringido == "M3",
    modelo_irrestricto == "M2"
  ) |>
  dplyr::pull(p_valor_mixto_una_varianza)

p_LB_M3_expectativa_12 <- tabla_diagnosticos |>
  dplyr::filter(
    codigo_modelo == "M3",
    ecuacion == "expectativa",
    diagnostico == "Ljung-Box",
    rezago == 12
  ) |>
  dplyr::pull(p_valor)

p_LB_M3_inflacion_12 <- tabla_diagnosticos |>
  dplyr::filter(
    codigo_modelo == "M3",
    ecuacion == "inflacion",
    diagnostico == "Ljung-Box",
    rezago == 12
  ) |>
  dplyr::pull(p_valor)

reglas_decision <- tibble::tibble(
  criterio = c(
    "M3 mejora al constante según AIC",
    "M3 mejora al constante según BIC",
    "LR referencial M3 vs M0 menor a 5%",
    "Q_c0 de M3 supera claramente el límite",
    "Q_c de M3 supera claramente el límite",
    "d_t adicional no es necesario: M2 no mejora a M3",
    "Residuos de expectativas M3 sin autocorrelación a 12 meses",
    "Residuos de inflación M3 sin autocorrelación a 12 meses",
    "M3 domina al TVP solo en inflación según BIC",
    "M3 domina al TVP irrestricto según BIC"
  ),
  cumple = c(
    fila_M3$AIC < fila_M0$AIC,
    fila_M3$BIC < fila_M0$BIC,
    is.finite(p_LR_M3_M0) && p_LR_M3_M0 < 0.05,
    Q_M3$Q_c0 > 100 * limite_varianza,
    Q_M3$Q_c > 100 * limite_varianza,
    !is.finite(p_LR_M2_M3_mixto) || p_LR_M2_M3_mixto >= 0.05,
    is.finite(p_LB_M3_expectativa_12) && p_LB_M3_expectativa_12 >= 0.05,
    is.finite(p_LB_M3_inflacion_12) && p_LB_M3_inflacion_12 >= 0.05,
    fila_M3$BIC < fila_M1$BIC,
    fila_M3$BIC < fila_M4$BIC
  ),
  interpretacion = c(
    "La dinámica de credibilidad mejora el ajuste penalizado por complejidad.",
    "La mejora se sostiene bajo una penalización más severa.",
    "Evidencia referencial de variación temporal en c0_t y c_t.",
    "El intercepto de expectativas presenta variación temporal identificable.",
    "La sensibilidad a inflación presenta variación temporal identificable.",
    "La restricción d constante es compatible con los datos.",
    "La ecuación que genera lambda_t está bien especificada dinámicamente.",
    "La ecuación de inflación no deja autocorrelación importante.",
    "La mejora no proviene únicamente de flexibilizar la inflación.",
    "La versión parsimoniosa es preferible al sistema totalmente irrestricto."
  )
)

readr::write_csv(
  reglas_decision,
  file.path(
    carpeta_salida,
    "22_reglas_decision_modelo_principal.csv"
  )
)


# 20. RESUMEN EJECUTIVO --------------------------------------------------------

obtener_valor <- function(tabla, codigo, variable) {
  tabla |>
    dplyr::filter(codigo_modelo == codigo) |>
    dplyr::pull({{ variable }})
}

Q_M3_vector <- c(
  Q_c0 = Q_M3$Q_c0,
  Q_c = Q_M3$Q_c,
  Q_d = Q_M3$Q_d
)

Q_M2_vector <- c(
  Q_c0 = Q_M2$Q_c0,
  Q_c = Q_M2$Q_c,
  Q_d = Q_M2$Q_d
)

lineas_resumen <- c(
  "",
  "======================================================================",
  "DEMERTZIS TVP RESTRINGIDO: PRUEBA DE ESPECIFICACIONES - GUATEMALA",
  "======================================================================",
  paste0("Archivo: ", archivo_entrada),
  paste0(
    "Muestra efectiva: ",
    format(min(datos_modelo$fecha), "%Y-%m"),
    " a ",
    format(max(datos_modelo$fecha), "%Y-%m")
  ),
  paste0("Observaciones: ", n_obs),
  "",
  "MODELO CANDIDATO M3",
  "Inflación: pi_t = a0 + a pi_(t-1) + b e_(t-1) + eps_pi,t",
  "Expectativas: e_t = c0_t + c_t pi_(t-1) + d e_(t-1) + eps_e,t",
  "c0_t y c_t siguen paseos aleatorios; d, a0, a y b son constantes.",
  "lambda_t = 1 - c_t/(1-d)",
  "pi_estrella_t = c0_t/(1-d-c_t)",
  "",
  "COMPARACIÓN DE MODELOS",
  paste0(
    "AIC M0/M1/M2/M3/M4: ",
    paste(
      formato_num_vector(c(
        obtener_valor(tabla_comparacion, "M0", AIC),
        obtener_valor(tabla_comparacion, "M1", AIC),
        obtener_valor(tabla_comparacion, "M2", AIC),
        obtener_valor(tabla_comparacion, "M3", AIC),
        obtener_valor(tabla_comparacion, "M4", AIC)
      )),
      collapse = " / "
    )
  ),
  paste0(
    "BIC M0/M1/M2/M3/M4: ",
    paste(
      formato_num_vector(c(
        obtener_valor(tabla_comparacion, "M0", BIC),
        obtener_valor(tabla_comparacion, "M1", BIC),
        obtener_valor(tabla_comparacion, "M2", BIC),
        obtener_valor(tabla_comparacion, "M3", BIC),
        obtener_valor(tabla_comparacion, "M4", BIC)
      )),
      collapse = " / "
    )
  ),
  paste0(
    "Mejor modelo por AIC: ",
    tabla_comparacion$codigo_modelo[which.min(tabla_comparacion$AIC)]
  ),
  paste0(
    "Mejor modelo por BIC: ",
    tabla_comparacion$codigo_modelo[which.min(tabla_comparacion$BIC)]
  ),
  "",
  "PRUEBAS CENTRALES",
  paste0(
    "LR M3 vs M0: ",
    formato_num(
      tabla_pruebas_LR |>
        dplyr::filter(
          modelo_restringido == "M0",
          modelo_irrestricto == "M3"
        ) |>
        dplyr::pull(estadistico_LR)
    ),
    " | p chi2(2) referencial = ",
    formato_cientifico_vector(p_LR_M3_M0)
  ),
  paste0(
    "LR M2 vs M3 para Q_d=0: ",
    formato_num(
      tabla_pruebas_LR |>
        dplyr::filter(
          modelo_restringido == "M3",
          modelo_irrestricto == "M2"
        ) |>
        dplyr::pull(estadistico_LR)
    ),
    " | p mezcla 0.5 chi1 = ",
    formato_cientifico_vector(p_LR_M2_M3_mixto)
  ),
  "",
  "HIPERPARÁMETROS DE EXPECTATIVAS",
  paste0(
    "M3 Q(c0,c,d): ",
    paste(formato_cientifico_vector(Q_M3_vector), collapse = ", ")
  ),
  paste0(
    "M2 Q(c0,c,d): ",
    paste(formato_cientifico_vector(Q_M2_vector), collapse = ", ")
  ),
  "",
  "PROXY DEL CANDIDATO M3",
  paste0("lambda promedio: ", formato_num(resumen_candidato$lambda_promedio)),
  paste0(
    "lambda inicial/final: ",
    formato_num(resumen_candidato$lambda_inicial),
    " / ",
    formato_num(resumen_candidato$lambda_final)
  ),
  paste0(
    "lambda mínimo/máximo: ",
    formato_num(resumen_candidato$lambda_minimo),
    " / ",
    formato_num(resumen_candidato$lambda_maximo)
  ),
  paste0(
    "% meses lambda fuera de [0,1]: ",
    formato_num(resumen_candidato$porcentaje_lambda_fuera_0_1, 2),
    "%"
  ),
  paste0(
    "pi estrella promedio: ",
    formato_num(resumen_candidato$pi_estrella_promedio)
  ),
  paste0(
    "MAE pi estrella vs meta: ",
    formato_num(resumen_candidato$MAE_pi_estrella_meta)
  ),
  paste0(
    "% meses ancla dentro de banda de meta: ",
    formato_num(resumen_candidato$porcentaje_ancla_dentro_banda_meta, 2),
    "%"
  ),
  "",
  "DIAGNÓSTICOS DEL CANDIDATO M3",
  paste0(
    "Ljung-Box expectativas lag 12, p: ",
    formato_cientifico_vector(p_LB_M3_expectativa_12)
  ),
  paste0(
    "Ljung-Box inflación lag 12, p: ",
    formato_cientifico_vector(p_LB_M3_inflacion_12)
  ),
  paste0(
    "Radio espectral del VAR constante: ",
    formato_num(radio_espectral_M0),
    " | estable: ",
    radio_espectral_M0 < 1
  ),
  "",
  "BOOTSTRAP M3 VS M0",
  paste0("Ejecutado: ", resultado_bootstrap$ejecutado),
  paste0("Réplicas válidas: ", resultado_bootstrap$replicas_validas),
  paste0(
    "p-valor bootstrap: ",
    formato_cientifico_vector(resultado_bootstrap$p_valor_bootstrap)
  ),
  "",
  "CRITERIOS DE DECISIÓN",
  "1. M3 debe mejorar a M0 para demostrar TVP en la credibilidad, no solo en inflación.",
  "2. M2 no debe mejorar materialmente a M3 para justificar d constante.",
  "3. Si M1 domina a M3, la evidencia TVP previa provenía sobre todo de inflación.",
  "4. Si M4 mejora a M3 solo por la ecuación de inflación, M4 se conserva como benchmark.",
  "5. La ecuación de expectativas debe presentar residuos sin autocorrelación relevante.",
  "6. lambda_t y pi_estrella_t no deben depender de denominadores cercanos a cero.",
  "",
  paste0("Resultados guardados en: ", carpeta_salida),
  "======================================================================"
)

writeLines(
  lineas_resumen,
  con = file.path(
    carpeta_salida,
    "00_resumen_ejecutivo.txt"
  )
)

cat(paste(lineas_resumen, collapse = "\n"), "\n")

