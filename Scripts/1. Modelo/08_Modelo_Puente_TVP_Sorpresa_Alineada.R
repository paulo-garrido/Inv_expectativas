# ==============================================================================
# OBJETIVO 2 - SCRIPT 08
# MODELO PUENTE TVP CON SORPRESA INFLACIONARIA ALINEADA
#
# Motivación:
#   Vincular la narrativa de formación de expectativas de Strohsal con el
#   modelo TVP de sensibilidad a información inflacionaria, sin construir
#   necesariamente el indicador A_t = 1 - theta1_t - theta2_t.
#
# ECUACIÓN EMPÍRICA:
#
#   Delta e24_t =
#       alpha
#     + lambda  (meta_t - e24_(t-1))
#     + gamma   (e12_(t-1) - meta_t)
#     + beta_t  N_t
#     + rho     Delta e24_(t-1)
#     + delta   D_sep2012,t
#     + epsilon_t
#
#   beta_t = beta_(t-1) + eta_t
#
# donde N_t es una noticia inflacionaria correctamente alineada:
#
#   N_t^exacta =
#       inflacion_t - E_(t-12)(inflacion_t)
#     = infl_gt_t - exp_inf_12m_(t-12)
#
# Dado que la encuesta del mes t podría levantarse antes de que la inflación
# del mismo mes esté disponible, también se estima una versión conservadora:
#
#   N_t^disponible =
#       inflacion_(t-1) - E_(t-13)(inflacion_(t-1))
#
# Las dos medidas se construyen mediante fechas objetivo, no solamente con
# posiciones de fila, para evitar errores si existieran meses faltantes.
#
# ESTRUCTURA DE ESTIMACIÓN:
#   A. Sorpresa alineada exacta:
#      1. beta constante
#      2. beta_t variante
#      3. beta constante + pulso septiembre 2012
#      4. beta_t variante + pulso septiembre 2012
#
#   B. Sorpresa alineada disponible:
#      1. beta constante
#      2. beta_t variante
#      3. beta constante + pulso septiembre 2012
#      4. beta_t variante + pulso septiembre 2012
#
# La muestra se homogeneiza para que ambas definiciones de noticia se estimen
# sobre exactamente las mismas fechas.
#
# INTERPRETACIÓN:
#   beta_t cercano a cero:
#       menor transmisión de la noticia inflacionaria a e24_t.
#
#   beta_t positivo y elevado:
#       mayor transmisión hacia la expectativa a 24 meses.
#
# IMPORTANTE:
#   Este modelo no estima el porcentaje estructural de anclaje de Strohsal.
#   Estima una sensibilidad dinámica con una narrativa semiestructural:
#   corrección hacia la meta, influencia de expectativas cortas, noticias
#   inflacionarias y dinámica de revisiones.
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
      "\nDirectorio actual: ",
      getwd()
    )
  )
}

carpeta_salida <- file.path(
  getwd(),
  "resultados_objetivo2_modelo_puente_sorpresa_alineada"
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
  "exp_inf_12m",
  "exp_inf_24m",
  "infl_gt"
)

fecha_pulso <- as.Date("2012-09-01")
nivel_confianza <- 0.95
z_critico <- qnorm(1 - (1 - nivel_confianza) / 2)

limite_varianza <- 1e-12


# 2. FUNCIONES GENERALES -------------------------------------------------------

verificar_columnas <- function(datos, columnas) {

  faltan <- setdiff(columnas, names(datos))

  if (length(faltan) > 0) {
    stop(
      "Faltan las siguientes variables: ",
      paste(faltan, collapse = ", ")
    )
  }
}


formato_num <- function(x, digitos = 5) {

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


formato_cientifico <- function(x, digitos = 4) {

  if (
    length(x) == 0 ||
      !is.finite(x[[1]])
  ) {
    return("NA")
  }

  formatC(
    x[[1]],
    digits = digitos,
    format = "e"
  )
}


calcular_AICc <- function(AIC, k, n) {

  if (n <= k + 1) {
    return(NA_real_)
  }

  AIC + (2 * k * (k + 1)) / (n - k - 1)
}


correlacion_segura <- function(x, y, metodo = "pearson") {

  validos <- is.finite(x) & is.finite(y)

  if (sum(validos) < 3) {
    return(NA_real_)
  }

  x_valido <- as.numeric(x[validos])
  y_valido <- as.numeric(y[validos])

  if (
    sd(x_valido) == 0 ||
      sd(y_valido) == 0
  ) {
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      x_valido,
      y_valido,
      method = metodo
    )
  )
}


# 3. CARGA Y ALINEACIÓN DE LA SORPRESA ----------------------------------------

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
  mutate(
    fecha = ymd(fecha)
  ) |>
  arrange(fecha)

if (anyNA(bd$fecha)) {
  stop("Existen fechas que no pudieron convertirse al formato Date.")
}

if (anyDuplicated(bd$fecha) > 0) {
  stop("Existen fechas duplicadas.")
}

# Verificación de continuidad mensual.
saltos_meses <- diff(
  year(bd$fecha) * 12 +
    month(bd$fecha)
)

if (any(saltos_meses != 1)) {
  warning(
    "La base presenta saltos mensuales. ",
    "La sorpresa se alineará por fecha objetivo para evitar errores."
  )
}

# Meta de inflación.
if (!"meta_inflacion" %in% names(bd)) {

  bd <- bd |>
    mutate(
      meta_inflacion = case_when(
        fecha < ymd("2012-01-01") ~ 5.0,
        fecha < ymd("2013-01-01") ~ 4.5,
        TRUE ~ 4.0
      )
    )

} else {

  message("Se utilizará la columna existente 'meta_inflacion'.")
}

# Cada expectativa de 12 meses se asigna a su fecha objetivo exacta.
mapa_pronosticos <- bd |>
  transmute(
    fecha_origen_pronostico = fecha,
    fecha_objetivo =
      fecha %m+% months(12),
    expectativa_12m_origen =
      exp_inf_12m
  )

# Realización y pronóstico para la misma fecha objetivo.
tabla_sorpresa_exacta <- bd |>
  transmute(
    fecha = fecha,
    inflacion_realizada = infl_gt
  ) |>
  left_join(
    mapa_pronosticos |>
      transmute(
        fecha = fecha_objetivo,
        fecha_origen_pronostico,
        expectativa_12m_origen
      ),
    by = "fecha"
  ) |>
  mutate(
    sorpresa_h12_exacta =
      inflacion_realizada -
        expectativa_12m_origen
  )

# La noticia exacta de t se considera disponible para la ecuación de t+1.
mapa_sorpresa_disponible <- tabla_sorpresa_exacta |>
  transmute(
    fecha =
      fecha %m+% months(1),
    fecha_realizacion_sorpresa_disponible =
      fecha %m-% months(1),
    fecha_origen_pronostico_disponible =
      fecha_origen_pronostico,
    expectativa_12m_origen_disponible =
      expectativa_12m_origen,
    inflacion_realizada_disponible =
      inflacion_realizada,
    sorpresa_h12_disponible =
      sorpresa_h12_exacta
  )

bd_modelo <- bd |>
  left_join(
    tabla_sorpresa_exacta |>
      select(
        fecha,
        fecha_origen_pronostico,
        expectativa_12m_origen,
        inflacion_realizada,
        sorpresa_h12_exacta
      ),
    by = "fecha"
  ) |>
  left_join(
    mapa_sorpresa_disponible,
    by = "fecha"
  ) |>
  arrange(fecha) |>
  mutate(
    exp24_lag1 =
      lag(exp_inf_24m, 1),

    exp12_lag1 =
      lag(exp_inf_12m, 1),

    d_exp24 =
      exp_inf_24m -
        exp24_lag1,

    d_exp24_lag1 =
      lag(d_exp24, 1),

    gap_correccion_meta =
      meta_inflacion -
        exp24_lag1,

    gap_expectativa_corta_meta =
      exp12_lag1 -
        meta_inflacion,

    dummy_sep2012 =
      as.numeric(fecha == fecha_pulso)
  )

# Muestra común entre ambas definiciones de sorpresa.
bd_comun <- bd_modelo |>
  filter(
    is.finite(d_exp24),
    is.finite(d_exp24_lag1),
    is.finite(gap_correccion_meta),
    is.finite(gap_expectativa_corta_meta),
    is.finite(sorpresa_h12_exacta),
    is.finite(sorpresa_h12_disponible),
    is.finite(dummy_sep2012)
  )

if (nrow(bd_comun) < 60) {
  stop(
    "La muestra común tiene menos de 60 observaciones. ",
    "Revise la disponibilidad de expectativas a 12 meses."
  )
}

readr::write_csv(
  bd_comun,
  file.path(
    carpeta_salida,
    "01_datos_modelo_puente_muestra_comun.csv"
  )
)

tabla_verificacion_alineacion <- bd_comun |>
  select(
    fecha,
    infl_gt,
    fecha_origen_pronostico,
    expectativa_12m_origen,
    sorpresa_h12_exacta,
    fecha_realizacion_sorpresa_disponible,
    fecha_origen_pronostico_disponible,
    expectativa_12m_origen_disponible,
    inflacion_realizada_disponible,
    sorpresa_h12_disponible
  )

readr::write_csv(
  tabla_verificacion_alineacion,
  file.path(
    carpeta_salida,
    "02_verificacion_alineacion_sorpresas.csv"
  )
)


# 4. DIAGNÓSTICO PRELIMINAR DE REGRESORES -------------------------------------

calcular_vif_manual <- function(datos, variables, nombre_timing) {

  X <- datos |>
    select(all_of(variables))

  filas_validas <- complete.cases(X)
  X <- X[filas_validas, , drop = FALSE]

  tabla_vif <- lapply(
    variables,
    function(variable_objetivo) {

      otras <- setdiff(
        variables,
        variable_objetivo
      )

      formula_aux <- reformulate(
        otras,
        response = variable_objetivo
      )

      modelo_aux <- lm(
        formula_aux,
        data = X
      )

      R2 <- summary(modelo_aux)$r.squared

      tibble(
        definicion_sorpresa =
          nombre_timing,
        variable =
          variable_objetivo,
        R2_auxiliar =
          R2,
        VIF =
          1 / (1 - R2)
      )
    }
  ) |>
    bind_rows()

  X_escalada <- scale(
    as.matrix(X)
  )

  numero_condicion <- tryCatch(
    kappa(X_escalada),
    error = function(e) NA_real_
  )

  tabla_vif |>
    mutate(
      numero_condicion_matriz =
        numero_condicion
    )
}

variables_base <- c(
  "gap_correccion_meta",
  "gap_expectativa_corta_meta",
  "d_exp24_lag1"
)

tabla_colinealidad <- bind_rows(
  calcular_vif_manual(
    datos = bd_comun,
    variables = c(
      variables_base,
      "sorpresa_h12_exacta"
    ),
    nombre_timing =
      "Sorpresa exacta"
  ),

  calcular_vif_manual(
    datos = bd_comun,
    variables = c(
      variables_base,
      "sorpresa_h12_disponible"
    ),
    nombre_timing =
      "Sorpresa disponible"
  )
)

readr::write_csv(
  tabla_colinealidad,
  file.path(
    carpeta_salida,
    "03_diagnostico_colinealidad_regresores.csv"
  )
)


# 5. CONSTRUCCIÓN DEL MODELO KFAS ---------------------------------------------

crear_Z <- function(
    datos,
    variable_sorpresa,
    incluir_pulso
) {

  n <- nrow(datos)
  m <- ifelse(incluir_pulso, 6L, 5L)

  Z_t <- array(
    0,
    dim = c(1, m, n)
  )

  # Estados:
  # 1 alpha: constante
  # 2 lambda: corrección hacia la meta
  # 3 gamma: brecha de expectativa corta frente a meta
  # 4 beta o beta_t: sensibilidad a la sorpresa alineada
  # 5 rho: dinámica de revisiones
  # 6 delta: pulso sep-2012, cuando aplica
  Z_t[1, 1, ] <- 1

  Z_t[1, 2, ] <-
    datos$gap_correccion_meta

  Z_t[1, 3, ] <-
    datos$gap_expectativa_corta_meta

  Z_t[1, 4, ] <-
    datos[[variable_sorpresa]]

  Z_t[1, 5, ] <-
    datos$d_exp24_lag1

  if (incluir_pulso) {
    Z_t[1, 6, ] <-
      datos$dummy_sep2012
  }

  Z_t
}


crear_modelo_constante <- function(
    datos,
    variable_sorpresa,
    incluir_pulso,
    estado_inicial
) {

  n <- nrow(datos)
  m <- length(estado_inicial)

  Z_t <- crear_Z(
    datos = datos,
    variable_sorpresa =
      variable_sorpresa,
    incluir_pulso =
      incluir_pulso
  )

  y_modelo <- as.numeric(
    datos$d_exp24
  )

  KFAS::SSModel(
    y_modelo ~ -1 +
      SSMcustom(
        Z = Z_t,
        T = diag(m),
        R = matrix(
          0,
          nrow = m,
          ncol = 1
        ),
        Q = matrix(
          0,
          nrow = 1,
          ncol = 1
        ),
        a1 = estado_inicial,
        P1 = matrix(
          0,
          nrow = m,
          ncol = m
        ),
        P1inf = diag(m),
        n = n
      ),
    H = matrix(
      NA_real_,
      nrow = 1,
      ncol = 1
    )
  )
}


crear_modelo_tvp <- function(
    datos,
    variable_sorpresa,
    incluir_pulso,
    estado_inicial
) {

  n <- nrow(datos)
  m <- length(estado_inicial)

  Z_t <- crear_Z(
    datos = datos,
    variable_sorpresa =
      variable_sorpresa,
    incluir_pulso =
      incluir_pulso
  )

  R_t <- matrix(
    0,
    nrow = m,
    ncol = 1
  )

  # Solo beta_t evoluciona.
  R_t[4, 1] <- 1

  y_modelo <- as.numeric(
    datos$d_exp24
  )

  KFAS::SSModel(
    y_modelo ~ -1 +
      SSMcustom(
        Z = Z_t,
        T = diag(m),
        R = R_t,
        Q = matrix(
          NA_real_,
          nrow = 1,
          ncol = 1
        ),
        a1 = estado_inicial,
        P1 = matrix(
          0,
          nrow = m,
          ncol = m
        ),
        P1inf = diag(m),
        n = n
      ),
    H = matrix(
      NA_real_,
      nrow = 1,
      ncol = 1
    )
  )
}


actualizar_H <- function(par, modelo) {

  modelo$H[1, 1, 1] <-
    exp(par[1])

  modelo
}


actualizar_H_Q <- function(par, modelo) {

  modelo$H[1, 1, 1] <-
    exp(par[1])

  modelo$Q[1, 1, 1] <-
    exp(par[2])

  modelo
}


# 6. INICIALIZACIÓN Y ESTIMACIÓN ----------------------------------------------

obtener_inicializacion <- function(
    datos,
    variable_sorpresa,
    incluir_pulso
) {

  regresores <- c(
    "gap_correccion_meta",
    "gap_expectativa_corta_meta",
    variable_sorpresa,
    "d_exp24_lag1"
  )

  if (incluir_pulso) {
    regresores <- c(
      regresores,
      "dummy_sep2012"
    )
  }

  formula_ols <- reformulate(
    regresores,
    response = "d_exp24"
  )

  modelo_ols <- lm(
    formula_ols,
    data = datos
  )

  coef_ols <- coef(
    modelo_ols
  )

  if (any(!is.finite(coef_ols))) {
    stop(
      "La inicialización OLS produjo coeficientes no finitos."
    )
  }

  nombres_estado <- c(
    "alpha",
    "lambda_meta",
    "gamma_expectativa_corta",
    "beta_sorpresa",
    "rho_revision"
  )

  estado_inicial <- c(
    unname(
      coef_ols["(Intercept)"]
    ),
    unname(
      coef_ols["gap_correccion_meta"]
    ),
    unname(
      coef_ols["gap_expectativa_corta_meta"]
    ),
    unname(
      coef_ols[variable_sorpresa]
    ),
    unname(
      coef_ols["d_exp24_lag1"]
    )
  )

  if (incluir_pulso) {

    estado_inicial <- c(
      estado_inicial,
      unname(
        coef_ols["dummy_sep2012"]
      )
    )

    nombres_estado <- c(
      nombres_estado,
      "delta_sep2012"
    )
  }

  names(estado_inicial) <-
    nombres_estado

  var_residuo <- max(
    var(
      residuals(modelo_ols),
      na.rm = TRUE
    ),
    1e-8
  )

  list(
    modelo_ols =
      modelo_ols,
    estado_inicial =
      estado_inicial,
    var_residuo =
      var_residuo
  )
}


estimar_constante <- function(
    modelo_sin_estimar,
    var_residuo,
    nombre_modelo
) {

  inicios <- log(
    pmax(
      c(
        var_residuo * 0.25,
        var_residuo * 0.50,
        var_residuo,
        var_residuo * 2,
        var_residuo * 4
      ),
      limite_varianza
    )
  )

  ajustes <- lapply(
    seq_along(inicios),
    function(i) {

      ajuste <- tryCatch(
        KFAS::fitSSM(
          model =
            modelo_sin_estimar,
          inits =
            inicios[i],
          updatefn =
            actualizar_H,
          method =
            "L-BFGS-B",
          lower =
            log(limite_varianza),
          upper =
            log(100),
          control = list(
            maxit = 10000,
            factr = 1e7
          )
        ),
        error = function(e) NULL
      )

      if (is.null(ajuste)) {
        return(NULL)
      }

      ll <- tryCatch(
        as.numeric(
          logLik(ajuste$model)
        ),
        error = function(e) NA_real_
      )

      list(
        intento = i,
        ajuste = ajuste,
        convergencia =
          ajuste$optim.out$convergence,
        mensaje = ifelse(
          is.null(
            ajuste$optim.out$message
          ),
          "",
          as.character(
            ajuste$optim.out$message
          )
        ),
        logLik = ll
      )
    }
  )

  ajustes <- Filter(
    Negate(is.null),
    ajustes
  )

  if (length(ajustes) == 0) {
    stop(
      "No fue posible estimar: ",
      nombre_modelo
    )
  }

  elegibles <- ajustes[
    vapply(
      ajustes,
      function(x) {
        x$convergencia == 0 &&
          is.finite(x$logLik)
      },
      logical(1)
    )
  ]

  if (length(elegibles) == 0) {
    warning(
      "Ningún intento convergió normalmente para ",
      nombre_modelo,
      ". Se usará la mayor log-verosimilitud."
    )
    elegibles <- ajustes
  }

  mejor <- elegibles[[
    which.max(
      vapply(
        elegibles,
        function(x) x$logLik,
        numeric(1)
      )
    )
  ]]

  tabla_intentos <- bind_rows(
    lapply(
      ajustes,
      function(x) {
        tibble(
          modelo =
            nombre_modelo,
          intento =
            x$intento,
          convergencia =
            x$convergencia,
          log_verosimilitud =
            x$logLik,
          mensaje =
            x$mensaje
        )
      }
    )
  )

  list(
    mejor =
      mejor,
    modelo =
      mejor$ajuste$model,
    tabla_intentos =
      tabla_intentos
  )
}


estimar_tvp <- function(
    modelo_sin_estimar,
    var_residuo,
    nombre_modelo
) {

  escalas_Q <- c(
    1e-10,
    1e-8,
    1e-7,
    1e-6,
    1e-5,
    1e-4,
    1e-3
  )

  inicios <- lapply(
    escalas_Q,
    function(q_ini) {
      log(
        c(
          var_residuo,
          q_ini
        )
      )
    }
  )

  inicios <- c(
    inicios,
    list(
      log(
        c(
          var_residuo * 0.5,
          1e-5
        )
      ),
      log(
        c(
          var_residuo * 2,
          1e-5
        )
      )
    )
  )

  ajustes <- lapply(
    seq_along(inicios),
    function(i) {

      ajuste <- tryCatch(
        KFAS::fitSSM(
          model =
            modelo_sin_estimar,
          inits =
            inicios[[i]],
          updatefn =
            actualizar_H_Q,
          method =
            "L-BFGS-B",
          lower =
            rep(
              log(limite_varianza),
              2
            ),
          upper =
            rep(
              log(100),
              2
            ),
          control = list(
            maxit = 20000,
            factr = 1e7
          )
        ),
        error = function(e) NULL
      )

      if (is.null(ajuste)) {
        return(NULL)
      }

      ll <- tryCatch(
        as.numeric(
          logLik(ajuste$model)
        ),
        error = function(e) NA_real_
      )

      list(
        intento = i,
        ajuste = ajuste,
        convergencia =
          ajuste$optim.out$convergence,
        mensaje = ifelse(
          is.null(
            ajuste$optim.out$message
          ),
          "",
          as.character(
            ajuste$optim.out$message
          )
        ),
        logLik = ll
      )
    }
  )

  ajustes <- Filter(
    Negate(is.null),
    ajustes
  )

  if (length(ajustes) == 0) {
    stop(
      "No fue posible estimar: ",
      nombre_modelo
    )
  }

  elegibles <- ajustes[
    vapply(
      ajustes,
      function(x) {
        x$convergencia == 0 &&
          is.finite(x$logLik)
      },
      logical(1)
    )
  ]

  if (length(elegibles) == 0) {
    warning(
      "Ningún intento convergió normalmente para ",
      nombre_modelo,
      ". Se usará la mayor log-verosimilitud."
    )
    elegibles <- ajustes
  }

  mejor <- elegibles[[
    which.max(
      vapply(
        elegibles,
        function(x) x$logLik,
        numeric(1)
      )
    )
  ]]

  tabla_intentos <- bind_rows(
    lapply(
      ajustes,
      function(x) {
        tibble(
          modelo =
            nombre_modelo,
          intento =
            x$intento,
          convergencia =
            x$convergencia,
          log_verosimilitud =
            x$logLik,
          mensaje =
            x$mensaje
        )
      }
    )
  )

  list(
    mejor =
      mejor,
    modelo =
      mejor$ajuste$model,
    tabla_intentos =
      tabla_intentos
  )
}


# 7. EXTRACCIÓN Y DIAGNÓSTICOS -------------------------------------------------

extraer_estado_constante <- function(
    kfs,
    indice_estado
) {

  estimaciones <- as.numeric(
    kfs$alphahat[, indice_estado]
  )

  varianzas <- as.numeric(
    kfs$V[
      indice_estado,
      indice_estado,
    ]
  )

  estimacion <- mean(
    estimaciones[
      is.finite(estimaciones)
    ],
    na.rm = TRUE
  )

  error_estandar <- mean(
    sqrt(
      pmax(
        varianzas[
          is.finite(varianzas)
        ],
        0
      )
    ),
    na.rm = TRUE
  )

  c(
    estimacion =
      estimacion,
    error_estandar =
      error_estandar
  )
}


extraer_residuos_estandarizados <- function(
    kfs,
    n_esperado
) {

  residuos <- tryCatch(
    as.numeric(
      rstandard(
        kfs,
        type = "recursive"
      )
    ),
    error = function(e) {
      rep(
        NA_real_,
        n_esperado
      )
    }
  )

  if (length(residuos) != n_esperado) {
    residuos <- rep(
      NA_real_,
      n_esperado
    )
  }

  residuos
}


extraer_innovaciones <- function(
    kfs,
    n_esperado
) {

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
    innovaciones <- rep(
      NA_real_,
      n_esperado
    )
  }

  innovaciones
}


extraer_varianza_innovaciones <- function(
    kfs,
    n_esperado
) {

  F_t <- tryCatch(
    as.numeric(kfs$F),
    error = function(e) numeric(0)
  )

  if (length(F_t) != n_esperado) {
    F_t <- rep(
      NA_real_,
      n_esperado
    )
  }

  F_t
}


calcular_metricas <- function(
    y,
    innovaciones,
    F_t,
    nombre_modelo
) {

  prediccion <- y - innovaciones

  validos <-
    is.finite(y) &
    is.finite(innovaciones) &
    is.finite(prediccion)

  validos_log <-
    validos &
    is.finite(F_t) &
    F_t > 0

  perdida_log <- if (
    sum(validos_log) > 0
  ) {
    mean(
      0.5 * (
        log(2 * pi) +
          log(F_t[validos_log]) +
          innovaciones[validos_log]^2 /
            F_t[validos_log]
      )
    )
  } else {
    NA_real_
  }

  tibble(
    modelo =
      nombre_modelo,
    observaciones_validas =
      sum(validos),
    RMSE_un_paso =
      sqrt(
        mean(
          innovaciones[validos]^2
        )
      ),
    MAE_un_paso =
      mean(
        abs(
          innovaciones[validos]
        )
      ),
    MedAE_un_paso =
      median(
        abs(
          innovaciones[validos]
        )
      ),
    correlacion_observado_predicho =
      correlacion_segura(
        y[validos],
        prediccion[validos]
      ),
    perdida_logaritmica_media =
      perdida_log
  )
}


calcular_diagnosticos <- function(
    residuos,
    nombre_modelo
) {

  residuos <- residuos[
    is.finite(residuos)
  ]

  rezagos <- c(
    6,
    12,
    18,
    24
  )

  rezagos <- rezagos[
    rezagos < length(residuos)
  ]

  tabla_lb <- bind_rows(
    lapply(
      rezagos,
      function(L) {

        prueba <- Box.test(
          residuos,
          lag = L,
          type = "Ljung-Box",
          fitdf = 0
        )

        tibble(
          modelo =
            nombre_modelo,
          diagnostico =
            "Ljung-Box",
          rezago =
            L,
          estadistico =
            unname(
              prueba$statistic
            ),
          grados_libertad =
            unname(
              prueba$parameter
            ),
          p_valor =
            prueba$p.value,
          conclusion_5pct = ifelse(
            prueba$p.value < 0.05,
            "Se rechaza ausencia de autocorrelación",
            "No se rechaza ausencia de autocorrelación"
          )
        )
      }
    )
  )

  media_r <- mean(residuos)
  desv_r <- sd(residuos)

  if (
    !is.finite(desv_r) ||
      desv_r <= 0
  ) {

    tabla_jb <- tibble(
      modelo =
        nombre_modelo,
      diagnostico =
        "Jarque-Bera aproximado",
      rezago =
        NA_real_,
      estadistico =
        NA_real_,
      grados_libertad =
        2,
      p_valor =
        NA_real_,
      conclusion_5pct =
        "No disponible"
    )

  } else {

    asimetria <- mean(
      (
        (residuos - media_r) /
          desv_r
      )^3
    )

    curtosis <- mean(
      (
        (residuos - media_r) /
          desv_r
      )^4
    )

    JB <- length(residuos) / 6 *
      (
        asimetria^2 +
          (
            (curtosis - 3)^2
          ) / 4
      )

    p_JB <- pchisq(
      JB,
      df = 2,
      lower.tail = FALSE
    )

    tabla_jb <- tibble(
      modelo =
        nombre_modelo,
      diagnostico =
        "Jarque-Bera aproximado",
      rezago =
        NA_real_,
      estadistico =
        JB,
      grados_libertad =
        2,
      p_valor =
        p_JB,
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


# 8. ESTIMACIÓN INTEGRAL POR ESPECIFICACIÓN -----------------------------------

estimar_especificacion <- function(
    datos,
    variable_sorpresa,
    etiqueta_timing,
    incluir_pulso
) {

  etiqueta_pulso <- ifelse(
    incluir_pulso,
    "con pulso",
    "sin pulso"
  )

  nombre_base <- paste(
    etiqueta_timing,
    etiqueta_pulso
  )

  inicializacion <- obtener_inicializacion(
    datos =
      datos,
    variable_sorpresa =
      variable_sorpresa,
    incluir_pulso =
      incluir_pulso
  )

  modelo_constante_sin_estimar <-
    crear_modelo_constante(
      datos =
        datos,
      variable_sorpresa =
        variable_sorpresa,
      incluir_pulso =
        incluir_pulso,
      estado_inicial =
        inicializacion$estado_inicial
    )

  modelo_tvp_sin_estimar <-
    crear_modelo_tvp(
      datos =
        datos,
      variable_sorpresa =
        variable_sorpresa,
      incluir_pulso =
        incluir_pulso,
      estado_inicial =
        inicializacion$estado_inicial
    )

  ajuste_constante <- estimar_constante(
    modelo_sin_estimar =
      modelo_constante_sin_estimar,
    var_residuo =
      inicializacion$var_residuo,
    nombre_modelo =
      paste(
        nombre_base,
        "- beta constante"
      )
  )

  ajuste_tvp <- estimar_tvp(
    modelo_sin_estimar =
      modelo_tvp_sin_estimar,
    var_residuo =
      inicializacion$var_residuo,
    nombre_modelo =
      paste(
        nombre_base,
        "- beta TVP"
      )
  )

  kfs_constante <- KFAS::KFS(
    ajuste_constante$modelo,
    filtering = c(
      "state",
      "mean"
    ),
    smoothing = c(
      "state",
      "mean",
      "disturbance"
    )
  )

  kfs_tvp <- KFAS::KFS(
    ajuste_tvp$modelo,
    filtering = c(
      "state",
      "mean"
    ),
    smoothing = c(
      "state",
      "mean",
      "disturbance"
    )
  )

  nombres_estados <- names(
    inicializacion$estado_inicial
  )

  # Coeficientes del modelo constante.
  tabla_constante <- bind_rows(
    lapply(
      seq_along(nombres_estados),
      function(j) {

        salida <- extraer_estado_constante(
          kfs_constante,
          j
        )

        tibble(
          definicion_sorpresa =
            etiqueta_timing,
          pulso =
            incluir_pulso,
          modelo =
            "Beta constante",
          parametro =
            nombres_estados[j],
          estimacion =
            unname(
              salida["estimacion"]
            ),
          error_estandar =
            unname(
              salida["error_estandar"]
            )
        )
      }
    )
  ) |>
    mutate(
      limite_inferior_95 =
        estimacion -
          z_critico *
            error_estandar,
      limite_superior_95 =
        estimacion +
          z_critico *
            error_estandar,
      estadistico_z =
        estimacion /
          error_estandar,
      p_valor =
        2 * pnorm(
          abs(estadistico_z),
          lower.tail = FALSE
        ),
      significativo_5pct =
        p_valor < 0.05
    )

  # Coeficientes fijos del TVP.
  indices_fijos <- setdiff(
    seq_along(nombres_estados),
    4L
  )

  tabla_fijos_tvp <- bind_rows(
    lapply(
      indices_fijos,
      function(j) {

        salida <- extraer_estado_constante(
          kfs_tvp,
          j
        )

        tibble(
          definicion_sorpresa =
            etiqueta_timing,
          pulso =
            incluir_pulso,
          modelo =
            "Beta TVP",
          parametro =
            nombres_estados[j],
          estimacion =
            unname(
              salida["estimacion"]
            ),
          error_estandar =
            unname(
              salida["error_estandar"]
            )
        )
      }
    )
  ) |>
    mutate(
      limite_inferior_95 =
        estimacion -
          z_critico *
            error_estandar,
      limite_superior_95 =
        estimacion +
          z_critico *
            error_estandar,
      estadistico_z =
        estimacion /
          error_estandar,
      p_valor =
        2 * pnorm(
          abs(estadistico_z),
          lower.tail = FALSE
        ),
      significativo_5pct =
        p_valor < 0.05
    )

  # Trayectoria beta_t.
  beta_t <- as.numeric(
    kfs_tvp$alphahat[, 4]
  )

  var_beta_t <- pmax(
    as.numeric(
      kfs_tvp$V[4, 4, ]
    ),
    0
  )

  se_beta_t <- sqrt(
    var_beta_t
  )

  beta_constante <- tabla_constante |>
    filter(
      parametro ==
        "beta_sorpresa"
    ) |>
    pull(estimacion)

  trayectoria <- datos |>
    transmute(
      fecha,
      definicion_sorpresa =
        etiqueta_timing,
      pulso =
        incluir_pulso,
      variable_sorpresa =
        variable_sorpresa,
      sorpresa =
        .data[[variable_sorpresa]],
      beta_t =
        beta_t,
      beta_t_se =
        se_beta_t,
      beta_t_li_95 =
        beta_t -
          z_critico * se_beta_t,
      beta_t_ls_95 =
        beta_t +
          z_critico * se_beta_t,
      beta_constante =
        beta_constante,
      beta_t_significativo_positivo =
        beta_t_li_95 > 0,
      beta_t_intervalo_incluye_cero =
        beta_t_li_95 <= 0 &
        beta_t_ls_95 >= 0
    )

  H_constante <- as.numeric(
    ajuste_constante$modelo$H[
      1,
      1,
      1
    ]
  )

  H_tvp <- as.numeric(
    ajuste_tvp$modelo$H[
      1,
      1,
      1
    ]
  )

  Q_beta <- as.numeric(
    ajuste_tvp$modelo$Q[
      1,
      1,
      1
    ]
  )

  tabla_varianzas <- bind_rows(
    tibble(
      definicion_sorpresa =
        etiqueta_timing,
      pulso =
        incluir_pulso,
      modelo =
        "Beta constante",
      H =
        H_constante,
      Q_beta =
        0,
      Q_beta_sobre_H =
        0,
      sd_innovacion_beta =
        0,
      Q_beta_en_limite =
        FALSE
    ),

    tibble(
      definicion_sorpresa =
        etiqueta_timing,
      pulso =
        incluir_pulso,
      modelo =
        "Beta TVP",
      H =
        H_tvp,
      Q_beta =
        Q_beta,
      Q_beta_sobre_H =
        Q_beta / H_tvp,
      sd_innovacion_beta =
        sqrt(Q_beta),
      Q_beta_en_limite =
        Q_beta <=
          limite_varianza * 10
    )
  )

  # Predicción.
  n_obs <- nrow(datos)
  y <- datos$d_exp24

  innov_const <- extraer_innovaciones(
    kfs_constante,
    n_obs
  )

  innov_tvp <- extraer_innovaciones(
    kfs_tvp,
    n_obs
  )

  F_const <- extraer_varianza_innovaciones(
    kfs_constante,
    n_obs
  )

  F_tvp <- extraer_varianza_innovaciones(
    kfs_tvp,
    n_obs
  )

  metricas <- bind_rows(
    calcular_metricas(
      y =
        y,
      innovaciones =
        innov_const,
      F_t =
        F_const,
      nombre_modelo =
        paste(
          nombre_base,
          "- beta constante"
        )
    ),

    calcular_metricas(
      y =
        y,
      innovaciones =
        innov_tvp,
      F_t =
        F_tvp,
      nombre_modelo =
        paste(
          nombre_base,
          "- beta TVP"
        )
    )
  ) |>
    mutate(
      definicion_sorpresa =
        etiqueta_timing,
      pulso =
        incluir_pulso,
      .before = 1
    )

  # Comparación de información.
  ll_const <- as.numeric(
    logLik(
      ajuste_constante$modelo
    )
  )

  ll_tvp <- as.numeric(
    logLik(
      ajuste_tvp$modelo
    )
  )

  # Parámetros:
  # Sin pulso:
  #   constante = alpha, lambda, gamma, beta, rho, H = 6
  #   TVP añade Q_beta = 7
  #
  # Con pulso:
  #   constante = 7
  #   TVP = 8
  k_const <- ifelse(
    incluir_pulso,
    7,
    6
  )

  k_tvp <- k_const + 1

  AIC_const <-
    -2 * ll_const +
    2 * k_const

  AIC_tvp <-
    -2 * ll_tvp +
    2 * k_tvp

  AICc_const <- calcular_AICc(
    AIC_const,
    k_const,
    n_obs
  )

  AICc_tvp <- calcular_AICc(
    AIC_tvp,
    k_tvp,
    n_obs
  )

  BIC_const <-
    -2 * ll_const +
    log(n_obs) * k_const

  BIC_tvp <-
    -2 * ll_tvp +
    log(n_obs) * k_tvp

  comparacion <- tibble(
    definicion_sorpresa =
      etiqueta_timing,
    pulso =
      incluir_pulso,
    modelo = c(
      "Beta constante",
      "Beta TVP"
    ),
    observaciones =
      n_obs,
    parametros_contados = c(
      k_const,
      k_tvp
    ),
    log_verosimilitud = c(
      ll_const,
      ll_tvp
    ),
    AIC = c(
      AIC_const,
      AIC_tvp
    ),
    AICc = c(
      AICc_const,
      AICc_tvp
    ),
    BIC = c(
      BIC_const,
      BIC_tvp
    ),
    H = c(
      H_constante,
      H_tvp
    ),
    Q_beta = c(
      0,
      Q_beta
    )
  ) |>
    mutate(
      delta_AIC =
        AIC - min(AIC),
      peso_Akaike =
        exp(-0.5 * delta_AIC) /
        sum(
          exp(-0.5 * delta_AIC)
        ),
      delta_BIC =
        BIC - min(BIC),
      peso_BIC_aproximado =
        exp(-0.5 * delta_BIC) /
        sum(
          exp(-0.5 * delta_BIC)
        )
    )

  LR <- max(
    0,
    2 * (
      ll_tvp -
        ll_const
    )
  )

  p_LR_convencional <- pchisq(
    LR,
    df = 1,
    lower.tail = FALSE
  )

  p_LR_frontera <- ifelse(
    LR > 0,
    0.5 *
      p_LR_convencional,
    1
  )

  tabla_LR <- tibble(
    definicion_sorpresa =
      etiqueta_timing,
    pulso =
      incluir_pulso,
    hipotesis_nula =
      "H0: Q_beta = 0",
    estadistico_LR =
      LR,
    p_valor_chi2_1 =
      p_LR_convencional,
    p_valor_mezcla_frontera_50_50 =
      p_LR_frontera,
    conclusion_5pct = ifelse(
      p_LR_frontera < 0.05,
      "Se rechaza H0; evidencia a favor de beta_t variante",
      "No se rechaza H0; evidencia insuficiente de variación temporal"
    )
  )

  # Diagnósticos.
  residuos_const <- extraer_residuos_estandarizados(
    kfs_constante,
    n_obs
  )

  residuos_tvp <- extraer_residuos_estandarizados(
    kfs_tvp,
    n_obs
  )

  diagnosticos <- bind_rows(
    calcular_diagnosticos(
      residuos_const,
      paste(
        nombre_base,
        "- beta constante"
      )
    ),

    calcular_diagnosticos(
      residuos_tvp,
      paste(
        nombre_base,
        "- beta TVP"
      )
    )
  ) |>
    mutate(
      definicion_sorpresa =
        etiqueta_timing,
      pulso =
        incluir_pulso,
      .before = 1
    )

  # Resumen beta.
  indice_min <- which.min(
    trayectoria$beta_t
  )

  indice_max <- which.max(
    trayectoria$beta_t
  )

  resumen_beta <- tibble(
    definicion_sorpresa =
      etiqueta_timing,
    pulso =
      incluir_pulso,
    observaciones =
      nrow(trayectoria),
    fecha_inicial =
      first(trayectoria$fecha),
    fecha_final =
      last(trayectoria$fecha),
    beta_constante =
      beta_constante,
    beta_t_promedio =
      mean(
        trayectoria$beta_t,
        na.rm = TRUE
      ),
    beta_t_inicial =
      first(
        trayectoria$beta_t
      ),
    beta_t_final =
      last(
        trayectoria$beta_t
      ),
    cambio_total_beta =
      last(
        trayectoria$beta_t
      ) -
      first(
        trayectoria$beta_t
      ),
    reduccion_porcentual = ifelse(
      first(
        trayectoria$beta_t
      ) != 0,
      100 * (
        1 -
          last(
            trayectoria$beta_t
          ) /
          first(
            trayectoria$beta_t
          )
      ),
      NA_real_
    ),
    beta_t_minimo =
      trayectoria$beta_t[
        indice_min
      ],
    fecha_beta_minimo =
      trayectoria$fecha[
        indice_min
      ],
    beta_t_maximo =
      trayectoria$beta_t[
        indice_max
      ],
    fecha_beta_maximo =
      trayectoria$fecha[
        indice_max
      ],
    porcentaje_beta_positivo =
      100 * mean(
        trayectoria$beta_t > 0,
        na.rm = TRUE
      ),
    porcentaje_beta_significativo_positivo =
      100 * mean(
        trayectoria$
          beta_t_significativo_positivo,
        na.rm = TRUE
      ),
    porcentaje_intervalo_incluye_cero =
      100 * mean(
        trayectoria$
          beta_t_intervalo_incluye_cero,
        na.rm = TRUE
      ),
    H =
      H_tvp,
    Q_beta =
      Q_beta,
    Q_beta_sobre_H =
      Q_beta / H_tvp
  )

  list(
    nombre_base =
      nombre_base,
    etiqueta_timing =
      etiqueta_timing,
    variable_sorpresa =
      variable_sorpresa,
    incluir_pulso =
      incluir_pulso,
    ajuste_constante =
      ajuste_constante,
    ajuste_tvp =
      ajuste_tvp,
    kfs_constante =
      kfs_constante,
    kfs_tvp =
      kfs_tvp,
    tabla_intentos =
      bind_rows(
        ajuste_constante$
          tabla_intentos,
        ajuste_tvp$
          tabla_intentos
      ),
    tabla_coeficientes =
      bind_rows(
        tabla_constante,
        tabla_fijos_tvp
      ),
    tabla_varianzas =
      tabla_varianzas,
    trayectoria =
      trayectoria,
    metricas =
      metricas,
    comparacion =
      comparacion,
    tabla_LR =
      tabla_LR,
    diagnosticos =
      diagnosticos,
    resumen_beta =
      resumen_beta,
    residuos_constante =
      residuos_const,
    residuos_tvp =
      residuos_tvp
  )
}


# 9. ESTIMAR LAS OCHO ESPECIFICACIONES ----------------------------------------

configuraciones <- tribble(
  ~variable_sorpresa,
  ~etiqueta_timing,
  ~incluir_pulso,

  "sorpresa_h12_exacta",
  "Sorpresa exacta",
  FALSE,

  "sorpresa_h12_exacta",
  "Sorpresa exacta",
  TRUE,

  "sorpresa_h12_disponible",
  "Sorpresa disponible",
  FALSE,

  "sorpresa_h12_disponible",
  "Sorpresa disponible",
  TRUE
)

resultados <- pmap(
  configuraciones,
  function(
      variable_sorpresa,
      etiqueta_timing,
      incluir_pulso
  ) {

    estimar_especificacion(
      datos =
        bd_comun,
      variable_sorpresa =
        variable_sorpresa,
      etiqueta_timing =
        etiqueta_timing,
      incluir_pulso =
        incluir_pulso
    )
  }
)


# 10. CONSOLIDAR Y EXPORTAR ----------------------------------------------------

tabla_intentos <- bind_rows(
  lapply(
    resultados,
    function(x) x$tabla_intentos
  )
)

tabla_coeficientes <- bind_rows(
  lapply(
    resultados,
    function(x) x$tabla_coeficientes
  )
)

tabla_varianzas <- bind_rows(
  lapply(
    resultados,
    function(x) x$tabla_varianzas
  )
)

tabla_trayectorias <- bind_rows(
  lapply(
    resultados,
    function(x) x$trayectoria
  )
)

tabla_metricas <- bind_rows(
  lapply(
    resultados,
    function(x) x$metricas
  )
)

tabla_comparacion <- bind_rows(
  lapply(
    resultados,
    function(x) x$comparacion
  )
)

tabla_LR <- bind_rows(
  lapply(
    resultados,
    function(x) x$tabla_LR
  )
)

tabla_diagnosticos <- bind_rows(
  lapply(
    resultados,
    function(x) x$diagnosticos
  )
)

tabla_resumen_beta <- bind_rows(
  lapply(
    resultados,
    function(x) x$resumen_beta
  )
)

readr::write_csv(
  tabla_intentos,
  file.path(
    carpeta_salida,
    "04_intentos_estimacion.csv"
  )
)

readr::write_csv(
  tabla_coeficientes,
  file.path(
    carpeta_salida,
    "05_coeficientes_fijos_y_constantes.csv"
  )
)

readr::write_csv(
  tabla_varianzas,
  file.path(
    carpeta_salida,
    "06_varianzas_estimadas.csv"
  )
)

readr::write_csv(
  tabla_comparacion,
  file.path(
    carpeta_salida,
    "07_comparacion_constante_vs_tvp.csv"
  )
)

readr::write_csv(
  tabla_LR,
  file.path(
    carpeta_salida,
    "08_prueba_LR_variacion_beta.csv"
  )
)

readr::write_csv(
  tabla_metricas,
  file.path(
    carpeta_salida,
    "09_metricas_prediccion_un_paso.csv"
  )
)

readr::write_csv(
  tabla_diagnosticos,
  file.path(
    carpeta_salida,
    "10_diagnosticos_residuales.csv"
  )
)

readr::write_csv(
  tabla_trayectorias,
  file.path(
    carpeta_salida,
    "11_trayectorias_beta_t.csv"
  )
)

readr::write_csv(
  tabla_resumen_beta,
  file.path(
    carpeta_salida,
    "12_resumen_trayectorias_beta.csv"
  )
)


# 11. COMPARACIONES ENTRE DEFINICIONES ----------------------------------------

trayectoria_exacta_base <- tabla_trayectorias |>
  filter(
    definicion_sorpresa ==
      "Sorpresa exacta",
    pulso == FALSE
  ) |>
  select(
    fecha,
    beta_exacta =
      beta_t
  )

trayectoria_disponible_base <-
  tabla_trayectorias |>
  filter(
    definicion_sorpresa ==
      "Sorpresa disponible",
    pulso == FALSE
  ) |>
  select(
    fecha,
    beta_disponible =
      beta_t
  )

comparacion_timing <- inner_join(
  trayectoria_exacta_base,
  trayectoria_disponible_base,
  by = "fecha"
) |>
  arrange(fecha) |>
  mutate(
    diferencia =
      beta_exacta -
        beta_disponible
  )

tabla_comparacion_timing <- tibble(
  observaciones =
    nrow(comparacion_timing),
  correlacion_pearson =
    correlacion_segura(
      comparacion_timing$beta_exacta,
      comparacion_timing$beta_disponible,
      "pearson"
    ),
  correlacion_spearman =
    correlacion_segura(
      comparacion_timing$beta_exacta,
      comparacion_timing$beta_disponible,
      "spearman"
    ),
  diferencia_absoluta_media =
    mean(
      abs(
        comparacion_timing$diferencia
      ),
      na.rm = TRUE
    ),
  diferencia_absoluta_maxima =
    max(
      abs(
        comparacion_timing$diferencia
      ),
      na.rm = TRUE
    )
)

readr::write_csv(
  tabla_comparacion_timing,
  file.path(
    carpeta_salida,
    "13_comparacion_timing_sorpresa.csv"
  )
)

# Robustez del pulso dentro de cada timing.
comparar_pulso <- function(
    tabla,
    timing
) {

  base <- tabla |>
    filter(
      definicion_sorpresa ==
        timing,
      pulso == FALSE
    ) |>
    select(
      fecha,
      beta_base =
        beta_t
    )

  pulso <- tabla |>
    filter(
      definicion_sorpresa ==
        timing,
      pulso == TRUE
    ) |>
    select(
      fecha,
      beta_pulso =
        beta_t
    )

  unidos <- inner_join(
    base,
    pulso,
    by = "fecha"
  )

  tibble(
    definicion_sorpresa =
      timing,
    observaciones =
      nrow(unidos),
    correlacion_pearson =
      correlacion_segura(
        unidos$beta_base,
        unidos$beta_pulso
      ),
    diferencia_absoluta_media =
      mean(
        abs(
          unidos$beta_base -
            unidos$beta_pulso
        ),
        na.rm = TRUE
      ),
    diferencia_absoluta_maxima =
      max(
        abs(
          unidos$beta_base -
            unidos$beta_pulso
        ),
        na.rm = TRUE
      )
  )
}

tabla_robustez_pulso <- bind_rows(
  comparar_pulso(
    tabla_trayectorias,
    "Sorpresa exacta"
  ),
  comparar_pulso(
    tabla_trayectorias,
    "Sorpresa disponible"
  )
)

readr::write_csv(
  tabla_robustez_pulso,
  file.path(
    carpeta_salida,
    "14_robustez_trayectoria_con_pulso.csv"
  )
)


# 12. GRÁFICAS -----------------------------------------------------------------

# 12.1 Trayectorias base exacta y disponible.
datos_grafica_timing <- tabla_trayectorias |>
  filter(
    pulso == FALSE
  )

p_timing <- ggplot(
  datos_grafica_timing,
  aes(
    x = fecha,
    y = beta_t,
    linetype =
      definicion_sorpresa
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    title =
      "Sensibilidad TVP a noticias inflacionarias alineadas",
    subtitle =
      "Comparación entre sorpresa exacta y sorpresa disponible",
    x = NULL,
    y = expression(beta[t]),
    linetype = NULL
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    plot.title.position =
      "plot",
    panel.grid.minor =
      element_blank(),
    legend.position =
      "bottom"
  )

ggsave(
  file.path(
    carpeta_salida,
    "15_beta_t_exacta_vs_disponible.png"
  ),
  plot =
    p_timing,
  width =
    10,
  height =
    6,
  dpi =
    300
)


# 12.2 Trayectoria exacta con bandas.
datos_exacta <- tabla_trayectorias |>
  filter(
    definicion_sorpresa ==
      "Sorpresa exacta",
    pulso == FALSE
  )

p_exacta <- ggplot(
  datos_exacta,
  aes(
    x = fecha,
    y = beta_t
  )
) +
  geom_ribbon(
    aes(
      ymin =
        beta_t_li_95,
      ymax =
        beta_t_ls_95
    ),
    alpha = 0.18
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title =
      "Modelo puente TVP: sorpresa alineada exacta",
    subtitle =
      expression(
        N[t] ==
          pi[t] -
          E[t-12](pi[t])
      ),
    x = NULL,
    y = expression(beta[t])
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    plot.title.position =
      "plot",
    panel.grid.minor =
      element_blank()
  )

ggsave(
  file.path(
    carpeta_salida,
    "16_beta_t_sorpresa_exacta.png"
  ),
  plot =
    p_exacta,
  width =
    10,
  height =
    6,
  dpi =
    300
)


# 12.3 Trayectoria disponible con bandas.
datos_disponible <- tabla_trayectorias |>
  filter(
    definicion_sorpresa ==
      "Sorpresa disponible",
    pulso == FALSE
  )

p_disponible <- ggplot(
  datos_disponible,
  aes(
    x = fecha,
    y = beta_t
  )
) +
  geom_ribbon(
    aes(
      ymin =
        beta_t_li_95,
      ymax =
        beta_t_ls_95
    ),
    alpha = 0.18
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title =
      "Modelo puente TVP: sorpresa alineada disponible",
    subtitle =
      expression(
        N[t] ==
          pi[t-1] -
          E[t-13](pi[t-1])
      ),
    x = NULL,
    y = expression(beta[t])
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    plot.title.position =
      "plot",
    panel.grid.minor =
      element_blank()
  )

ggsave(
  file.path(
    carpeta_salida,
    "17_beta_t_sorpresa_disponible.png"
  ),
  plot =
    p_disponible,
  width =
    10,
  height =
    6,
  dpi =
    300
)


# 12.4 Robustez al pulso.
p_pulso <- ggplot(
  tabla_trayectorias,
  aes(
    x = fecha,
    y = beta_t,
    linetype =
      interaction(
        definicion_sorpresa,
        pulso,
        sep = " | pulso = "
      )
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  labs(
    title =
      "Robustez del modelo puente al pulso de septiembre de 2012",
    x = NULL,
    y = expression(beta[t]),
    linetype = NULL
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    plot.title.position =
      "plot",
    panel.grid.minor =
      element_blank(),
    legend.position =
      "bottom"
  )

ggsave(
  file.path(
    carpeta_salida,
    "18_robustez_beta_t_pulso.png"
  ),
  plot =
    p_pulso,
  width =
    10,
  height =
    6,
  dpi =
    300
)


# 12.5 ACF de residuos TVP por especificación.
for (i in seq_along(resultados)) {

  resultado <- resultados[[i]]

  etiqueta_archivo <- paste0(
    ifelse(
      resultado$etiqueta_timing ==
        "Sorpresa exacta",
      "exacta",
      "disponible"
    ),
    "_",
    ifelse(
      resultado$incluir_pulso,
      "con_pulso",
      "sin_pulso"
    )
  )

  residuos_validos <- resultado$
    residuos_tvp[
      is.finite(
        resultado$residuos_tvp
      )
    ]

  png(
    filename = file.path(
      carpeta_salida,
      paste0(
        "19_acf_residuos_tvp_",
        etiqueta_archivo,
        ".png"
      )
    ),
    width = 1800,
    height = 1100,
    res = 180
  )

  acf(
    residuos_validos,
    lag.max = 24,
    main = paste(
      "ACF residuos TVP:",
      resultado$nombre_base
    )
  )

  dev.off()
}


# 13. GUARDAR MODELOS ----------------------------------------------------------

for (i in seq_along(resultados)) {

  resultado <- resultados[[i]]

  etiqueta_archivo <- paste0(
    ifelse(
      resultado$etiqueta_timing ==
        "Sorpresa exacta",
      "exacta",
      "disponible"
    ),
    "_",
    ifelse(
      resultado$incluir_pulso,
      "con_pulso",
      "sin_pulso"
    )
  )

  saveRDS(
    resultado$ajuste_constante$modelo,
    file.path(
      carpeta_salida,
      paste0(
        "20_modelo_constante_",
        etiqueta_archivo,
        ".rds"
      )
    )
  )

  saveRDS(
    resultado$ajuste_tvp$modelo,
    file.path(
      carpeta_salida,
      paste0(
        "21_modelo_tvp_",
        etiqueta_archivo,
        ".rds"
      )
    )
  )
}


# 14. RESUMEN EJECUTIVO --------------------------------------------------------

lineas_resumen <- c(
  "",
  "======================================================================",
  "MODELO PUENTE TVP CON SORPRESA INFLACIONARIA ALINEADA",
  "======================================================================",
  paste0(
    "Archivo: ",
    archivo_entrada
  ),
  paste0(
    "Muestra común: ",
    format(
      min(bd_comun$fecha),
      "%Y-%m"
    ),
    " a ",
    format(
      max(bd_comun$fecha),
      "%Y-%m"
    )
  ),
  paste0(
    "Observaciones: ",
    nrow(bd_comun)
  ),
  "",
  "ECUACIÓN",
  paste(
    "Delta e24_t = alpha",
    "+ lambda(meta_t-e24_t-1)",
    "+ gamma(e12_t-1-meta_t)",
    "+ beta_t N_t",
    "+ rho Delta e24_t-1",
    "+ pulso opcional",
    "+ error."
  ),
  "",
  "DEFINICIONES DE NOTICIA",
  paste(
    "Exacta: inflacion_t - expectativa 12m formulada en t-12."
  ),
  paste(
    "Disponible: inflacion_t-1 - expectativa 12m formulada en t-13."
  ),
  ""
)

for (resultado in resultados) {

  comp_const <- resultado$comparacion |>
    filter(
      modelo ==
        "Beta constante"
    )

  comp_tvp <- resultado$comparacion |>
    filter(
      modelo ==
        "Beta TVP"
    )

  lr <- resultado$tabla_LR
  resumen <- resultado$resumen_beta

  coef_tvp <- resultado$tabla_coeficientes |>
    filter(
      modelo ==
        "Beta TVP"
    )

  alpha <- coef_tvp |>
    filter(
      parametro ==
        "alpha"
    )

  lambda <- coef_tvp |>
    filter(
      parametro ==
        "lambda_meta"
    )

  gamma <- coef_tvp |>
    filter(
      parametro ==
        "gamma_expectativa_corta"
    )

  rho <- coef_tvp |>
    filter(
      parametro ==
        "rho_revision"
    )

  lineas_resumen <- c(
    lineas_resumen,
    "----------------------------------------------------------------------",
    resultado$nombre_base,
    "----------------------------------------------------------------------",
    paste0(
      "Convergencia TVP: ",
      resultado$ajuste_tvp$
        mejor$convergencia
    ),
    paste0(
      "LogLik constante/TVP: ",
      formato_num(
        comp_const$
          log_verosimilitud
      ),
      " / ",
      formato_num(
        comp_tvp$
          log_verosimilitud
      )
    ),
    paste0(
      "AIC constante/TVP: ",
      formato_num(
        comp_const$AIC
      ),
      " / ",
      formato_num(
        comp_tvp$AIC
      )
    ),
    paste0(
      "BIC constante/TVP: ",
      formato_num(
        comp_const$BIC
      ),
      " / ",
      formato_num(
        comp_tvp$BIC
      )
    ),
    paste0(
      "LR: ",
      formato_num(
        lr$estadistico_LR
      )
    ),
    paste0(
      "p frontera: ",
      formato_cientifico(
        lr$
          p_valor_mezcla_frontera_50_50
      )
    ),
    paste0(
      "alpha: ",
      formato_num(
        alpha$estimacion
      ),
      " | p = ",
      formato_cientifico(
        alpha$p_valor
      )
    ),
    paste0(
      "lambda meta: ",
      formato_num(
        lambda$estimacion
      ),
      " | p = ",
      formato_cientifico(
        lambda$p_valor
      )
    ),
    paste0(
      "gamma expectativa corta: ",
      formato_num(
        gamma$estimacion
      ),
      " | p = ",
      formato_cientifico(
        gamma$p_valor
      )
    ),
    paste0(
      "rho revision: ",
      formato_num(
        rho$estimacion
      ),
      " | p = ",
      formato_cientifico(
        rho$p_valor
      )
    ),
    paste0(
      "Q_beta: ",
      formato_cientifico(
        resumen$Q_beta
      )
    ),
    paste0(
      "Q_beta/H: ",
      formato_cientifico(
        resumen$
          Q_beta_sobre_H
      )
    ),
    paste0(
      "Beta inicial/final: ",
      formato_num(
        resumen$
          beta_t_inicial
      ),
      " / ",
      formato_num(
        resumen$
          beta_t_final
      )
    ),
    paste0(
      "Reducción beta: ",
      formato_num(
        resumen$
          reduccion_porcentual,
        2
      ),
      "%"
    ),
    paste0(
      "% meses beta significativamente positivo: ",
      formato_num(
        resumen$
          porcentaje_beta_significativo_positivo,
        2
      ),
      "%"
    ),
    ""
  )
}

lineas_resumen <- c(
  lineas_resumen,
  "COMPARACIÓN ENTRE DEFINICIONES",
  paste0(
    "Correlación Pearson beta exacta/disponible: ",
    formato_num(
      tabla_comparacion_timing$
        correlacion_pearson
    )
  ),
  paste0(
    "Correlación Spearman beta exacta/disponible: ",
    formato_num(
      tabla_comparacion_timing$
        correlacion_spearman
    )
  ),
  "",
  "CRITERIO DE INTERPRETACIÓN",
  paste(
    "La definición exacta alinea correctamente pronóstico y realización."
  ),
  paste(
    "La definición disponible evita usar inflación no observable",
    "al momento de la encuesta, si la publicación ocurre después."
  ),
  paste(
    "La selección final debe considerar resultados econométricos",
    "y la cronología documental de la encuesta y del IPC."
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
