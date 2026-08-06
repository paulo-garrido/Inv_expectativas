# ==============================================================================
# OBJETIVO 2 - SCRIPT 07 - VERSIÓN 2
# MODELO DE STROHSAL, MELNICK Y NAUTZ ADAPTADO A GUATEMALA
#
# Fuente metodológica:
# Strohsal, T., Melnick, R. y Nautz, D. (2016),
# "The Time-Varying Degree of Inflation Expectations Anchoring",
# Journal of Macroeconomics, 48, 62-71.
#
# ECUACIÓN ORIGINAL:
#
# Delta e^L_t =
#   (1 - theta1_t - theta2_t) (meta_t - e^L_(t-1))
#   + theta1_t (inflacion_(t-1) - e^L_(t-1))
#   + theta2_t (e^S_(t-1) - e^L_(t-1))
#   + beta' X_t
#   + u_t
#
# theta1_t = theta1_(t-1) + epsilon1_t
# theta2_t = theta2_(t-1) + epsilon2_t
# u_t      = alpha0 + alpha1 u_(t-1) + nu_t
#
# Var(epsilon1_t, epsilon2_t, nu_t)' = diagonal.
#
# PARÁMETRO DE ANCLAJE:
#
# A_t = 1 - theta1_t - theta2_t
#
# Multiplicado por 100, A_t puede interpretarse, dentro de la estructura del
# modelo, como el porcentaje de anclaje en el período t.
#
# REPRESENTACIÓN EQUIVALENTE UTILIZADA EN EL CÓDIGO:
#
# y*_t = Delta e24_t - (meta_t - e24_(t-1))
#
# y*_t =
#   theta1_t (inflacion_(t-1) - meta_t)
#   + theta2_t (e12_(t-1) - meta_t)
#   + beta' X_t
#   + u_t
#
# Esta transformación es algebraicamente equivalente a la ecuación anterior.
#
# ADAPTACIÓN DE DATOS:
# - Expectativa de largo plazo: exp_inf_24m.
# - Expectativa de corto plazo: exp_inf_12m.
# - Inflación observada: infl_gt.
# - Meta: meta_inflacion; si no existe, se construye:
#       2010-2011 = 5.0
#       2012      = 4.5
#       2013-...  = 4.0
#
# ESPECIFICACIONES:
# 1. "Núcleo Strohsal": sin controles adicionales X_t.
# 2. "Strohsal + TC": incluye tcdep_var rezagada como control fijo.
#
# El artículo utiliza controles de riesgo y liquidez porque trabaja con
# expectativas de mercado. Esos controles no están disponibles ni son
# necesariamente pertinentes para expectativas de encuesta. La segunda
# especificación conserva la forma beta'X_t usando un control relevante
# para Guatemala, pero debe presentarse como adaptación.
#
# NOTAS:
# - No se imponen restricciones theta1 >= 0, theta2 >= 0 ni A_t entre 0 y 1,
#   porque el artículo tampoco las impone en la estimación.
# - H = 0: el error u_t forma parte del vector de estados, como en el paper.
# - El modelo se estima por máxima verosimilitud mediante filtro de Kalman.
# ==============================================================================


# 0. PAQUETES ------------------------------------------------------------------

paquetes <- c(
  "dplyr",
  "tidyr",
  "readr",
  "lubridate",
  "ggplot2",
  "tibble",
  "purrr",
  "sandwich",
  "lmtest",
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

archivos_datos_candidatos <- c(
  "bd_sample_modelo.csv",
  "bd_sample_modelo(1).csv"
)

archivo_datos <- archivos_datos_candidatos[
  file.exists(archivos_datos_candidatos)
][1]

if (is.na(archivo_datos)) {
  stop(
    paste0(
      "No se encontró la base de datos.\n",
      "Rutas revisadas:\n",
      paste0("- ", archivos_datos_candidatos, collapse = "\n"),
      "\nDirectorio actual: ",
      getwd()
    )
  )
}

carpeta_salida <- file.path(
  getwd(),
  "resultados_objetivo2_modelo_strohsal"
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

candidatas_tc <- c(
  "tcdep_var",
  "depreciacion_tc",
  "tc_var",
  "q_t",
  "tipo_cambio_var",
  "dep_tc"
)

nivel_confianza <- 0.95
z_critico <- qnorm(1 - (1 - nivel_confianza) / 2)

limite_varianza <- 1e-12
phi_maximo <- 0.995

# Se utiliza la meta contemporánea meta_t, porque la ecuación estática
# subyacente define e^L_t como promedio de la meta vigente en t,
# inflación t-1 y expectativas cortas t-1.
fecha_inicio_opcional <- as.Date("2011-01-01")


# 2. FUNCIONES GENERALES -------------------------------------------------------

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


# 3. CARGA Y PREPARACIÓN DE DATOS ---------------------------------------------

bd_original <- readr::read_csv(
  archivo_datos,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", ".", "null")
)

verificar_columnas(
  bd_original,
  columnas_requeridas
)

col_tc <- primera_columna_disponible(
  bd_original,
  candidatas_tc
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

bd_modelo <- bd |>
  mutate(
    exp24_lag1 = lag(exp_inf_24m, 1),
    exp12_lag1 = lag(exp_inf_12m, 1),
    infl_lag1 = lag(infl_gt, 1),

    d_exp24 =
      exp_inf_24m - exp24_lag1,

    # Componente de ajuste hacia la meta.
    gap_meta_largo =
      meta_inflacion - exp24_lag1,

    # Regresores de los parámetros de desanclaje.
    x_inflacion =
      infl_lag1 - meta_inflacion,

    x_expectativa_corta =
      exp12_lag1 - meta_inflacion,

    # Dependiente transformada equivalente a la ecuación de Strohsal.
    y_ajustada =
      d_exp24 - gap_meta_largo
  )

if (!is.na(col_tc)) {

  bd_modelo <- bd_modelo |>
    mutate(
      control_tc =
        lag(.data[[col_tc]], 1)
    )

} else {

  bd_modelo <- bd_modelo |>
    mutate(
      control_tc = NA_real_
    )

  warning(
    "No se encontró variable cambiaria. ",
    "Se estimará únicamente el núcleo Strohsal."
  )
}

bd_modelo <- bd_modelo |>
  filter(fecha >= fecha_inicio_opcional) |>
  filter(
    is.finite(d_exp24),
    is.finite(gap_meta_largo),
    is.finite(x_inflacion),
    is.finite(x_expectativa_corta),
    is.finite(y_ajustada)
  )

if (nrow(bd_modelo) < 60) {
  stop("La muestra efectiva tiene menos de 60 observaciones.")
}

saltos_meses <- diff(
  year(bd_modelo$fecha) * 12 +
    month(bd_modelo$fecha)
)

if (any(saltos_meses != 1)) {
  warning(
    "La muestra presenta saltos mensuales. ",
    "Revise observaciones faltantes."
  )
}

readr::write_csv(
  bd_modelo,
  file.path(
    carpeta_salida,
    "01_datos_utilizados_modelo_strohsal.csv"
  )
)


# 4. MODELO CONSTANTE DEL PAPER CON HAC ----------------------------------------

estimar_constante_hac <- function(
    datos,
    usar_control_tc = FALSE,
    nombre_especificacion
) {

  datos_estimacion <- datos

  if (usar_control_tc) {
    datos_estimacion <- datos_estimacion |>
      filter(is.finite(control_tc))
  }

  formula_cp <- if (usar_control_tc) {
    y_ajustada ~
      x_inflacion +
      x_expectativa_corta +
      control_tc
  } else {
    y_ajustada ~
      x_inflacion +
      x_expectativa_corta
  }

  modelo <- lm(
    formula_cp,
    data = datos_estimacion
  )

  # El paper emplea Newey-West con selección automática de bandwidth.
  bandwidth_nw <- sandwich::bwNeweyWest(
    modelo,
    prewhite = FALSE
  )

  lag_nw <- floor(bandwidth_nw)

  vcov_hac <- sandwich::NeweyWest(
    modelo,
    lag = lag_nw,
    prewhite = FALSE,
    adjust = TRUE
  )

  tabla_hac <- lmtest::coeftest(
    modelo,
    vcov. = vcov_hac
  )

  tabla_coeficientes <- tibble(
    especificacion = nombre_especificacion,
    parametro = rownames(tabla_hac),
    estimacion = tabla_hac[, "Estimate"],
    error_estandar_hac = tabla_hac[, "Std. Error"],
    estadistico_t = tabla_hac[, "t value"],
    p_valor = tabla_hac[, "Pr(>|t|)"],
    bandwidth_newey_west = bandwidth_nw,
    lag_newey_west = lag_nw
  )

  nombres_theta <- c(
    "x_inflacion",
    "x_expectativa_corta"
  )

  b_theta <- coef(modelo)[nombres_theta]

  V_theta <- vcov_hac[
    nombres_theta,
    nombres_theta,
    drop = FALSE
  ]

  W <- tryCatch(
    as.numeric(
      t(b_theta) %*%
        solve(V_theta) %*%
        b_theta
    ),
    error = function(e) NA_real_
  )

  p_W <- ifelse(
    is.finite(W),
    pchisq(W, df = 2, lower.tail = FALSE),
    NA_real_
  )

  theta1_cp <- unname(
    coef(modelo)["x_inflacion"]
  )

  theta2_cp <- unname(
    coef(modelo)["x_expectativa_corta"]
  )

  anclaje_cp <-
    1 - theta1_cp - theta2_cp

  gradiente_A <- matrix(
    c(-1, -1),
    ncol = 1
  )

  var_A <- tryCatch(
    as.numeric(
      t(gradiente_A) %*%
        V_theta %*%
        gradiente_A
    ),
    error = function(e) NA_real_
  )

  se_A <- ifelse(
    is.finite(var_A) && var_A >= 0,
    sqrt(var_A),
    NA_real_
  )

  z_perfecto <- ifelse(
    is.finite(se_A) && se_A > 0,
    (anclaje_cp - 1) / se_A,
    NA_real_
  )

  p_perfecto <- ifelse(
    is.finite(z_perfecto),
    2 * pnorm(abs(z_perfecto), lower.tail = FALSE),
    NA_real_
  )

  tabla_conjunta <- tibble(
    especificacion = nombre_especificacion,
    theta1 = theta1_cp,
    theta2 = theta2_cp,
    anclaje = anclaje_cp,
    anclaje_porcentaje = 100 * anclaje_cp,
    se_anclaje = se_A,
    anclaje_li_95 =
      anclaje_cp - z_critico * se_A,
    anclaje_ls_95 =
      anclaje_cp + z_critico * se_A,
    estadistico_Wald_theta1_theta2_cero = W,
    grados_libertad = 2,
    p_valor_Wald = p_W,
    p_valor_anclaje_igual_a_uno = p_perfecto,
    observaciones = nobs(modelo),
    r_squared = summary(modelo)$r.squared,
    adj_r_squared = summary(modelo)$adj.r.squared
  )

  list(
    modelo = modelo,
    datos = datos_estimacion,
    vcov_hac = vcov_hac,
    tabla_coeficientes = tabla_coeficientes,
    tabla_conjunta = tabla_conjunta
  )
}


cp_hac_nucleo <- estimar_constante_hac(
  datos = bd_modelo,
  usar_control_tc = FALSE,
  nombre_especificacion = "Núcleo Strohsal"
)

lista_cp_hac <- list(cp_hac_nucleo)

if (!is.na(col_tc)) {

  cp_hac_tc <- estimar_constante_hac(
    datos = bd_modelo,
    usar_control_tc = TRUE,
    nombre_especificacion = "Strohsal + tipo de cambio"
  )

  lista_cp_hac <- c(
    lista_cp_hac,
    list(cp_hac_tc)
  )
}

tabla_cp_coeficientes <- bind_rows(
  lapply(
    lista_cp_hac,
    function(x) x$tabla_coeficientes
  )
)

tabla_cp_conjunta <- bind_rows(
  lapply(
    lista_cp_hac,
    function(x) x$tabla_conjunta
  )
)

readr::write_csv(
  tabla_cp_coeficientes,
  file.path(
    carpeta_salida,
    "02_modelo_constante_coeficientes_HAC.csv"
  )
)

readr::write_csv(
  tabla_cp_conjunta,
  file.path(
    carpeta_salida,
    "03_modelo_constante_anclaje_y_prueba_conjunta.csv"
  )
)


# 5. CONSTRUCCIÓN DE LOS MODELOS EN ESPACIO DE ESTADOS ------------------------

crear_matrices_strohsal <- function(
    datos,
    usar_control_tc = FALSE,
    tvp = TRUE,
    estado_inicial
) {

  n <- nrow(datos)
  n_controles <- ifelse(usar_control_tc, 1L, 0L)

  # Estados:
  # 1 theta1_t
  # 2 theta2_t
  # 3 u_t
  # 4 constante_t = 1
  # 5 beta_tc, cuando aplica
  m <- 4L + n_controles

  Z_t <- array(
    0,
    dim = c(1, m, n)
  )

  Z_t[1, 1, ] <- datos$x_inflacion
  Z_t[1, 2, ] <- datos$x_expectativa_corta
  Z_t[1, 3, ] <- 1
  Z_t[1, 4, ] <- 0

  if (usar_control_tc) {
    Z_t[1, 5, ] <- datos$control_tc
  }

  T_mat <- diag(m)

  # Valores provisionales; se actualizan durante optimización.
  T_mat[3, 3] <- 0.50
  T_mat[3, 4] <- 0.00

  # Innovaciones de theta1, theta2 y u.
  R_mat <- matrix(
    0,
    nrow = m,
    ncol = 3
  )

  R_mat[1, 1] <- 1
  R_mat[2, 2] <- 1
  R_mat[3, 3] <- 1

  if (tvp) {

    Q_mat <- diag(
      rep(NA_real_, 3)
    )

  } else {

    # En el modelo constante theta1 y theta2 no evolucionan.
    Q_mat <- diag(
      c(0, 0, NA_real_)
    )
  }

  P1inf_diag <- c(
    1, # theta1
    1, # theta2
    1, # u
    0  # constante fija
  )

  if (usar_control_tc) {
    P1inf_diag <- c(
      P1inf_diag,
      1 # beta_tc fijo pero desconocido
    )
  }

  # SSModel reconoce SSMcustom como un término especial de la fórmula.
  # No debe escribirse KFAS::SSMcustom dentro de la fórmula, porque
  # model.frame lo interpreta como una variable de tipo lista.
  y_strohsal <- as.numeric(datos$y_ajustada)

  modelo <- KFAS::SSModel(
    y_strohsal ~ -1 +
      SSMcustom(
        Z = Z_t,
        T = T_mat,
        R = R_mat,
        Q = Q_mat,
        a1 = estado_inicial,
        P1 = matrix(
          0,
          nrow = m,
          ncol = m
        ),
        P1inf = diag(P1inf_diag),
        n = n
      ),
    H = matrix(
      0,
      nrow = 1,
      ncol = 1
    )
  )

  modelo
}


actualizar_modelo_cp <- function(par, modelo) {

  phi <- phi_maximo * tanh(par[1])
  alpha0 <- par[2]
  Q_nu <- exp(par[3])

  modelo$T[3, 3, 1] <- phi
  modelo$T[3, 4, 1] <- alpha0

  modelo$Q[1, 1, 1] <- 0
  modelo$Q[2, 2, 1] <- 0
  modelo$Q[3, 3, 1] <- Q_nu

  modelo
}


actualizar_modelo_tvp <- function(par, modelo) {

  phi <- phi_maximo * tanh(par[1])
  alpha0 <- par[2]

  Q_theta1 <- exp(par[3])
  Q_theta2 <- exp(par[4])
  Q_nu <- exp(par[5])

  modelo$T[3, 3, 1] <- phi
  modelo$T[3, 4, 1] <- alpha0

  modelo$Q[1, 1, 1] <- Q_theta1
  modelo$Q[2, 2, 1] <- Q_theta2
  modelo$Q[3, 3, 1] <- Q_nu

  modelo
}


# 6. INICIALIZACIÓN ------------------------------------------------------------

obtener_inicializacion <- function(
    datos,
    usar_control_tc = FALSE
) {

  formula_ols <- if (usar_control_tc) {
    y_ajustada ~
      x_inflacion +
      x_expectativa_corta +
      control_tc
  } else {
    y_ajustada ~
      x_inflacion +
      x_expectativa_corta
  }

  ols <- lm(
    formula_ols,
    data = datos
  )

  theta1_ini <- unname(
    coef(ols)["x_inflacion"]
  )

  theta2_ini <- unname(
    coef(ols)["x_expectativa_corta"]
  )

  beta_tc_ini <- if (usar_control_tc) {
    unname(coef(ols)["control_tc"])
  } else {
    numeric(0)
  }

  residuos <- residuals(ols)

  if (length(residuos) >= 3) {

    modelo_ar <- lm(
      residuos[-1] ~
        residuos[-length(residuos)]
    )

    alpha0_ini <- unname(
      coef(modelo_ar)[1]
    )

    phi_ini <- unname(
      coef(modelo_ar)[2]
    )

    phi_ini <- max(
      min(phi_ini, 0.90),
      -0.90
    )

    innovaciones_ar <- residuals(modelo_ar)

    q_nu_ini <- max(
      var(innovaciones_ar, na.rm = TRUE),
      1e-6
    )

  } else {

    alpha0_ini <- 0
    phi_ini <- 0.5
    q_nu_ini <- max(
      var(residuos, na.rm = TRUE),
      1e-6
    )
  }

  estado_inicial <- c(
    theta1_ini,
    theta2_ini,
    0,
    1,
    beta_tc_ini
  )

  list(
    ols = ols,
    estado_inicial = estado_inicial,
    theta1_ini = theta1_ini,
    theta2_ini = theta2_ini,
    beta_tc_ini = beta_tc_ini,
    alpha0_ini = alpha0_ini,
    phi_ini = phi_ini,
    q_nu_ini = q_nu_ini
  )
}


# 7. ESTIMACIÓN CON MÚLTIPLES INICIOS -----------------------------------------

estimar_modelo_cp_ssm <- function(
    modelo_sin_estimar,
    inicializacion,
    nombre_especificacion
) {

  raw_phi_ini <- atanh(
    inicializacion$phi_ini /
      phi_maximo
  )

  inicios <- list(
    c(
      raw_phi_ini,
      inicializacion$alpha0_ini,
      log(inicializacion$q_nu_ini)
    ),
    c(
      atanh(0.20 / phi_maximo),
      0,
      log(inicializacion$q_nu_ini)
    ),
    c(
      atanh(0.60 / phi_maximo),
      inicializacion$alpha0_ini,
      log(inicializacion$q_nu_ini * 0.5)
    ),
    c(
      atanh(0.85 / phi_maximo),
      inicializacion$alpha0_ini,
      log(inicializacion$q_nu_ini * 2)
    ),
    c(
      atanh(-0.20 / phi_maximo),
      0,
      log(inicializacion$q_nu_ini)
    )
  )

  ajustes <- lapply(
    seq_along(inicios),
    function(i) {

      ajuste <- tryCatch(
        KFAS::fitSSM(
          model = modelo_sin_estimar,
          inits = inicios[[i]],
          updatefn = actualizar_modelo_cp,
          method = "L-BFGS-B",
          lower = c(
            -4,
            -5,
            log(limite_varianza)
          ),
          upper = c(
            4,
            5,
            2
          ),
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
        as.numeric(logLik(ajuste$model)),
        error = function(e) NA_real_
      )

      list(
        intento = i,
        ajuste = ajuste,
        convergencia =
          ajuste$optim.out$convergence,
        mensaje = ifelse(
          is.null(ajuste$optim.out$message),
          "",
          as.character(ajuste$optim.out$message)
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
      "No fue posible estimar el modelo constante SSM: ",
      nombre_especificacion
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
      "Ningún intento CP convergió normalmente para ",
      nombre_especificacion,
      ". Se usará la mayor log-verosimilitud disponible."
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
          especificacion =
            nombre_especificacion,
          modelo = "Constante SSM",
          intento = x$intento,
          convergencia = x$convergencia,
          log_verosimilitud = x$logLik,
          mensaje = x$mensaje
        )
      }
    )
  )

  list(
    mejor = mejor,
    modelo = mejor$ajuste$model,
    tabla_intentos = tabla_intentos
  )
}


estimar_modelo_tvp_ssm <- function(
    modelo_sin_estimar,
    inicializacion,
    nombre_especificacion
) {

  raw_phi_ini <- atanh(
    inicializacion$phi_ini /
      phi_maximo
  )

  escalas_Q <- c(
    1e-8,
    1e-6,
    1e-5,
    1e-4,
    1e-3,
    1e-2
  )

  inicios <- lapply(
    escalas_Q,
    function(q_theta) {
      c(
        raw_phi_ini,
        inicializacion$alpha0_ini,
        log(q_theta),
        log(q_theta),
        log(inicializacion$q_nu_ini)
      )
    }
  )

  inicios <- c(
    inicios,
    list(
      c(
        atanh(0.20 / phi_maximo),
        0,
        log(1e-5),
        log(1e-4),
        log(inicializacion$q_nu_ini)
      ),
      c(
        atanh(0.70 / phi_maximo),
        inicializacion$alpha0_ini,
        log(1e-4),
        log(1e-5),
        log(inicializacion$q_nu_ini)
      ),
      c(
        atanh(0.90 / phi_maximo),
        inicializacion$alpha0_ini,
        log(1e-3),
        log(1e-3),
        log(inicializacion$q_nu_ini * 0.5)
      ),
      c(
        atanh(-0.20 / phi_maximo),
        0,
        log(1e-5),
        log(1e-5),
        log(inicializacion$q_nu_ini)
      )
    )
  )

  ajustes <- lapply(
    seq_along(inicios),
    function(i) {

      ajuste <- tryCatch(
        KFAS::fitSSM(
          model = modelo_sin_estimar,
          inits = inicios[[i]],
          updatefn = actualizar_modelo_tvp,
          method = "L-BFGS-B",
          lower = c(
            -4,
            -5,
            rep(log(limite_varianza), 3)
          ),
          upper = c(
            4,
            5,
            rep(2, 3)
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
        as.numeric(logLik(ajuste$model)),
        error = function(e) NA_real_
      )

      list(
        intento = i,
        ajuste = ajuste,
        convergencia =
          ajuste$optim.out$convergence,
        mensaje = ifelse(
          is.null(ajuste$optim.out$message),
          "",
          as.character(ajuste$optim.out$message)
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
      "No fue posible estimar el modelo TVP: ",
      nombre_especificacion
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
      "Ningún intento TVP convergió normalmente para ",
      nombre_especificacion,
      ". Se usará la mayor log-verosimilitud disponible."
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
          especificacion =
            nombre_especificacion,
          modelo = "TVP Strohsal",
          intento = x$intento,
          convergencia = x$convergencia,
          log_verosimilitud = x$logLik,
          mensaje = x$mensaje
        )
      }
    )
  )

  list(
    mejor = mejor,
    modelo = mejor$ajuste$model,
    tabla_intentos = tabla_intentos
  )
}


# 8. FUNCIÓN INTEGRAL POR ESPECIFICACIÓN --------------------------------------

estimar_especificacion_strohsal <- function(
    datos,
    usar_control_tc = FALSE,
    nombre_especificacion,
    etiqueta_archivo
) {

  datos_estimacion <- datos

  if (usar_control_tc) {

    datos_estimacion <- datos_estimacion |>
      filter(is.finite(control_tc))
  }

  inicializacion <- obtener_inicializacion(
    datos = datos_estimacion,
    usar_control_tc = usar_control_tc
  )

  modelo_cp_sin_estimar <- crear_matrices_strohsal(
    datos = datos_estimacion,
    usar_control_tc = usar_control_tc,
    tvp = FALSE,
    estado_inicial =
      inicializacion$estado_inicial
  )

  modelo_tvp_sin_estimar <- crear_matrices_strohsal(
    datos = datos_estimacion,
    usar_control_tc = usar_control_tc,
    tvp = TRUE,
    estado_inicial =
      inicializacion$estado_inicial
  )

  ajuste_cp <- estimar_modelo_cp_ssm(
    modelo_sin_estimar =
      modelo_cp_sin_estimar,
    inicializacion = inicializacion,
    nombre_especificacion =
      nombre_especificacion
  )

  ajuste_tvp <- estimar_modelo_tvp_ssm(
    modelo_sin_estimar =
      modelo_tvp_sin_estimar,
    inicializacion = inicializacion,
    nombre_especificacion =
      nombre_especificacion
  )

  kfs_cp <- KFAS::KFS(
    ajuste_cp$modelo,
    filtering = c("state", "mean"),
    smoothing = c(
      "state",
      "mean",
      "disturbance"
    )
  )

  kfs_tvp <- KFAS::KFS(
    ajuste_tvp$modelo,
    filtering = c("state", "mean"),
    smoothing = c(
      "state",
      "mean",
      "disturbance"
    )
  )

  # Trayectorias TVP.
  theta1 <- as.numeric(
    kfs_tvp$alphahat[, 1]
  )

  theta2 <- as.numeric(
    kfs_tvp$alphahat[, 2]
  )

  var_theta1 <- pmax(
    as.numeric(kfs_tvp$V[1, 1, ]),
    0
  )

  var_theta2 <- pmax(
    as.numeric(kfs_tvp$V[2, 2, ]),
    0
  )

  cov_theta12 <- as.numeric(
    kfs_tvp$V[1, 2, ]
  )

  se_theta1 <- sqrt(var_theta1)
  se_theta2 <- sqrt(var_theta2)

  anclaje <- 1 - theta1 - theta2

  var_anclaje <-
    var_theta1 +
    var_theta2 +
    2 * cov_theta12

  var_anclaje <- pmax(
    var_anclaje,
    0
  )

  se_anclaje <- sqrt(
    var_anclaje
  )

  resultados_tvp <- datos_estimacion |>
    transmute(
      fecha,
      exp_inf_24m,
      exp_inf_12m,
      infl_gt,
      meta_inflacion,
      d_exp24,
      gap_meta_largo,
      x_inflacion,
      x_expectativa_corta,
      control_tc,

      theta1_inflacion = theta1,
      theta1_se = se_theta1,
      theta1_li_95 =
        theta1 - z_critico * se_theta1,
      theta1_ls_95 =
        theta1 + z_critico * se_theta1,

      theta2_expectativa_corta = theta2,
      theta2_se = se_theta2,
      theta2_li_95 =
        theta2 - z_critico * se_theta2,
      theta2_ls_95 =
        theta2 + z_critico * se_theta2,

      parametro_anclaje = anclaje,
      anclaje_porcentaje = 100 * anclaje,
      anclaje_se = se_anclaje,
      anclaje_li_95 =
        anclaje - z_critico * se_anclaje,
      anclaje_ls_95 =
        anclaje + z_critico * se_anclaje,

      # Prueba puntual H0: A_t = 1.
      z_anclaje_igual_uno = ifelse(
        se_anclaje > 0,
        (anclaje - 1) / se_anclaje,
        NA_real_
      ),

      p_anclaje_igual_uno =
        2 * pnorm(
          abs(z_anclaje_igual_uno),
          lower.tail = FALSE
        ),

      anclaje_significativamente_menor_uno =
        anclaje_ls_95 < 1,

      anclaje_significativamente_mayor_cero =
        anclaje_li_95 > 0,

      theta1_significativamente_positivo =
        theta1_li_95 > 0,

      theta2_significativamente_positivo =
        theta2_li_95 > 0,

      fuera_rango_0_1 =
        parametro_anclaje < 0 |
        parametro_anclaje > 1,

      residuo_estado_u =
        as.numeric(kfs_tvp$alphahat[, 3])
    )

  # Estados constantes del comparador.
  theta1_cp <- mean(
    as.numeric(kfs_cp$alphahat[, 1]),
    na.rm = TRUE
  )

  theta2_cp <- mean(
    as.numeric(kfs_cp$alphahat[, 2]),
    na.rm = TRUE
  )

  beta_tc_cp <- if (usar_control_tc) {
    mean(
      as.numeric(kfs_cp$alphahat[, 5]),
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  beta_tc_tvp <- if (usar_control_tc) {
    mean(
      as.numeric(kfs_tvp$alphahat[, 5]),
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  phi_cp <- as.numeric(
    ajuste_cp$modelo$T[3, 3, 1]
  )

  alpha0_cp <- as.numeric(
    ajuste_cp$modelo$T[3, 4, 1]
  )

  phi_tvp <- as.numeric(
    ajuste_tvp$modelo$T[3, 3, 1]
  )

  alpha0_tvp <- as.numeric(
    ajuste_tvp$modelo$T[3, 4, 1]
  )

  Q_cp <- diag(
    ajuste_cp$modelo$Q[, , 1]
  )

  Q_tvp <- diag(
    ajuste_tvp$modelo$Q[, , 1]
  )

  tabla_hiperparametros <- bind_rows(
    tibble(
      especificacion =
        nombre_especificacion,
      modelo = "Constante SSM",
      theta1_promedio = theta1_cp,
      theta2_promedio = theta2_cp,
      anclaje_promedio =
        1 - theta1_cp - theta2_cp,
      beta_control_tc = beta_tc_cp,
      alpha0 = alpha0_cp,
      alpha1 = phi_cp,
      Q_theta1 = Q_cp[1],
      Q_theta2 = Q_cp[2],
      Q_nu = Q_cp[3],
      H = 0,
      convergencia =
        ajuste_cp$mejor$convergencia,
      log_verosimilitud =
        as.numeric(
          logLik(ajuste_cp$modelo)
        )
    ),

    tibble(
      especificacion =
        nombre_especificacion,
      modelo = "TVP Strohsal",
      theta1_promedio =
        mean(theta1, na.rm = TRUE),
      theta2_promedio =
        mean(theta2, na.rm = TRUE),
      anclaje_promedio =
        mean(anclaje, na.rm = TRUE),
      beta_control_tc = beta_tc_tvp,
      alpha0 = alpha0_tvp,
      alpha1 = phi_tvp,
      Q_theta1 = Q_tvp[1],
      Q_theta2 = Q_tvp[2],
      Q_nu = Q_tvp[3],
      H = 0,
      convergencia =
        ajuste_tvp$mejor$convergencia,
      log_verosimilitud =
        as.numeric(
          logLik(ajuste_tvp$modelo)
        )
    )
  )

  # Comparación CP vs TVP con la misma estructura de error AR(1).
  ll_cp <- as.numeric(
    logLik(ajuste_cp$modelo)
  )

  ll_tvp <- as.numeric(
    logLik(ajuste_tvp$modelo)
  )

  n_obs <- nrow(datos_estimacion)
  n_controles <- ifelse(
    usar_control_tc,
    1,
    0
  )

  # Conteo:
  # CP: theta1, theta2, controles, alpha0, alpha1, Q_nu.
  # TVP añade Q_theta1 y Q_theta2.
  k_cp <- 5 + n_controles
  k_tvp <- k_cp + 2

  AIC_cp <- -2 * ll_cp + 2 * k_cp
  AIC_tvp <- -2 * ll_tvp + 2 * k_tvp

  BIC_cp <- -2 * ll_cp + log(n_obs) * k_cp
  BIC_tvp <- -2 * ll_tvp + log(n_obs) * k_tvp

  AICc_cp <- calcular_AICc(
    AIC_cp,
    k_cp,
    n_obs
  )

  AICc_tvp <- calcular_AICc(
    AIC_tvp,
    k_tvp,
    n_obs
  )

  tabla_comparacion <- tibble(
    especificacion =
      nombre_especificacion,
    modelo = c(
      "Constante SSM",
      "TVP Strohsal"
    ),
    observaciones = n_obs,
    parametros_contados = c(
      k_cp,
      k_tvp
    ),
    log_verosimilitud = c(
      ll_cp,
      ll_tvp
    ),
    AIC = c(
      AIC_cp,
      AIC_tvp
    ),
    AICc = c(
      AICc_cp,
      AICc_tvp
    ),
    BIC = c(
      BIC_cp,
      BIC_tvp
    )
  ) |>
    mutate(
      delta_AIC =
        AIC - min(AIC),
      peso_Akaike =
        exp(-0.5 * delta_AIC) /
        sum(exp(-0.5 * delta_AIC)),
      delta_BIC =
        BIC - min(BIC),
      peso_BIC_aproximado =
        exp(-0.5 * delta_BIC) /
        sum(exp(-0.5 * delta_BIC))
    )

  LR <- max(
    0,
    2 * (ll_tvp - ll_cp)
  )

  p_LR_chi2_2 <- pchisq(
    LR,
    df = 2,
    lower.tail = FALSE
  )

  # Mezcla referencial para dos componentes de varianza en frontera:
  # 0.25 chi2_0 + 0.50 chi2_1 + 0.25 chi2_2.
  # No sustituye un bootstrap paramétrico.
  p_LR_mezcla_referencial <- ifelse(
    LR > 0,
    0.50 * pchisq(
      LR,
      df = 1,
      lower.tail = FALSE
    ) +
      0.25 * pchisq(
        LR,
        df = 2,
        lower.tail = FALSE
      ),
    1
  )

  tabla_LR <- tibble(
    especificacion =
      nombre_especificacion,
    hipotesis_nula =
      "Q_theta1 = Q_theta2 = 0",
    estadistico_LR = LR,
    p_valor_chi2_2_referencial =
      p_LR_chi2_2,
    p_valor_mezcla_frontera_referencial =
      p_LR_mezcla_referencial,
    advertencia = paste(
      "Dos varianzas están en la frontera bajo H0.",
      "Los p-valores son referenciales;",
      "un bootstrap paramétrico es preferible para inferencia definitiva."
    )
  )

  # Residuos estandarizados recursivos.
  residuos_cp <- tryCatch(
    as.numeric(
      rstandard(
        kfs_cp,
        type = "recursive"
      )
    ),
    error = function(e) {
      rep(NA_real_, n_obs)
    }
  )

  residuos_tvp <- tryCatch(
    as.numeric(
      rstandard(
        kfs_tvp,
        type = "recursive"
      )
    ),
    error = function(e) {
      rep(NA_real_, n_obs)
    }
  )

  calcular_diagnosticos <- function(
      residuos,
      nombre_modelo
  ) {

    residuos <- residuos[
      is.finite(residuos)
    ]

    rezagos <- c(
      1,
      2,
      5,
      10,
      12,
      18,
      24
    )

    rezagos <- rezagos[
      rezagos < length(residuos)
    ]

    bind_rows(
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
            especificacion =
              nombre_especificacion,
            modelo = nombre_modelo,
            diagnostico = "Ljung-Box",
            rezago = L,
            estadistico =
              unname(prueba$statistic),
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
  }

  diagnosticos <- bind_rows(
    calcular_diagnosticos(
      residuos_cp,
      "Constante SSM"
    ),
    calcular_diagnosticos(
      residuos_tvp,
      "TVP Strohsal"
    )
  )

  # Resumen de trayectoria.
  indice_max_A <- which.max(
    resultados_tvp$parametro_anclaje
  )

  indice_min_A <- which.min(
    resultados_tvp$parametro_anclaje
  )

  resumen_anclaje <- tibble(
    especificacion =
      nombre_especificacion,
    fecha_inicial =
      first(resultados_tvp$fecha),
    fecha_final =
      last(resultados_tvp$fecha),
    observaciones =
      nrow(resultados_tvp),

    anclaje_promedio =
      mean(
        resultados_tvp$parametro_anclaje,
        na.rm = TRUE
      ),

    anclaje_promedio_porcentaje =
      100 * mean(
        resultados_tvp$parametro_anclaje,
        na.rm = TRUE
      ),

    anclaje_inicial =
      first(
        resultados_tvp$parametro_anclaje
      ),

    anclaje_final =
      last(
        resultados_tvp$parametro_anclaje
      ),

    anclaje_minimo =
      resultados_tvp$parametro_anclaje[
        indice_min_A
      ],

    fecha_anclaje_minimo =
      resultados_tvp$fecha[
        indice_min_A
      ],

    anclaje_maximo =
      resultados_tvp$parametro_anclaje[
        indice_max_A
      ],

    fecha_anclaje_maximo =
      resultados_tvp$fecha[
        indice_max_A
      ],

    porcentaje_meses_A_significativamente_menor_1 =
      100 * mean(
        resultados_tvp$
          anclaje_significativamente_menor_uno,
        na.rm = TRUE
      ),

    porcentaje_meses_A_significativamente_mayor_0 =
      100 * mean(
        resultados_tvp$
          anclaje_significativamente_mayor_cero,
        na.rm = TRUE
      ),

    porcentaje_meses_theta1_significativo_positivo =
      100 * mean(
        resultados_tvp$
          theta1_significativamente_positivo,
        na.rm = TRUE
      ),

    porcentaje_meses_theta2_significativo_positivo =
      100 * mean(
        resultados_tvp$
          theta2_significativamente_positivo,
        na.rm = TRUE
      ),

    porcentaje_meses_A_fuera_0_1 =
      100 * mean(
        resultados_tvp$fuera_rango_0_1,
        na.rm = TRUE
      )
  )

  # Exportaciones por especificación.
  readr::write_csv(
    resultados_tvp,
    file.path(
      carpeta_salida,
      paste0(
        "04_trayectoria_tvp_",
        etiqueta_archivo,
        ".csv"
      )
    )
  )

  saveRDS(
    ajuste_cp$modelo,
    file.path(
      carpeta_salida,
      paste0(
        "20_modelo_constante_ssm_",
        etiqueta_archivo,
        ".rds"
      )
    )
  )

  saveRDS(
    ajuste_tvp$modelo,
    file.path(
      carpeta_salida,
      paste0(
        "21_modelo_tvp_strohsal_",
        etiqueta_archivo,
        ".rds"
      )
    )
  )

  # Gráfica del parámetro de anclaje.
  p_anclaje <- ggplot(
    resultados_tvp,
    aes(
      x = fecha,
      y = parametro_anclaje
    )
  ) +
    geom_ribbon(
      aes(
        ymin = anclaje_li_95,
        ymax = anclaje_ls_95
      ),
      alpha = 0.18
    ) +
    geom_line(
      linewidth = 0.9
    ) +
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
        paste(
          "Parámetro dinámico de anclaje:",
          nombre_especificacion
        ),
      subtitle =
        expression(
          A[t] == 1 - theta[1*t] - theta[2*t]
        ),
      x = NULL,
      y = "Parámetro de anclaje",
      caption = paste(
        "La línea en 1 representa anclaje perfecto bajo el modelo.",
        "Bandas del 95%."
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title.position = "plot",
      panel.grid.minor = element_blank()
    )

  ggsave(
    file.path(
      carpeta_salida,
      paste0(
        "10_parametro_anclaje_",
        etiqueta_archivo,
        ".png"
      )
    ),
    plot = p_anclaje,
    width = 10,
    height = 6,
    dpi = 300
  )

  # Fuentes de desanclaje.
  datos_theta <- resultados_tvp |>
    select(
      fecha,
      `Inflación pasada: theta1` =
        theta1_inflacion,
      `Expectativa de corto plazo: theta2` =
        theta2_expectativa_corta
    ) |>
    pivot_longer(
      cols = -fecha,
      names_to = "fuente",
      values_to = "coeficiente"
    )

  p_fuentes <- ggplot(
    datos_theta,
    aes(
      x = fecha,
      y = coeficiente,
      linetype = fuente
    )
  ) +
    geom_line(
      linewidth = 0.85
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    labs(
      title =
        paste(
          "Fuentes dinámicas de desanclaje:",
          nombre_especificacion
        ),
      x = NULL,
      y = "Coeficiente variante",
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
      paste0(
        "11_fuentes_desanclaje_",
        etiqueta_archivo,
        ".png"
      )
    ),
    plot = p_fuentes,
    width = 10,
    height = 6,
    dpi = 300
  )

  # ACF de residuos TVP.
  residuos_tvp_validos <- residuos_tvp[
    is.finite(residuos_tvp)
  ]

  png(
    filename = file.path(
      carpeta_salida,
      paste0(
        "12_acf_residuos_tvp_",
        etiqueta_archivo,
        ".png"
      )
    ),
    width = 1800,
    height = 1100,
    res = 180
  )

  acf(
    residuos_tvp_validos,
    lag.max = 24,
    main = paste(
      "ACF residuos TVP Strohsal:",
      nombre_especificacion
    )
  )

  dev.off()

  list(
    nombre_especificacion =
      nombre_especificacion,
    etiqueta_archivo =
      etiqueta_archivo,
    datos =
      datos_estimacion,
    ajuste_cp =
      ajuste_cp,
    ajuste_tvp =
      ajuste_tvp,
    kfs_cp =
      kfs_cp,
    kfs_tvp =
      kfs_tvp,
    tabla_intentos =
      bind_rows(
        ajuste_cp$tabla_intentos,
        ajuste_tvp$tabla_intentos
      ),
    tabla_hiperparametros =
      tabla_hiperparametros,
    tabla_comparacion =
      tabla_comparacion,
    tabla_LR =
      tabla_LR,
    diagnosticos =
      diagnosticos,
    resultados_tvp =
      resultados_tvp,
    resumen_anclaje =
      resumen_anclaje
  )
}


# 9. ESTIMACIÓN DE LAS ESPECIFICACIONES ---------------------------------------

resultado_nucleo <- estimar_especificacion_strohsal(
  datos = bd_modelo,
  usar_control_tc = FALSE,
  nombre_especificacion =
    "Núcleo Strohsal",
  etiqueta_archivo =
    "nucleo"
)

resultados_especificaciones <- list(
  resultado_nucleo
)

if (!is.na(col_tc)) {

  resultado_tc <- estimar_especificacion_strohsal(
    datos = bd_modelo,
    usar_control_tc = TRUE,
    nombre_especificacion =
      "Strohsal + tipo de cambio",
    etiqueta_archivo =
      "control_tc"
  )

  resultados_especificaciones <- c(
    resultados_especificaciones,
    list(resultado_tc)
  )
}


# 10. CONSOLIDACIÓN ------------------------------------------------------------

tabla_intentos <- bind_rows(
  lapply(
    resultados_especificaciones,
    function(x) x$tabla_intentos
  )
)

tabla_hiperparametros <- bind_rows(
  lapply(
    resultados_especificaciones,
    function(x) x$tabla_hiperparametros
  )
)

tabla_comparacion <- bind_rows(
  lapply(
    resultados_especificaciones,
    function(x) x$tabla_comparacion
  )
)

tabla_LR <- bind_rows(
  lapply(
    resultados_especificaciones,
    function(x) x$tabla_LR
  )
)

tabla_diagnosticos <- bind_rows(
  lapply(
    resultados_especificaciones,
    function(x) x$diagnosticos
  )
)

tabla_resumen_anclaje <- bind_rows(
  lapply(
    resultados_especificaciones,
    function(x) x$resumen_anclaje
  )
)

readr::write_csv(
  tabla_intentos,
  file.path(
    carpeta_salida,
    "05_intentos_estimacion.csv"
  )
)

readr::write_csv(
  tabla_hiperparametros,
  file.path(
    carpeta_salida,
    "06_hiperparametros_modelos_ssm.csv"
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
    "08_prueba_LR_variacion_theta.csv"
  )
)

readr::write_csv(
  tabla_diagnosticos,
  file.path(
    carpeta_salida,
    "09_diagnosticos_residuales.csv"
  )
)

readr::write_csv(
  tabla_resumen_anclaje,
  file.path(
    carpeta_salida,
    "13_resumen_parametro_anclaje.csv"
  )
)


# 11. RESUMEN EJECUTIVO --------------------------------------------------------

lineas_resumen <- c(
  "",
  "======================================================================",
  "MODELO DE STROHSAL APLICADO A GUATEMALA",
  "======================================================================",
  paste0("Archivo: ", archivo_datos),
  paste0(
    "Muestra base efectiva: ",
    format(min(bd_modelo$fecha), "%Y-%m"),
    " a ",
    format(max(bd_modelo$fecha), "%Y-%m")
  ),
  paste0(
    "Observaciones base: ",
    nrow(bd_modelo)
  ),
  paste0(
    "Variable cambiaria disponible: ",
    ifelse(
      is.na(col_tc),
      "No",
      col_tc
    )
  ),
  "",
  "MODELO",
  paste(
    "theta1_t y theta2_t evolucionan como paseos aleatorios;"
  ),
  paste(
    "u_t sigue un AR(1) con intercepto; H = 0;"
  ),
  paste(
    "A_t = 1 - theta1_t - theta2_t."
  ),
  ""
)

for (resultado in resultados_especificaciones) {

  nombre <- resultado$nombre_especificacion

  resumen <- resultado$resumen_anclaje

  hiper_tvp <- resultado$tabla_hiperparametros |>
    filter(modelo == "TVP Strohsal")

  comparacion_tvp <- resultado$tabla_comparacion |>
    filter(modelo == "TVP Strohsal")

  comparacion_cp <- resultado$tabla_comparacion |>
    filter(modelo == "Constante SSM")

  prueba_lr <- resultado$tabla_LR

  diagnosticos_tvp <- resultado$diagnosticos |>
    filter(modelo == "TVP Strohsal")

  p_lb_10 <- diagnosticos_tvp$p_valor[
    diagnosticos_tvp$rezago == 10
  ]

  lineas_resumen <- c(
    lineas_resumen,
    "----------------------------------------------------------------------",
    nombre,
    "----------------------------------------------------------------------",
    paste0(
      "Convergencia TVP: ",
      hiper_tvp$convergencia
    ),
    paste0(
      "LogLik CP / TVP: ",
      formato_num(
        comparacion_cp$log_verosimilitud
      ),
      " / ",
      formato_num(
        comparacion_tvp$log_verosimilitud
      )
    ),
    paste0(
      "AIC CP / TVP: ",
      formato_num(
        comparacion_cp$AIC
      ),
      " / ",
      formato_num(
        comparacion_tvp$AIC
      )
    ),
    paste0(
      "BIC CP / TVP: ",
      formato_num(
        comparacion_cp$BIC
      ),
      " / ",
      formato_num(
        comparacion_tvp$BIC
      )
    ),
    paste0(
      "LR variación theta: ",
      formato_num(
        prueba_lr$estadistico_LR
      )
    ),
    paste0(
      "p mezcla frontera referencial: ",
      formato_cientifico(
        prueba_lr$
          p_valor_mezcla_frontera_referencial
      )
    ),
    paste0(
      "Q_theta1: ",
      formato_cientifico(
        hiper_tvp$Q_theta1
      )
    ),
    paste0(
      "Q_theta2: ",
      formato_cientifico(
        hiper_tvp$Q_theta2
      )
    ),
    paste0(
      "Q_nu: ",
      formato_cientifico(
        hiper_tvp$Q_nu
      )
    ),
    paste0(
      "alpha0: ",
      formato_num(
        hiper_tvp$alpha0
      )
    ),
    paste0(
      "alpha1: ",
      formato_num(
        hiper_tvp$alpha1
      )
    ),
    paste0(
      "Anclaje promedio: ",
      formato_num(
        resumen$anclaje_promedio,
        4
      ),
      " (",
      formato_num(
        resumen$anclaje_promedio_porcentaje,
        2
      ),
      "%)"
    ),
    paste0(
      "Anclaje inicial/final: ",
      formato_num(
        resumen$anclaje_inicial,
        4
      ),
      " / ",
      formato_num(
        resumen$anclaje_final,
        4
      )
    ),
    paste0(
      "Mínimo: ",
      formato_num(
        resumen$anclaje_minimo,
        4
      ),
      " en ",
      format(
        resumen$fecha_anclaje_minimo,
        "%Y-%m"
      )
    ),
    paste0(
      "% meses A_t significativamente menor que 1: ",
      formato_num(
        resumen$
          porcentaje_meses_A_significativamente_menor_1,
        2
      ),
      "%"
    ),
    paste0(
      "% meses theta1 significativamente positivo: ",
      formato_num(
        resumen$
          porcentaje_meses_theta1_significativo_positivo,
        2
      ),
      "%"
    ),
    paste0(
      "% meses theta2 significativamente positivo: ",
      formato_num(
        resumen$
          porcentaje_meses_theta2_significativo_positivo,
        2
      ),
      "%"
    ),
    paste0(
      "% meses A_t fuera de [0,1]: ",
      formato_num(
        resumen$
          porcentaje_meses_A_fuera_0_1,
        2
      ),
      "%"
    ),
    paste0(
      "Ljung-Box TVP, lag 10, p: ",
      formato_cientifico(p_lb_10)
    ),
    ""
  )
}

lineas_resumen <- c(
  lineas_resumen,
  "ADVERTENCIAS DE INTERPRETACIÓN",
  paste(
    "1. La estructura econométrica replica el núcleo de Strohsal,"
  ),
  paste(
    "   pero 24 meses no equivale a las expectativas forward de 10 años."
  ),
  paste(
    "2. La meta de Guatemala cambia al inicio de la muestra;"
  ),
  paste(
    "   se usa la meta contemporánea en la ecuación estática."
  ),
  paste(
    "3. La prueba LR tiene dos varianzas en frontera;"
  ),
  paste(
    "   el p-valor reportado es referencial y puede requerir bootstrap."
  ),
  paste(
    "4. Valores de A_t fuera de [0,1] no se recortan artificialmente."
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
