# ==============================================================================
# OBJETIVO 2 - DEMERTZIS TVP-VAR RESTRINGIDO (M3) - VERSION FINAL DEPURADA
# Guatemala: expectativas de inflación a 24 meses
# ==============================================================================
#
# MODELO PRINCIPAL
#
#   Inflación:
#     pi_t = a0 + a*pi_(t-1) + b*e_(t-1) + eps_pi,t
#
#   Expectativas a 24 meses:
#     e_t = c0_t + c_t*pi_(t-1) + d*e_(t-1) + eps_e,t
#
#   Estados variantes:
#     c0_t = c0_(t-1) + eta_c0,t
#     c_t  = c_(t-1)  + eta_c,t
#
#   Parámetros constantes:
#     a0, a, b, d
#
# PROXY DE CREDIBILIDAD
#
#     lambda_t = 1 - c_t/(1-d)
#
# ANCLA INFLACIONARIA IMPLÍCITA
#
#     pi_estrella_t = c0_t/(1-d-c_t)
#
# NOTA METODOLÓGICA
# -----------------
# El modelo es una adaptación restringida del TVP-VAR de Demertzis, Marcellino
# y Viegi. El modelo original permite variación temporal en todos los
# coeficientes; aquí la variación se concentra en c0_t y c_t, de acuerdo con la
# especificación M3 seleccionada para Guatemala.
#
# Este script estima:
#   1. M0: VAR(1) con parámetros constantes, como benchmark.
#   2. M3: TVP-VAR restringido, modelo principal.
#
# No vuelve a estimar M1, M2 ni M4. Las pruebas de especificaciones alternativas,
# pulsos, leave-one-out y medidas alternativas de inflación deben conservarse
# como ejercicios de robustez separados, no como parte del script central.
#
# SALIDA
# ------
# Carpeta: resultados_objetivo2_TVP_VAR
#
# Principales archivos:
#   00_resumen_ejecutivo.txt
#   01_datos_utilizados.csv
#   02_coeficientes_OLS_inicializacion.csv
#   03_intentos_estimacion.csv
#   04_hiperparametros_M0_M3.csv
#   05_comparacion_M0_M3.csv
#   06_trayectoria_M3.csv
#   07_resumen_M3.csv
#   08_resumen_subperiodos_M3.csv
#   09_diagnosticos_residuales_M3.csv
#   10_estabilidad_local_M3.csv
#   11_residuos_M3.csv
#   figura_01_c0_t.png
#   figura_02_c_t.png
#   figura_03_lambda_t.png
#   figura_04_ancla_implicita.png
#   figura_05_brecha_ancla_meta.png
#   figura_06_ACF_residuos_inflacion.png
#   figura_07_ACF_residuos_expectativa.png
#   modelo_M0.rds
#   modelo_M3.rds
#   kfs_M3.rds
# ==============================================================================

setwd("/Users/paulogarridogrijalva/Documents/PES/Seminario/Inv_expectativas/Data y resultados")
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
getwd()

# 1. CONFIGURACIÓN -------------------------------------------------------------

archivo_entrada <- "bd_sample_modelo.csv"

if (!file.exists(archivo_entrada)) {
  stop(
    "No se encontró 'bd_sample_modelo.csv' en el directorio actual: ",
    getwd()
  )
}

carpeta_salida <- file.path(
  getwd(),
  "resultados_objetivo2_TVP_VAR"
)

dir.create(
  carpeta_salida,
  recursive = TRUE,
  showWarnings = FALSE
)


fecha_fin_muestra <- as.Date("2026-06-01")

columnas_requeridas <- c(
  "fecha",
  "infl_gt",
  "exp_inf_24m"
)

nivel_confianza <- 0.95
z_critico <- qnorm(1 - (1 - nivel_confianza) / 2)

limite_varianza <- 1e-12
tolerancia_denominador <- 0.05

# Diagnósticos principales.
rezagos_ljung_box <- c(6, 12, 24)

# Número máximo de iteraciones del optimizador.
maxit_M0 <- 20000
maxit_M3 <- 30000


# 2. FUNCIONES AUXILIARES ------------------------------------------------------

verificar_columnas <- function(datos, columnas) {

  faltan <- setdiff(columnas, names(datos))

  if (length(faltan) > 0) {
    stop(
      "Faltan las siguientes variables requeridas: ",
      paste(faltan, collapse = ", ")
    )
  }
}


jarque_bera_manual <- function(x) {

  x <- x[is.finite(x)]
  n <- length(x)

  if (n < 8 || !is.finite(stats::sd(x)) || stats::sd(x) == 0) {
    return(
      tibble::tibble(
        estadistico = NA_real_,
        p_valor = NA_real_
      )
    )
  }

  z <- (x - mean(x)) / stats::sd(x)
  asimetria <- mean(z^3)
  curtosis <- mean(z^4)

  JB <- n / 6 * (
    asimetria^2 +
      (curtosis - 3)^2 / 4
  )

  tibble::tibble(
    estadistico = JB,
    p_valor = stats::pchisq(
      JB,
      df = 2,
      lower.tail = FALSE
    )
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
    fecha = as.Date(fecha)
  ) |>
  dplyr::arrange(fecha) |>
  dplyr::filter(
    fecha <= fecha_fin_muestra
  )

if (anyNA(bd$fecha)) {
  stop("Existen fechas que no pudieron convertirse al formato Date.")
}

if (anyDuplicated(bd$fecha) > 0) {
  stop("Existen fechas duplicadas en la base.")
}

# La meta NO entra en la estimación del TVP-VAR.
# Se incorpora únicamente para contrastar posteriormente el ancla implícita.
if (!"meta_inflacion" %in% names(bd)) {

  bd <- bd |>
    dplyr::mutate(
      meta_inflacion = dplyr::case_when(
        fecha < as.Date("2012-01-01") ~ 5.0,
        fecha < as.Date("2013-01-01") ~ 4.5,
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
  stop("La muestra efectiva contiene menos de 80 observaciones.")
}

# Verificación de continuidad mensual.
indice_mensual <- lubridate::year(datos_niveles$fecha) * 12 +
  lubridate::month(datos_niveles$fecha)

if (any(diff(indice_mensual) != 1)) {
  stop(
    "La muestra presenta saltos mensuales después del filtrado. ",
    "Revise observaciones faltantes en inflación o expectativas."
  )
}

# VAR(1): la primera observación se pierde por el rezago.
datos_modelo <- datos_niveles |>
  dplyr::mutate(
    inflacion_lag1 = dplyr::lag(inflacion),
    expectativa_lag1 = dplyr::lag(expectativa)
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
    "01_datos_utilizados.csv"
  )
)


# 4. VAR(1) OLS PARA INICIALIZACIÓN -------------------------------------------

inicializar_desde_OLS <- function(datos) {

  eq_pi <- stats::lm(
    inflacion ~ inflacion_lag1 + expectativa_lag1,
    data = datos
  )

  eq_e <- stats::lm(
    expectativa ~ inflacion_lag1 + expectativa_lag1,
    data = datos
  )

  beta_pi <- stats::coef(eq_pi)
  beta_e <- stats::coef(eq_e)

  estado_inicial <- c(
    a0 = unname(beta_pi["(Intercept)"]),
    a = unname(beta_pi["inflacion_lag1"]),
    b = unname(beta_pi["expectativa_lag1"]),
    c0 = unname(beta_e["(Intercept)"]),
    c = unname(beta_e["inflacion_lag1"]),
    d = unname(beta_e["expectativa_lag1"])
  )

  if (any(!is.finite(estado_inicial))) {
    stop(
      "La inicialización OLS produjo coeficientes no finitos. ",
      "Revise la base y la posible colinealidad."
    )
  }
#significa parallel maximum: compara elemento por elemento y se queda con el mayor valor.
  #<- H es la matriz de varianzas de los errores de observación.
  H_inicial <- pmax(
    c(
      stats::var(stats::residuals(eq_pi)),
      stats::var(stats::residuals(eq_e))
    ),
    1e-8
  )

  list(
    estado_inicial = estado_inicial,
    H_inicial = H_inicial,
    eq_pi = eq_pi,
    eq_e = eq_e
  )
}

init <- inicializar_desde_OLS(datos_modelo)

estado_inicial <- init$estado_inicial
H_inicial <- init$H_inicial

tabla_OLS <- dplyr::bind_rows(
  tibble::tibble(
    ecuacion = "Inflación",
    parametro = names(stats::coef(init$eq_pi)),
    estimacion = as.numeric(stats::coef(init$eq_pi))
  ),
  tibble::tibble(
    ecuacion = "Expectativa 24m",
    parametro = names(stats::coef(init$eq_e)),
    estimacion = as.numeric(stats::coef(init$eq_e))
  )
)

readr::write_csv(
  tabla_OLS,
  file.path(
    carpeta_salida,
    "02_coeficientes_OLS_inicializacion.csv"
  )
)


# 5. CONSTRUCCIÓN DEL MODELO EN ESPACIO DE ESTADOS ----------------------------

# Estados:
# 1 = a0
# 2 = a
# 3 = b
# 4 = c0_t
# 5 = c_t
# 6 = d
#
# M0: Q = 0 para todos los estados.
# M3: Q_c0 > 0 y Q_c > 0; los demás Q = 0.

crear_modelo <- function(
    datos,
    Y,
    estado_inicial,
    modelo = c("M0", "M3")
) {

  modelo <- match.arg(modelo) #diseñamos dos modelos M0 y M3.

  n <- nrow(datos) #T
  m <- 6L #6 estados, L = entero.

  Z_t <- array(
    0,
    dim = c(2, m, n) #<-dimensión 2 (variables obseradas), * m estados * n variables, 
  ) 
  
  #la Ecuación de estado espacio es yt = Z_t*alpha_t + error, donde yt es el vector de variables 
  #observables, y alpha_t es el vector de estados
  
  # Ecuación de inflación:
  # pi_t = a0 + a*pi_(t-1) + b*e_(t-1) + eps_pi,t
  Z_t[1, 1, ] <- 1
  Z_t[1, 2, ] <- datos$inflacion_lag1 #asignando en cada elemento de la matriz, el valor de las variables observadas
  Z_t[1, 3, ] <- datos$expectativa_lag1

  # Ecuación de expectativas:
  # e_t = c0_t + c_t*pi_(t-1) + d*e_(t-1) + eps_e,t
  Z_t[2, 4, ] <- 1
  Z_t[2, 5, ] <- datos$inflacion_lag1
  Z_t[2, 6, ] <- datos$expectativa_lag1

  # básicamente al realizar el producto entre el vector de estados por esta matriz z recuperamos
  #nuestras dos ecuaciones de inflación y expectativas.
    Q_mat <- matrix( # se construye q(var(n_t)), Q = cuanto pueden moverse los parametros.
    0,
    nrow = m,
    ncol = m
  )

  if (modelo == "M3") {
    Q_mat[4, 4] <- NA_real_
    Q_mat[5, 5] <- NA_real_
  }

  KFAS::SSModel(
    Y ~ -1 +
      SSMcustom(
        Z = Z_t,
        T = diag(m),
        R = diag(m),
        Q = Q_mat,
        a1 = estado_inicial,
        P1 = matrix(0, m, m),
        P1inf = diag(m), #<- los seis estados se inicializan como difusos, 
        n = n #es decir, con incertidumbre inicial esencialmente no informativa.
      ),
    H = diag(rep(NA_real_, 2)) #La documentación de KFAS define H como la matriz de covarianzas de los errores de observación y
  ) # Q como la matriz de covarianzas de las innovaciones de estado
}



crear_updatefn <- function(modelo = c("M0", "M3")) {

  modelo <- match.arg(modelo)

  function(par, objeto) {

    # Varianzas de observación.
    objeto$H[1, 1, 1] <- exp(par[1]) #par es el vector de parámetros que el optimizador está probando en cada iteración
    objeto$H[2, 2, 1] <- exp(par[2])
    objeto$H[1, 2, 1] <- 0
    objeto$H[2, 1, 1] <- 0

    # M3: varianzas de transición de c0_t y c_t.
    if (modelo == "M3") {
      objeto$Q[4, 4, 1] <- exp(par[3])
      objeto$Q[5, 5, 1] <- exp(par[4])
    } #<- M3 tiene cuatro hiperparámetros de varianza que el optimizador debe encontrar
#Cuando hablamos de hiperparámetros, nos referimos a parámetros que no aparecen directamente como coeficientes económicos de las ecuaciones, 
 #sino que gobiernan el comportamiento estadístico del modelo.
    objeto
  }
}

#hiperaparámetros: H dice qué tan ruidosas son las ecuaciones;
#Q dice qué tan móviles son los parámetros.
crear_inicios <- function(
    H_inicial,
    modelo = c("M0", "M3")
) {

  modelo <- match.arg(modelo)

  if (modelo == "M0") {

    return(
      list(
        log(H_inicial * 0.5),
        log(H_inicial),
        log(H_inicial * 2)
      )
    )
  }

  # Distintos órdenes de magnitud para reducir el riesgo de máximos locales.
  Q_inicios <- list(
    c(1e-8, 1e-8),
    c(1e-6, 1e-6),
    c(1e-4, 1e-4),
    c(1e-3, 1e-6),
    c(1e-6, 1e-3),
    c(1e-3, 1e-4),
    c(1e-4, 1e-3)
  )

  inicios <- lapply(
    Q_inicios,
    function(q) {
      log(
        c(
          H_inicial,
          q
        )
      )
    }
  )

  inicios <- c(
    inicios,
    list(
      log(c(H_inicial * 0.5, 1e-5, 1e-5)),
      log(c(H_inicial * 2.0, 1e-5, 1e-5))
    )
  )

  inicios
}


estimar_modelo <- function(
    codigo,
    datos,
    Y,
    estado_inicial,
    H_inicial,
    maxit
) {

  modelo_base <- crear_modelo(
    datos = datos,
    Y = Y,
    estado_inicial = estado_inicial,
    modelo = codigo
  )

  updatefn <- crear_updatefn(codigo)
  inicios <- crear_inicios(H_inicial, codigo)

  numero_parametros <- if (codigo == "M0") 2L else 4L

  resultados <- lapply(
    seq_along(inicios),
    function(i) {

      fit <- tryCatch(
        KFAS::fitSSM( #por maximoverosimilitud en modelos estado-espacio
          model = modelo_base,
          inits = inicios[[i]],
          updatefn = updatefn,
          method = "L-BFGS-B", #para encontrar el máximo
          lower = rep( #permite colocar limites abajo y arriba
            log(limite_varianza),
            numero_parametros
          ),
          upper = rep(
            log(100),
            numero_parametros
          ),
          control = list(
            maxit = maxit,
            factr = 1e7
          )
        ),
        error = function(e) NULL
      )

      if (is.null(fit)) {

        return(
          list(
            intento = i,
            convergencia = NA_integer_,
            logLik = NA_real_,
            fit = NULL
          )
        )
      }

      ll <- tryCatch(
        as.numeric(stats::logLik(fit$model)),
        error = function(e) NA_real_
      )

      list(
        intento = i,
        convergencia = fit$optim.out$convergence,
        logLik = ll,
        fit = fit
      )
    }
  )

  tabla_intentos <- dplyr::bind_rows(
    lapply(
      resultados,
      function(x) {
        tibble::tibble(
          modelo = codigo,
          intento = x$intento,
          convergencia = x$convergencia,
          log_verosimilitud = x$logLik
        )
      }
    )
  )

  validos <- resultados[
    vapply(
      resultados,
      function(x) {
        !is.null(x$fit) && is.finite(x$logLik)
      },
      logical(1)
    )
  ]

  if (length(validos) == 0) {
    stop(
      "No fue posible estimar ", codigo,
      ". Revise los valores iniciales y los datos."
    )
  }

  convergentes <- validos[
    vapply(
      validos,
      function(x) {
        isTRUE(x$convergencia == 0)
      },
      logical(1)
    )
  ]

  candidatos <- if (length(convergentes) > 0) {
    convergentes
  } else {
    warning(
      "No hubo convergencia formal en ", codigo,
      "; se conserva el intento con mayor log-verosimilitud."
    )
    validos
  }

  mejor <- candidatos[[
    which.max(
      vapply(
        candidatos,
        function(x) x$logLik,
        numeric(1)
      )
    )
  ]]

  list(
    codigo = codigo,
    fit = mejor$fit,
    modelo = mejor$fit$model,
    intento = mejor$intento,
    convergencia = mejor$convergencia,
    logLik = mejor$logLik,
    tabla_intentos = tabla_intentos
  )
}


# 6. ESTIMACIÓN: SOLO M0 Y M3 --------------------------------------------------

cat("\nEstimando M0: benchmark de parámetros constantes...\n")

ajuste_M0 <- estimar_modelo(
  codigo = "M0",
  datos = datos_modelo,
  Y = Y_modelo,
  estado_inicial = estado_inicial,
  H_inicial = H_inicial,
  maxit = maxit_M0
)

cat("\nEstimando M3: TVP-VAR restringido principal...\n")

ajuste_M3 <- estimar_modelo(
  codigo = "M3",
  datos = datos_modelo,
  Y = Y_modelo,
  estado_inicial = estado_inicial,
  H_inicial = H_inicial,
  maxit = maxit_M3
)

tabla_intentos <- dplyr::bind_rows(
  ajuste_M0$tabla_intentos,
  ajuste_M3$tabla_intentos
)

readr::write_csv(
  tabla_intentos,
  file.path(
    carpeta_salida,
    "03_intentos_estimacion.csv"
  )
)


# 7. FILTRO Y SUAVIZADOR -------------------------------------------------------

kfs_M0 <- KFAS::KFS( #¿qué creía el modelo sobre ct usando únicamente la información disponible hasta t?
  ajuste_M0$modelo,
  filtering = c("state", "mean"),
  smoothing = c("state", "mean", "disturbance") 
)

#el suavizado ¿cuál es la mejor estimación de ct, utilizando toda la muestra, incluso información posterior a t?
kfs_M3 <- KFAS::KFS(
  ajuste_M3$modelo,
  filtering = c("state", "mean"),
  smoothing = c("state", "mean", "disturbance")
)


# 8. HIPERPARÁMETROS Y COMPARACIÓN M0-M3 --------------------------------------

extraer_hiperparametros <- function(
    codigo,
    ajuste
) {

  H <- diag(ajuste$modelo$H[, , 1])
  Q <- diag(ajuste$modelo$Q[, , 1])

  tibble::tibble(
    modelo = codigo,
    H_inflacion = H[1],
    H_expectativa = H[2],
    Q_a0 = Q[1],
    Q_a = Q[2],
    Q_b = Q[3],
    Q_c0 = Q[4],
    Q_c = Q[5],
    Q_d = Q[6],
    convergencia = ajuste$convergencia,
    intento_seleccionado = ajuste$intento,
    log_verosimilitud = ajuste$logLik
  )
}

tabla_hiperparametros <- dplyr::bind_rows(
  extraer_hiperparametros("M0", ajuste_M0),
  extraer_hiperparametros("M3", ajuste_M3)
)

readr::write_csv(
  tabla_hiperparametros,
  file.path(
    carpeta_salida,
    "04_hiperparametros_M0_M3.csv"
  )
)

# Conteo para criterios de información:
# 6 coeficientes del VAR + 2 H + Q activas.
k_M0 <- 8
k_M3 <- 10

ll_M0 <- ajuste_M0$logLik
ll_M3 <- ajuste_M3$logLik

AIC_M0 <- -2 * ll_M0 + 2 * k_M0
AIC_M3 <- -2 * ll_M3 + 2 * k_M3

BIC_M0 <- -2 * ll_M0 + log(n_obs) * k_M0
BIC_M3 <- -2 * ll_M3 + log(n_obs) * k_M3

LR_M3_M0 <- 2 * (ll_M3 - ll_M0)

if (LR_M3_M0 < 0) {
  warning(
    "M3 presenta una log-verosimilitud menor que M0. ",
    "Revise convergencia e inicialización."
  )
}

LR_referencial <- max(0, LR_M3_M0)

# Este p-valor es SOLO referencial:
# Q_c0 = Q_c = 0 está en la frontera del espacio paramétrico bajo H0.
p_LR_chi2_referencial <- stats::pchisq(
  LR_referencial,
  df = 2,
  lower.tail = FALSE
)

tabla_comparacion <- tibble::tibble(
  modelo = c("VAR", "TVP-VAR"),
  descripcion = c(
    "VAR(1) constante",
    "TVP-VAR restringido: c0_t y c_t variantes"
  ),
  log_verosimilitud = c(ll_M0, ll_M3),
  parametros_contados = c(k_M0, k_M3),
  AIC = c(AIC_M0, AIC_M3),
  BIC = c(BIC_M0, BIC_M3)
)

readr::write_csv(
  tabla_comparacion,
  file.path(
    carpeta_salida,
    "05_comparacion_M0_M3.csv"
  )
)


# 9. TRAYECTORIA DE ESTADOS M3 ------------------------------------------------

nombres_estados <- c(
  "a0_t",
  "a_t",
  "b_t",
  "c0_t",
  "c_t",
  "d_t"
)

estados_M3 <- as.data.frame(kfs_M3$alphahat) #estados suavizados
names(estados_M3) <- nombres_estados

trayectoria_M3 <- dplyr::bind_cols(
  datos_modelo |>
    dplyr::select(
      fecha,
      inflacion,
      expectativa,
      meta_inflacion
    ),
  estados_M3
)

# Errores estándar de c0_t y c_t a partir de la matriz suavizada V.
trayectoria_M3$c0_se <- sqrt(
  pmax(
    kfs_M3$V[4, 4, ],
    0
  )
)

trayectoria_M3$c_se <- sqrt(
  pmax(
    kfs_M3$V[5, 5, ],
    0
  )
)

trayectoria_M3 <- trayectoria_M3 |>
  dplyr::mutate(
    c0_li_95 = c0_t - z_critico * c0_se,
    c0_ls_95 = c0_t + z_critico * c0_se,
    c_li_95 = c_t - z_critico * c_se,
    c_ls_95 = c_t + z_critico * c_se
  )


# 10. PROXY LAMBDA Y ANCLA IMPLÍCITA ------------------------------------------

calcular_proxy_t <- function(t) { #p/cada período t, toma los parámetros estimados y calcula la proxy de credibilidad y el ancla implícita.

  c0 <- trayectoria_M3$c0_t[t]
  c_coef <- trayectoria_M3$c_t[t]
  d_coef <- trayectoria_M3$d_t[t]

  den_lambda <- 1 - d_coef
  den_ancla <- 1 - d_coef - c_coef

  if (
    !is.finite(den_lambda) ||
      !is.finite(den_ancla) ||
      abs(den_lambda) < tolerancia_denominador ||
      abs(den_ancla) < tolerancia_denominador
  ) {

    return(
      tibble::tibble(
        lambda_t = NA_real_,
        lambda_se = NA_real_,
        pi_estrella_t = NA_real_,
        pi_estrella_se = NA_real_,
        denominador_lambda = den_lambda,
        denominador_ancla = den_ancla
      )
    )
  }

  # Covarianza suavizada de (c0_t, c_t, d).
  V_sub <- kfs_M3$V[ #matriz de covarianza de los estados
    c(4, 5, 6),
    c(4, 5, 6),
    t,
    drop = FALSE
  ]

  V_sub <- matrix(
    V_sub,
    nrow = 3,
    ncol = 3
  )

  lambda <- 1 - c_coef / den_lambda #indicador de anclaje

  pi_estrella <- c0 / den_ancla #indicador de ancla ímplicita

  # Método delta:
  # lambda_t y pi_estrella_t son funciones no lineales de los parámetros
  # estimados. El método delta aproxima su varianza utilizando la matriz de
  # covarianzas de los estados y el gradiente de cada función respecto de
  # (c0_t, c_t, d).
  
  grad_lambda <- matrix(
    c(
      0,
      -1 / den_lambda,
      -c_coef / den_lambda^2
    ),
    ncol = 1
  )
  # Gradientes:
  # Se calculan las derivadas parciales de lambda_t y pi_estrella_t respecto de
  # (c0_t, c_t, d). Estas derivadas indican cómo cambia cada medida ante pequeñas
  # variaciones en los parámetros y permiten propagar su incertidumbre mediante
  # el método delta.
  grad_pi <- matrix(
    c(
      1 / den_ancla,
      c0 / den_ancla^2,
      c0 / den_ancla^2
    ),
    ncol = 1
  )
  
  #estimación de varianzas

  var_lambda <- as.numeric(
    t(grad_lambda) %*%
      V_sub %*%
      grad_lambda
  )

  var_pi <- as.numeric(
    t(grad_pi) %*%
      V_sub %*%
      grad_pi
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
    denominador_lambda = den_lambda,
    denominador_ancla = den_ancla
  )
}

proxy_M3 <- dplyr::bind_rows(
  lapply(
    seq_len(nrow(trayectoria_M3)),
    calcular_proxy_t
  )
)

trayectoria_M3 <- dplyr::bind_cols(
  trayectoria_M3,
  proxy_M3
) |>
  dplyr::mutate(
    lambda_li_95 = lambda_t - z_critico * lambda_se,
    lambda_ls_95 = lambda_t + z_critico * lambda_se,
    pi_estrella_li_95 = pi_estrella_t - z_critico * pi_estrella_se,
    pi_estrella_ls_95 = pi_estrella_t + z_critico * pi_estrella_se,

    lambda_fuera_0_1 =
      lambda_t < 0 |
      lambda_t > 1,

    lambda_no_distinta_de_1 =
      lambda_li_95 <= 1 &
      lambda_ls_95 >= 1,

    brecha_ancla_meta =
      pi_estrella_t -
      meta_inflacion,

    banda_inferior =
      meta_inflacion - 1,

    banda_superior =
      meta_inflacion + 1,

    ancla_dentro_banda_meta =
      pi_estrella_t >= banda_inferior &
      pi_estrella_t <= banda_superior
  )


# 11. ESTABILIDAD LOCAL DEL TVP-VAR M3 ----------------------------------------

# Para un VAR(1), en cada t:
#
#       [ a     b ]
# A_t = [         ]
#       [ c_t   d ]
#
# El sistema es localmente estable si todas las raíces tienen módulo < 1.

calcular_estabilidad_t <- function(t) {

  A_t <- matrix(
    c(
      trayectoria_M3$a_t[t],
      trayectoria_M3$b_t[t],
      trayectoria_M3$c_t[t],
      trayectoria_M3$d_t[t]
    ),
    nrow = 2,
    byrow = TRUE
  )

  raices <- eigen(
    A_t,
    only.values = TRUE
  )$values

  tibble::tibble(
    fecha = trayectoria_M3$fecha[t],
    raiz_1_real = Re(raices[1]),
    raiz_1_imaginaria = Im(raices[1]),
    raiz_2_real = Re(raices[2]),
    raiz_2_imaginaria = Im(raices[2]),
    radio_espectral = max(Mod(raices)),
    estable = max(Mod(raices)) < 1
  )
}

tabla_estabilidad <- dplyr::bind_rows(
  lapply(
    seq_len(nrow(trayectoria_M3)),
    calcular_estabilidad_t
  )
)

trayectoria_M3 <- trayectoria_M3 |>
  dplyr::left_join(
    tabla_estabilidad |>
      dplyr::select(
        fecha,
        radio_espectral,
        estable
      ),
    by = "fecha"
  )

readr::write_csv(
  trayectoria_M3,
  file.path(
    carpeta_salida,
    "06_trayectoria_M3.csv"
  )
)

readr::write_csv(
  tabla_estabilidad,
  file.path(
    carpeta_salida,
    "10_estabilidad_local_M3.csv"
  )
)


# 12. RESUMEN DEL MODELO M3 ---------------------------------------------------

indice_min_lambda <- which.min(trayectoria_M3$lambda_t)
indice_max_lambda <- which.max(trayectoria_M3$lambda_t)

resumen_M3 <- tibble::tibble(
  fecha_inicial = dplyr::first(trayectoria_M3$fecha),
  fecha_final = dplyr::last(trayectoria_M3$fecha),
  observaciones = nrow(trayectoria_M3),

  lambda_promedio = mean(
    trayectoria_M3$lambda_t,
    na.rm = TRUE
  ),

  lambda_inicial = dplyr::first(
    trayectoria_M3$lambda_t
  ),

  lambda_final = dplyr::last(
    trayectoria_M3$lambda_t
  ),

  lambda_minimo =
    trayectoria_M3$lambda_t[
      indice_min_lambda
    ],

  fecha_lambda_minimo =
    trayectoria_M3$fecha[
      indice_min_lambda
    ],

  lambda_maximo =
    trayectoria_M3$lambda_t[
      indice_max_lambda
    ],

  fecha_lambda_maximo =
    trayectoria_M3$fecha[
      indice_max_lambda
    ],

  porcentaje_lambda_fuera_0_1 =
    100 * mean(
      trayectoria_M3$lambda_fuera_0_1,
      na.rm = TRUE
    ),

  pi_estrella_promedio =
    mean(
      trayectoria_M3$pi_estrella_t,
      na.rm = TRUE
    ),

  pi_estrella_inicial =
    dplyr::first(
      trayectoria_M3$pi_estrella_t
    ),

  pi_estrella_final =
    dplyr::last(
      trayectoria_M3$pi_estrella_t
    ),

  MAE_pi_estrella_meta =
    mean(
      abs(
        trayectoria_M3$brecha_ancla_meta
      ),
      na.rm = TRUE
    ),

  RMSE_pi_estrella_meta =
    sqrt(
      mean(
        trayectoria_M3$brecha_ancla_meta^2,
        na.rm = TRUE
      )
    ),

  porcentaje_ancla_dentro_banda_meta =
    100 * mean(
      trayectoria_M3$ancla_dentro_banda_meta,
      na.rm = TRUE
    ),

  porcentaje_meses_estables =
    100 * mean(
      trayectoria_M3$estable,
      na.rm = TRUE
    ),

  radio_espectral_maximo =
    max(
      trayectoria_M3$radio_espectral,
      na.rm = TRUE
    ),

  c0_promedio =
    mean(
      trayectoria_M3$c0_t,
      na.rm = TRUE
    ),

  c_promedio =
    mean(
      trayectoria_M3$c_t,
      na.rm = TRUE
    ),

  d_estimado =
    mean(
      trayectoria_M3$d_t,
      na.rm = TRUE
    )
)

readr::write_csv(
  resumen_M3,
  file.path(
    carpeta_salida,
    "07_resumen_M3.csv"
  )
)

# Subperíodos descriptivos.
resumen_subperiodos <- trayectoria_M3 |>
  dplyr::mutate(
    subperiodo = dplyr::case_when(
      fecha <= as.Date("2014-12-01") ~ "2010-2014",
      fecha <= as.Date("2019-12-01") ~ "2015-2019",
      fecha <= as.Date("2022-12-01") ~ "2020-2022",
      TRUE ~ "2023-2026"
    )
  ) |>
  dplyr::group_by(subperiodo) |>
  dplyr::summarise(
    fecha_inicial = min(fecha),
    fecha_final = max(fecha),
    observaciones = dplyr::n(),
    lambda_promedio = mean(lambda_t, na.rm = TRUE),
    pi_estrella_promedio = mean(pi_estrella_t, na.rm = TRUE),
    meta_promedio = mean(meta_inflacion, na.rm = TRUE),
    brecha_ancla_meta_promedio = mean(
      brecha_ancla_meta,
      na.rm = TRUE
    ),
    c0_promedio = mean(c0_t, na.rm = TRUE),
    c_promedio = mean(c_t, na.rm = TRUE),
    d_promedio = mean(d_t, na.rm = TRUE),
    radio_espectral_maximo = max(
      radio_espectral,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

readr::write_csv(
  resumen_subperiodos,
  file.path(
    carpeta_salida,
    "08_resumen_subperiodos_M3.csv"
  )
)


# 13. DIAGNÓSTICOS RESIDUALES M3 ----------------------------------------------

extraer_residuos_estandarizados <- function(
    kfs,
    n_esperado
) {

  r <- tryCatch(
    as.matrix(
      stats::rstandard(
        kfs,
        type = "recursive"
      )
    ),
    error = function(e) NULL
  )

  if (
    is.null(r) ||
      nrow(r) != n_esperado ||
      ncol(r) != 2
  ) {

    stop(
      "No fue posible extraer correctamente los residuos ",
      "estandarizados recursivos del M3."
    )
  }

  colnames(r) <- c(
    "inflacion",
    "expectativa"
  )

  r
}

residuos_M3 <- extraer_residuos_estandarizados(
  kfs_M3,
  n_obs
)

tabla_residuos <- tibble::tibble(
  fecha = datos_modelo$fecha,
  residuo_inflacion = residuos_M3[, "inflacion"],
  residuo_expectativa = residuos_M3[, "expectativa"]
)

readr::write_csv(
  tabla_residuos,
  file.path(
    carpeta_salida,
    "11_residuos_M3.csv"
  )
)

crear_diagnosticos_ecuacion <- function(
    residuos,
    nombre_ecuacion
) {

  r <- residuos[is.finite(residuos)]

  filas_LB <- lapply(
    rezagos_ljung_box,
    function(L) {

      if (length(r) <= L + 5) {
        return(NULL)
      }

      prueba <- stats::Box.test(
        r,
        lag = L,
        type = "Ljung-Box",
        fitdf = 0
      )

      tibble::tibble(
        ecuacion = nombre_ecuacion,
        diagnostico = "Ljung-Box",
        rezago = L,
        estadistico = unname(prueba$statistic),
        p_valor = prueba$p.value
      )
    }
  )

  JB <- jarque_bera_manual(r)

  fila_JB <- tibble::tibble(
    ecuacion = nombre_ecuacion,
    diagnostico = "Jarque-Bera",
    rezago = NA_real_,
    estadistico = JB$estadistico,
    p_valor = JB$p_valor
  )

  dplyr::bind_rows(
    filas_LB,
    fila_JB
  )
}

tabla_diagnosticos <- dplyr::bind_rows(
  crear_diagnosticos_ecuacion(
    residuos_M3[, "inflacion"],
    "Inflación"
  ),
  crear_diagnosticos_ecuacion(
    residuos_M3[, "expectativa"],
    "Expectativa 24m"
  )
)

# Correlación contemporánea de los residuos.
rho_residuos <- correlacion_segura(
  residuos_M3[, "inflacion"],
  residuos_M3[, "expectativa"]
)

n_rho <- sum(
  is.finite(residuos_M3[, "inflacion"]) &
    is.finite(residuos_M3[, "expectativa"])
)

fisher_z <- if (
  is.finite(rho_residuos) &&
    abs(rho_residuos) < 1 &&
    n_rho > 3
) {
  atanh(rho_residuos) *
    sqrt(n_rho - 3)
} else {
  NA_real_
}

p_fisher <- if (is.finite(fisher_z)) {
  2 * stats::pnorm(
    abs(fisher_z),
    lower.tail = FALSE
  )
} else {
  NA_real_
}

tabla_diagnosticos <- dplyr::bind_rows(
  tabla_diagnosticos,
  tibble::tibble(
    ecuacion = "Sistema",
    diagnostico = "Correlación contemporánea Fisher",
    rezago = NA_real_,
    estadistico = fisher_z,
    p_valor = p_fisher
  )
) |>
  dplyr::mutate(
    conclusion_5pct = dplyr::case_when(
      !is.finite(p_valor) ~ "No disponible",
      p_valor < 0.05 ~ "Se rechaza H0",
      TRUE ~ "No se rechaza H0"
    )
  )

readr::write_csv(
  tabla_diagnosticos,
  file.path(
    carpeta_salida,
    "09_diagnosticos_residuales_M3.csv"
  )
)


# ==============================================================================
# 14. GRÁFICAS CENTRALES - LÍNEA GRÁFICA DEL WORKING PAPER
# ==============================================================================
#
# Este bloque reemplaza por completo la sección 14 del script base.
# Conserva los mismos nombres de archivo para no alterar el flujo del proyecto.
#
# Principios gráficos:
# - sin títulos ni subtítulos internos: se administran desde Quarto;
# - línea principal azul grisáceo;
# - bandas de confianza azul claro;
# - fondo blanco y cuadrícula tenue;
# - eje temporal con marcas cada 2 años y menores cada 1 año;
# - tamaño 8.5 x 4.3, 300 dpi, consistente con rolling regressions.
# ==============================================================================


# 14.0 Estilo común ------------------------------------------------------------

col_serie <- "#2F4A67"
col_banda <- "#AFC4D6"
col_ref   <- "#4D4D4D"

# Colores secundarios para figuras con varias series.
col_expectativa <- "#6F7D8C"
col_meta        <- "#3A3A3A"
col_inflacion   <- "#A0A0A0"

tema_wp <- ggplot2::theme_minimal(base_size = 10.5) +
  ggplot2::theme(
    panel.grid.minor.x = ggplot2::element_line(
      colour = "grey94",
      linewidth = 0.25
    ),
    panel.grid.minor.y = ggplot2::element_blank(),
    
    panel.grid.major.x = ggplot2::element_line(
      colour = "grey90",
      linewidth = 0.30
    ),
    panel.grid.major.y = ggplot2::element_line(
      colour = "grey87",
      linewidth = 0.35
    ),
    
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(
      size = 10.5,
      colour = "black",
      margin = ggplot2::margin(r = 8)
    ),
    
    axis.text = ggplot2::element_text(
      size = 9,
      colour = "grey25"
    ),
    
    axis.ticks = ggplot2::element_blank(),
    
    # Los títulos y notas se colocan desde Quarto.
    plot.title = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank(),
    
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(
      size = 9,
      colour = "grey25"
    ),
    
    plot.margin = ggplot2::margin(
      t = 4,
      r = 8,
      b = 4,
      l = 4
    )
  )


escala_x_wp <- ggplot2::scale_x_date(
  date_breaks = "2 years",
  date_labels = "%Y",
  date_minor_breaks = "1 year",
  expand = ggplot2::expansion(
    mult = c(0.01, 0.015)
  )
)


# 14.1 Intercepto variante c0_t -----------------------------------------------

p_c0 <- ggplot2::ggplot(
  trayectoria_M3,
  ggplot2::aes(
    x = fecha,
    y = c0_t
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = c0_li_95,
      ymax = c0_ls_95
    ),
    fill = col_banda,
    alpha = 0.28,
    linewidth = 0
  ) +
  ggplot2::geom_line(
    colour = col_serie,
    linewidth = 0.80,
    lineend = "round"
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey30",
    linewidth = 0.45
  ) +
  escala_x_wp +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.04, 0.06)
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = expression(hat(c)[0*t])
  ) +
  tema_wp

ggplot2::ggsave(
  file.path(
    carpeta_salida,
    "figura_01_c0_t.png"
  ),
  plot = p_c0,
  width = 8.5,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

print(p_c0)


# 14.2 Sensibilidad cambiante c_t ---------------------------------------------

p_c <- ggplot2::ggplot(
  trayectoria_M3,
  ggplot2::aes(
    x = fecha,
    y = c_t
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = c_li_95,
      ymax = c_ls_95
    ),
    fill = col_banda,
    alpha = 0.28,
    linewidth = 0
  ) +
  ggplot2::geom_line(
    colour = col_serie,
    linewidth = 0.80,
    lineend = "round"
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey30",
    linewidth = 0.45
  ) +
  escala_x_wp +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.04, 0.06)
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = expression(hat(c)[t])
  ) +
  tema_wp

ggplot2::ggsave(
  file.path(
    carpeta_salida,
    "figura_02_c_t.png"
  ),
  plot = p_c,
  width = 8.5,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

print(p_c)

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}

library(patchwork)

p_c0_panel <- p_c0 +
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 4,
      r = 6,
      b = 4,
      l = 4
    )
  )

p_c_panel <- p_c +
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 4,
      r = 4,
      b = 4,
      l = 6
    )
  )

panel_parametros_tvp <- (
  p_c0_panel +
    p_c_panel
) +
  patchwork::plot_layout(
    ncol = 2,
    widths = c(1, 1)
  ) +
  patchwork::plot_annotation(
    tag_levels = "a"
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      size = 10,
      face = "bold",
      colour = "grey25"
    ),
    plot.tag.position = c(0.02, 0.98)
  )

ggplot2::ggsave(
  file.path(
    carpeta_salida,
    "figura_02_panel_parametros_tvp.png"
  ),
  plot = panel_parametros_tvp,
  width = 8.5,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

print(panel_parametros_tvp)

# 14.3 Proxy dinámica de credibilidad lambda_t --------------------------------
#
# Se omite la línea horizontal en cero porque no aporta información visual
# en la región empírica estimada. La referencia relevante es lambda = 1.

p_lambda <- ggplot2::ggplot(
  trayectoria_M3,
  ggplot2::aes(
    x = fecha,
    y = lambda_t
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = lambda_li_95,
      ymax = lambda_ls_95
    ),
    fill = col_banda,
    alpha = 0.28,
    linewidth = 0
  ) +
  ggplot2::geom_line(
    colour = col_serie,
    linewidth = 0.80,
    lineend = "round"
  ) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype = "dotted",
    colour = "grey45",
    linewidth = 0.45
  ) +
  escala_x_wp +
  ggplot2::scale_y_continuous(
    breaks = seq(
      0.6,
      1.0,
      by = 0.1
    ),
    expand = ggplot2::expansion(
      mult = c(0.03, 0.05)
    )
  ) +
  ggplot2::coord_cartesian(
    ylim = c(
      min(
        trayectoria_M3$lambda_li_95,
        na.rm = TRUE
      ) - 0.015,
      1.015
    ),
    clip = "off"
  ) +
  ggplot2::labs(
    x = NULL,
    y = expression(hat(lambda)[t])
  ) +
  tema_wp

ggplot2::ggsave(
  file.path(
    carpeta_salida,
    "figura_03_lambda_t.png"
  ),
  plot = p_lambda,
  width = 8.5,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

print(p_lambda)


# 14.4 Ancla implícita, expectativas y meta -----------------------------------

p_ancla <- ggplot2::ggplot(
  trayectoria_M3,
  ggplot2::aes(x = fecha)
) +
  
  # Rango meta oficial: referencia institucional de fondo
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = banda_inferior,
      ymax = banda_superior
    ),
    fill = col_banda,
    alpha = 0.10,
    linewidth = 0
  ) +
  
  # IC 95% del ancla implícita:
  # visible, pero deliberadamente tenue
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = pi_estrella_li_95,
      ymax = pi_estrella_ls_95
    ),
    fill = col_serie,
    alpha = 0.07,
    linewidth = 0
  ) +
  
  # Expectativa a 24 meses:
  # serie secundaria frente al ancla implícita
  ggplot2::geom_line(
    ggplot2::aes(
      y = expectativa,
      colour = "Expectativa a 24 meses",
      linetype = "Expectativa a 24 meses"
    ),
    linewidth = 0.55,
    alpha = 0.65,
    lineend = "round"
  ) +
  
  # Meta central
  ggplot2::geom_line(
    ggplot2::aes(
      y = meta_inflacion,
      colour = "Meta central",
      linetype = "Meta central"
    ),
    linewidth = 0.55,
    alpha = 0.90
  ) +
  
  # Ancla implícita: objeto principal de la figura
  ggplot2::geom_line(
    ggplot2::aes(
      y = pi_estrella_t,
      colour = "Ancla inflacionaria implícita",
      linetype = "Ancla inflacionaria implícita"
    ),
    linewidth = 0.95,
    lineend = "round"
  ) +
  
  ggplot2::scale_colour_manual(
    values = c(
      "Ancla inflacionaria implícita" = col_serie,
      "Expectativa a 24 meses" = col_expectativa,
      "Meta central" = col_meta
    ),
    breaks = c(
      "Ancla inflacionaria implícita",
      "Expectativa a 24 meses",
      "Meta central"
    )
  ) +
  
  ggplot2::scale_linetype_manual(
    values = c(
      "Ancla inflacionaria implícita" = "solid",
      "Expectativa a 24 meses" = "longdash",
      "Meta central" = "dotted"
    ),
    breaks = c(
      "Ancla inflacionaria implícita",
      "Expectativa a 24 meses",
      "Meta central"
    )
  ) +
  
  escala_x_wp +
  
  ggplot2::scale_y_continuous(
    limits = c(2, 7),
    breaks = seq(2, 7, by = 1),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "Porcentaje",
    colour = NULL,
    linetype = NULL
  ) +
  
  tema_wp +
  
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.spacing.x = grid::unit(0.18, "cm"),
    legend.key.width = grid::unit(1.2, "cm"),
    legend.text = ggplot2::element_text(size = 8.5)
  )


ggplot2::ggsave(
  file.path(
    carpeta_salida,
    "figura_04_ancla_implicita.png"
  ),
  plot = p_ancla,
  width = 8.5,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

print(p_ancla)

# 14.5 Dinámica conjunta: inflación, expectativas, ancla implícita y credibilidad ----

# Figura sintética inspirada en la lógica de Demertzis:
# eje izquierdo  -> inflación observada, expectativa a 24 meses, ancla implícita y meta
# eje derecho    -> proxy de credibilidad lambda_t
#
# La intención de esta figura no es reemplazar las dos anteriores,
# sino condensar en una sola lectura la relación entre:
#   (i) la inflación efectiva,
#   (ii) el nivel hacia el cual parecen converger las expectativas, y
#   (iii) el grado de desacoplamiento de las expectativas respecto de la inflación.

# -------------------------------------------------------------------
# 1. Parámetros de escala para el doble eje
# -------------------------------------------------------------------

lambda_min <- 0.70
lambda_max <- 1.00

y_min <- 0
y_max <- 8

# Transformación lineal para dibujar lambda_t sobre el eje izquierdo
trayectoria_M3 <- trayectoria_M3 |>
  dplyr::mutate(
    lambda_plot = y_min +
      (lambda_t - lambda_min) / (lambda_max - lambda_min) * (y_max - y_min)
  )

# -------------------------------------------------------------------
# 2. Colores auxiliares
# -------------------------------------------------------------------

col_inflacion <- "black"
col_lambda <- "#2C6BE0"   # puedes ajustarlo si ya tienes un azul institucional

# -------------------------------------------------------------------
# 3. Gráfica
# -------------------------------------------------------------------

p_sintesis_ancla <- ggplot2::ggplot(
  trayectoria_M3,
  ggplot2::aes(x = fecha)
) +
  
  # Rango meta oficial
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = banda_inferior,
      ymax = banda_superior
    ),
    fill = col_banda,
    alpha = 0.10,
    linewidth = 0
  ) +
  
  # Ancla implícita: banda de confianza tenue
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = pi_estrella_li_95,
      ymax = pi_estrella_ls_95
    ),
    fill = col_serie,
    alpha = 0.06,
    linewidth = 0
  ) +
  
  # Inflación observada
  ggplot2::geom_line(
    ggplot2::aes(
      y = inflacion,
      colour = "Inflación observada",
      linetype = "Inflación observada"
    ),
    linewidth = 0.70,
    alpha = 0.85,
    lineend = "round"
  ) +
  
  # Expectativa a 24 meses
  ggplot2::geom_line(
    ggplot2::aes(
      y = expectativa,
      colour = "Expectativa a 24 meses",
      linetype = "Expectativa a 24 meses"
    ),
    linewidth = 0.60,
    alpha = 0.65,
    lineend = "round"
  ) +
  
  # Meta central
  ggplot2::geom_line(
    ggplot2::aes(
      y = meta_inflacion,
      colour = "Meta central",
      linetype = "Meta central"
    ),
    linewidth = 0.55,
    alpha = 0.95
  ) +
  
  # Ancla implícita
  ggplot2::geom_line(
    ggplot2::aes(
      y = pi_estrella_t,
      colour = "Ancla inflacionaria implícita",
      linetype = "Ancla inflacionaria implícita"
    ),
    linewidth = 0.95,
    lineend = "round"
  ) +
  
  # Lambda_t reescalado para el eje derecho
  ggplot2::geom_line(
    ggplot2::aes(
      y = lambda_plot,
      colour = "Proxy de credibilidad",
      linetype = "Proxy de credibilidad"
    ),
    linewidth = 0.75,
    alpha = 0.95,
    lineend = "round"
  ) +
  
  ggplot2::scale_colour_manual(
    values = c(
      "Inflación observada" = col_inflacion,
      "Expectativa a 24 meses" = col_expectativa,
      "Meta central" = col_meta,
      "Ancla inflacionaria implícita" = col_serie,
      "Proxy de credibilidad" = col_lambda
    ),
    breaks = c(
      "Ancla inflacionaria implícita",
      "Expectativa a 24 meses",
      "Inflación observada",
      "Meta central",
      "Proxy de credibilidad"
    )
  ) +
  
  ggplot2::scale_linetype_manual(
    values = c(
      "Inflación observada" = "solid",
      "Expectativa a 24 meses" = "longdash",
      "Meta central" = "dotted",
      "Ancla inflacionaria implícita" = "solid",
      "Proxy de credibilidad" = "solid"
    ),
    breaks = c(
      "Ancla inflacionaria implícita",
      "Expectativa a 24 meses",
      "Inflación observada",
      "Meta central",
      "Proxy de credibilidad"
    )
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      nrow = 2,
      byrow = TRUE
    ),
    linetype = ggplot2::guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  escala_x_wp +
  
  ggplot2::scale_y_continuous(
    breaks = seq(0, 10, by = 1),
    expand = ggplot2::expansion(mult = c(0, 0)),
    name = "Porcentaje",
    sec.axis = ggplot2::sec_axis(
      trans = ~ lambda_min + (. - y_min) *
        (lambda_max - lambda_min) / (y_max - y_min),
      name = expression(hat(lambda)[t]),
      breaks = seq(0.70, 1.00, by = 0.05)
    )
  ) +
  ggplot2::coord_cartesian(
    ylim = c(0, 10)
  ) +
  
  ggplot2::labs(
    x = NULL,
    colour = NULL,
    linetype = NULL
  ) +
  
  tema_wp +
  
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.spacing.x = grid::unit(0.18, "cm"),
    legend.key.width = grid::unit(1.25, "cm"),
    legend.text = ggplot2::element_text(size = 8.3),
    axis.title.y.right = ggplot2::element_text(
      colour = col_lambda,
      size = 13,
      margin = ggplot2::margin(l = 8)
    ),
    axis.text.y.right = ggplot2::element_text(
      colour = col_lambda
    ),
    axis.line.y.right = ggplot2::element_line(
      colour = col_lambda
    ),
    axis.ticks.y.right = ggplot2::element_line(
      colour = col_lambda
    )
  )

ggplot2::ggsave(
  file.path(
    carpeta_salida,
    "figura_05_dinamica_conjunta_ancla_credibilidad.png"
  ),
  plot = p_sintesis_ancla,
  width = 8.8,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

print(p_sintesis_ancla)
# 14.5 Brecha del ancla implícita respecto de la meta -------------------------

p_brecha <- ggplot2::ggplot(
  trayectoria_M3,
  ggplot2::aes(
    x = fecha,
    y = brecha_ancla_meta
  )
) +
  ggplot2::geom_line(
    colour = col_serie,
    linewidth = 0.80,
    lineend = "round"
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey30",
    linewidth = 0.45
  ) +
  escala_x_wp +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.04, 0.06)
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Puntos porcentuales"
  ) +
  tema_wp

ggplot2::ggsave(
  file.path(
    carpeta_salida,
    "figura_05_brecha_ancla_meta.png"
  ),
  plot = p_brecha,
  width = 8.5,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

print(p_brecha)


# 14.6 ACF de residuos con la misma línea gráfica -----------------------------
#
# Sustituye stats::acf(..., plot = TRUE) por una versión ggplot2
# consistente con el resto del working paper.

crear_grafica_acf_wp <- function(
    residuos,
    nombre_eje_y = "Autocorrelación",
    lag_max = 24
) {
  
  r <- residuos[
    is.finite(residuos)
  ]
  
  acf_obj <- stats::acf(
    r,
    lag.max = lag_max,
    plot = FALSE
  )
  
  datos_acf <- tibble::tibble(
    rezago = as.numeric(
      acf_obj$lag
    ),
    acf = as.numeric(
      acf_obj$acf
    )
  ) |>
    dplyr::filter(
      rezago > 0
    )
  
  limite_95 <- 1.96 / sqrt(
    length(r)
  )
  
  ggplot2::ggplot(
    datos_acf,
    ggplot2::aes(
      x = rezago,
      y = acf
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "grey35",
      linewidth = 0.40
    ) +
    ggplot2::geom_hline(
      yintercept = c(
        -limite_95,
        limite_95
      ),
      linetype = "dashed",
      colour = "grey55",
      linewidth = 0.40
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = rezago,
        y = 0,
        yend = acf
      ),
      colour = col_serie,
      linewidth = 0.65
    ) +
    ggplot2::geom_point(
      colour = col_serie,
      size = 1.45
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(
        0,
        lag_max,
        by = 4
      ),
      minor_breaks = seq(
        0,
        lag_max,
        by = 2
      ),
      expand = ggplot2::expansion(
        mult = c(0.01, 0.02)
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(0.05, 0.08)
      )
    ) +
    ggplot2::labs(
      x = "Rezago",
      y = nombre_eje_y
    ) +
    ggplot2::theme_minimal(
      base_size = 10.5
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        colour = "grey92",
        linewidth = 0.30
      ),
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey88",
        linewidth = 0.35
      ),
      
      axis.title.x = ggplot2::element_text(
        size = 10,
        colour = "black",
        margin = ggplot2::margin(t = 7)
      ),
      axis.title.y = ggplot2::element_text(
        size = 10.5,
        colour = "black",
        margin = ggplot2::margin(r = 8)
      ),
      axis.text = ggplot2::element_text(
        size = 9,
        colour = "grey25"
      ),
      axis.ticks = ggplot2::element_blank(),
      
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      
      plot.margin = ggplot2::margin(
        t = 4,
        r = 8,
        b = 4,
        l = 4
      )
    )
}


p_acf_inflacion <- crear_grafica_acf_wp(
  residuos = residuos_M3[, "inflacion"]
)

ggplot2::ggsave(
  file.path(
    carpeta_salida,
    "figura_06_ACF_residuos_inflacion.png"
  ),
  plot = p_acf_inflacion,
  width = 8.5,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

print(p_acf_inflacion)


p_acf_expectativa <- crear_grafica_acf_wp(
  residuos = residuos_M3[, "expectativa"]
)

ggplot2::ggsave(
  file.path(
    carpeta_salida,
    "figura_07_ACF_residuos_expectativa.png"
  ),
  plot = p_acf_expectativa,
  width = 8.5,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

print(p_acf_expectativa)


# 14.4 Ancla implícita, expectativas, lambda y meta -----------------------------------

p_ancla_lambda <- ggplot2::ggplot(
  trayectoria_M3,
  ggplot2::aes(x = fecha)
) +
  
  # Rango meta oficial: referencia institucional de fondo
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = banda_inferior,
      ymax = banda_superior
    ),
    fill = col_banda,
    alpha = 0.10,
    linewidth = 0
  ) +
  
  # IC 95% del ancla implícita:
  # visible, pero deliberadamente tenue
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = pi_estrella_li_95,
      ymax = pi_estrella_ls_95
    ),
    fill = col_serie,
    alpha = 0.07,
    linewidth = 0
  ) +
  
  # Expectativa a 24 meses:
  # serie secundaria frente al ancla implícita
  ggplot2::geom_line(
    ggplot2::aes(
      y = expectativa,
      colour = "Expectativa a 24 meses",
      linetype = "Expectativa a 24 meses"
    ),
    linewidth = 0.55,
    alpha = 0.65,
    lineend = "round"
  ) +
  
  # Meta central
  ggplot2::geom_line(
    ggplot2::aes(
      y = meta_inflacion,
      colour = "Meta central",
      linetype = "Meta central"
    ),
    linewidth = 0.55,
    alpha = 0.90
  ) +
  
  # Ancla implícita: objeto principal de la figura
  ggplot2::geom_line(
    ggplot2::aes(
      y = pi_estrella_t,
      colour = "Ancla inflacionaria implícita",
      linetype = "Ancla inflacionaria implícita"
    ),
    linewidth = 0.95,
    lineend = "round"
  ) +
  
  ggplot2::scale_colour_manual(
    values = c(
      "Ancla inflacionaria implícita" = col_serie,
      "Expectativa a 24 meses" = col_expectativa,
      "Meta central" = col_meta
    ),
    breaks = c(
      "Ancla inflacionaria implícita",
      "Expectativa a 24 meses",
      "Meta central"
    )
  ) +
  
  ggplot2::scale_linetype_manual(
    values = c(
      "Ancla inflacionaria implícita" = "solid",
      "Expectativa a 24 meses" = "longdash",
      "Meta central" = "dotted"
    ),
    breaks = c(
      "Ancla inflacionaria implícita",
      "Expectativa a 24 meses",
      "Meta central"
    )
  ) +
  
  escala_x_wp +
  
  ggplot2::scale_y_continuous(
    limits = c(2, 7),
    breaks = seq(2, 7, by = 1),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "Porcentaje",
    colour = NULL,
    linetype = NULL
  ) +
  
  tema_wp +
  
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.spacing.x = grid::unit(0.18, "cm"),
    legend.key.width = grid::unit(1.2, "cm"),
    legend.text = ggplot2::element_text(size = 8.5)
  )


ggplot2::ggsave(
  file.path(
    carpeta_salida,
    "figura_04_ancla_implicita.png"
  ),
  plot = p_ancla,
  width = 8.5,
  height = 4.3,
  dpi = 300,
  bg = "white"
)

print(p_ancla)



# 15. GUARDAR OBJETOS ----------------------------------------------------------

saveRDS(
  ajuste_M0$modelo,
  file.path(
    carpeta_salida,
    "modelo_M0.rds"
  )
)

saveRDS(
  ajuste_M3$modelo,
  file.path(
    carpeta_salida,
    "modelo_M3.rds"
  )
)

saveRDS(
  kfs_M3,
  file.path(
    carpeta_salida,
    "kfs_M3.rds"
  )
)


# 16. RESUMEN EJECUTIVO --------------------------------------------------------

Q_M3 <- tabla_hiperparametros |>
  dplyr::filter(
    modelo == "M3"
  )

p_LB_pi_12 <- tabla_diagnosticos |>
  dplyr::filter(
    ecuacion == "Inflación",
    diagnostico == "Ljung-Box",
    rezago == 12
  ) |>
  dplyr::pull(p_valor)

p_LB_e_12 <- tabla_diagnosticos |>
  dplyr::filter(
    ecuacion == "Expectativa 24m",
    diagnostico == "Ljung-Box",
    rezago == 12
  ) |>
  dplyr::pull(p_valor)

fmt <- function(x, digits = 5) {

  if (length(x) == 0 || !is.finite(x[1])) {
    return("NA")
  }

  formatC(
    x[1],
    digits = digits,
    format = "f"
  )
}

fmt_sci <- function(x, digits = 4) {

  if (length(x) == 0 || !is.finite(x[1])) {
    return("NA")
  }

  formatC(
    x[1],
    digits = digits,
    format = "e"
  )
}

lineas_resumen <- c(
  "",
  "======================================================================",
  "DEMERTZIS TVP-VAR RESTRINGIDO M3 - GUATEMALA",
  "======================================================================",
  paste0(
    "Archivo: ",
    archivo_entrada
  ),
  paste0(
    "Muestra efectiva: ",
    format(
      min(datos_modelo$fecha),
      "%Y-%m"
    ),
    " a ",
    format(
      max(datos_modelo$fecha),
      "%Y-%m"
    )
  ),
  paste0(
    "Observaciones: ",
    n_obs
  ),
  "",
  "ESPECIFICACIÓN PRINCIPAL",
  "pi_t = a0 + a*pi_(t-1) + b*e_(t-1) + eps_pi,t",
  "e_t  = c0_t + c_t*pi_(t-1) + d*e_(t-1) + eps_e,t",
  "c0_t y c_t siguen paseos aleatorios.",
  "a0, a, b y d permanecen constantes.",
  "",
  "lambda_t = 1 - c_t/(1-d)",
  "pi_estrella_t = c0_t/(1-d-c_t)",
  "",
  "HIPERPARÁMETROS M3",
  paste0(
    "H inflación: ",
    fmt_sci(Q_M3$H_inflacion)
  ),
  paste0(
    "H expectativas: ",
    fmt_sci(Q_M3$H_expectativa)
  ),
  paste0(
    "Q_c0: ",
    fmt_sci(Q_M3$Q_c0)
  ),
  paste0(
    "Q_c: ",
    fmt_sci(Q_M3$Q_c)
  ),
  "",
  "COMPARACIÓN CON PARÁMETROS CONSTANTES",
  paste0(
    "logLik M0 / M3: ",
    fmt(ll_M0),
    " / ",
    fmt(ll_M3)
  ),
  paste0(
    "AIC M0 / M3: ",
    fmt(AIC_M0),
    " / ",
    fmt(AIC_M3)
  ),
  paste0(
    "BIC M0 / M3: ",
    fmt(BIC_M0),
    " / ",
    fmt(BIC_M3)
  ),
  paste0(
    "LR M3 vs M0: ",
    fmt(LR_referencial),
    " | p chi2(2) SOLO REFERENCIAL: ",
    fmt_sci(p_LR_chi2_referencial)
  ),
  "Nota: las Q están en la frontera bajo H0; no interpretar el p chi2 como prueba definitiva.",
  "",
  "PROXY DE CREDIBILIDAD",
  paste0(
    "lambda promedio: ",
    fmt(resumen_M3$lambda_promedio)
  ),
  paste0(
    "lambda inicial / final: ",
    fmt(resumen_M3$lambda_inicial),
    " / ",
    fmt(resumen_M3$lambda_final)
  ),
  paste0(
    "lambda mínimo: ",
    fmt(resumen_M3$lambda_minimo),
    " en ",
    format(
      resumen_M3$fecha_lambda_minimo,
      "%Y-%m"
    )
  ),
  paste0(
    "lambda máximo: ",
    fmt(resumen_M3$lambda_maximo),
    " en ",
    format(
      resumen_M3$fecha_lambda_maximo,
      "%Y-%m"
    )
  ),
  paste0(
    "% meses lambda fuera de [0,1]: ",
    fmt(
      resumen_M3$porcentaje_lambda_fuera_0_1,
      2
    ),
    "%"
  ),
  "",
  "ANCLA IMPLÍCITA",
  paste0(
    "pi estrella promedio: ",
    fmt(
      resumen_M3$pi_estrella_promedio
    )
  ),
  paste0(
    "pi estrella inicial / final: ",
    fmt(
      resumen_M3$pi_estrella_inicial
    ),
    " / ",
    fmt(
      resumen_M3$pi_estrella_final
    )
  ),
  paste0(
    "MAE pi estrella vs meta: ",
    fmt(
      resumen_M3$MAE_pi_estrella_meta
    )
  ),
  paste0(
    "% meses ancla dentro de banda meta: ",
    fmt(
      resumen_M3$porcentaje_ancla_dentro_banda_meta,
      2
    ),
    "%"
  ),
  "",
  "ESTABILIDAD LOCAL",
  paste0(
    "% meses con radio espectral < 1: ",
    fmt(
      resumen_M3$porcentaje_meses_estables,
      2
    ),
    "%"
  ),
  paste0(
    "Radio espectral máximo: ",
    fmt(
      resumen_M3$radio_espectral _maximo
    )
  ),
  "",
  "DIAGNÓSTICOS RESIDUALES M3",
  paste0(
    "Ljung-Box inflación, lag 12, p: ",
    fmt_sci(
      p_LB_pi_12
    )
  ),
  paste0(
    "Ljung-Box expectativas, lag 12, p: ",
    fmt_sci(
      p_LB_e_12
    )
  ),
  paste0(
    "Correlación contemporánea de residuos: ",
    fmt(
      rho_residuos
    ),
    " | p Fisher: ",
    fmt_sci(
      p_fisher
    )
  ),
  "",
  "INTERPRETACIÓN",
  "1. lambda_t resume el desacoplamiento de expectativas respecto de la inflación.",
  "2. pi_estrella_t estima el nivel del ancla implícita.",
  "3. Para un régimen creíble bajo metas de inflación interesa tanto un lambda_t elevado",
  "   como un ancla implícita próxima a la meta oficial.",
  "4. La expectativa a 24 meses se interpreta como anclaje de mediano plazo.",
  "",
  paste0(
    "Resultados guardados en: ",
    carpeta_salida
  ),
  "======================================================================"
)

writeLines(
  lineas_resumen,
  file.path(
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


