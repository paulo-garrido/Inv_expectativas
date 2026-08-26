# ==============================================================================
# OBJETIVO 2 - ROBUSTEZ: TVP-VAR RESTRINGIDO EN PRIMERAS DIFERENCIAS
# Guatemala: inflación y expectativas de inflación a 24 meses
# ==============================================================================
#
# PROPÓSITO
# ---------
# Reestimar la MISMA arquitectura bivariada del modelo TVP-VAR central,
# pero utilizando inflación y expectativas a 24 meses en primeras diferencias.
#
# La finalidad es evaluar si la reducción de la sensibilidad de las
# expectativas frente a la inflación se mantiene cuando se elimina la elevada
# persistencia contenida en los niveles de las series.
#
#
# MODELO DE ROBUSTEZ
# ------------------
#
# Inflación:
#
#   Δpi_t =
#       a0
#       + a * Δpi_(t-1)
#       + b * Δe_(t-1)
#       + eps_pi,t
#
#
# Expectativas a 24 meses:
#
#   Δe_t =
#       c0_t
#       + c_t * Δpi_(t-1)
#       + d * Δe_(t-1)
#       + eps_e,t
#
#
# Estados variantes:
#
#   c0_t = c0_(t-1) + eta_c0,t
#
#   c_t  = c_(t-1)  + eta_c,t
#
#
# Parámetros constantes:
#
#   a0, a, b, d
#
#
# INTERPRETACIÓN
# --------------
#
# c_t mide la sensibilidad cambiante de los CAMBIOS en las expectativas
# frente a los CAMBIOS rezagados en la inflación.
#
# Una reducción de c_t a lo largo de la muestra constituye evidencia
# complementaria de una menor sensibilidad de las expectativas.
#
#
# IMPORTANTE
# ----------
#
# En esta especificación NO se calculan:
#
#   lambda_t
#   pi_estrella_t
#
# porque ambas magnitudes se derivan de una relación de equilibrio en niveles
# y no deben trasladarse mecánicamente a una representación en diferencias.
#
#
# MODELOS ESTIMADOS
# -----------------
#
# M0_D:
#   VAR(1) en primeras diferencias con parámetros constantes.
#
# M3_D:
#   TVP-VAR(1) restringido en primeras diferencias:
#   únicamente c0_t y c_t varían en el tiempo.
#
#
# SALIDA
# ------
#
# Carpeta:
#
#   resultados_objetivo2_TVP_VAR_primeras_diferencias
#
# Archivos:
#
#   00_resumen_ejecutivo.txt
#   01_datos_diferencias.csv
#   02_coeficientes_OLS_inicializacion.csv
#   03_intentos_estimacion.csv
#   04_hiperparametros_M0_M3.csv
#   05_comparacion_M0_M3.csv
#   06_trayectoria_TVP_VAR_diferencias.csv
#   07_resumen_c_t.csv
#   08_resumen_subperiodos_c_t.csv
#   09_diagnosticos_residuales.csv
#   10_estabilidad_local.csv
#   11_residuos_estandarizados.csv
#
# Figuras:
#
#   figura_01_c_t_diferencias.png
#   figura_02_c0_t_diferencias.png
#
# Objetos:
#
#   modelo_M0_D.rds
#   modelo_M3_D.rds
#   kfs_M3_D.rds
#
# ==============================================================================


# ==============================================================================
# 0. DIRECTORIO, LIMPIEZA Y PAQUETES
# ==============================================================================

setwd(
    "/Users/paulogarridogrijalva/Documents/PES/Seminario/Inv_expectativas/Data y resultados"
)

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
    !vapply(
        paquetes,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]


if (length(faltantes) > 0) {
    
    install.packages(
        faltantes,
        dependencies = TRUE
    )
}


invisible(
    lapply(
        paquetes,
        library,
        character.only = TRUE
    )
)


cat(
    "\nDirectorio de trabajo:\n",
    getwd(),
    "\n"
)


# ==============================================================================
# 1. CONFIGURACIÓN
# ==============================================================================

archivo_entrada <-
    "bd_sample_modelo.csv"


if (!file.exists(archivo_entrada)) {
    
    stop(
        "No se encontró 'bd_sample_modelo.csv' en:\n",
        getwd()
    )
}


carpeta_salida <- file.path(
    getwd(),
    "resultados_objetivo2_TVP_VAR_primeras_diferencias"
)


dir.create(
    carpeta_salida,
    recursive = TRUE,
    showWarnings = FALSE
)


fecha_fin_muestra <-
    as.Date(
        "2026-06-01"
    )


columnas_requeridas <- c(
    "fecha",
    "infl_gt",
    "exp_inf_24m"
)


nivel_confianza <- 0.95


z_critico <- qnorm(
    1 -
        (
            1 -
                nivel_confianza
        ) / 2
)


# Límite inferior de las varianzas estimadas.
limite_varianza <- 1e-12


# Diagnósticos residuales.
rezagos_ljung_box <- c(
    6,
    12,
    24
)


# Iteraciones máximas.
maxit_M0 <- 20000
maxit_M3 <- 30000


# ==============================================================================
# 2. FUNCIONES AUXILIARES
# ==============================================================================


verificar_columnas <- function(
        datos,
        columnas
) {
    
    faltan <- setdiff(
        columnas,
        names(datos)
    )
    
    
    if (length(faltan) > 0) {
        
        stop(
            "Faltan las siguientes variables requeridas: ",
            paste(
                faltan,
                collapse = ", "
            )
        )
    }
}


# ------------------------------------------------------------------------------
# Jarque-Bera manual
# ------------------------------------------------------------------------------

jarque_bera_manual <- function(x) {
    
    x <- x[
        is.finite(x)
    ]
    
    
    n <- length(x)
    
    
    if (
        n < 8 ||
        !is.finite(
            stats::sd(x)
        ) ||
        stats::sd(x) == 0
    ) {
        
        return(
            tibble::tibble(
                estadistico = NA_real_,
                p_valor = NA_real_
            )
        )
    }
    
    
    z <- (
        x -
            mean(x)
    ) /
        stats::sd(x)
    
    
    asimetria <-
        mean(
            z^3
        )
    
    
    curtosis <-
        mean(
            z^4
        )
    
    
    JB <- n / 6 *
        (
            asimetria^2 +
                (
                    curtosis - 3
                )^2 / 4
        )
    
    
    tibble::tibble(
        estadistico = JB,
        p_valor =
            stats::pchisq(
                JB,
                df = 2,
                lower.tail = FALSE
            )
    )
}


# ------------------------------------------------------------------------------
# Correlación segura
# ------------------------------------------------------------------------------

correlacion_segura <- function(
        x,
        y
) {
    
    validos <-
        is.finite(x) &
        is.finite(y)
    
    
    if (
        sum(validos) <
        4
    ) {
        
        return(
            NA_real_
        )
    }
    
    
    x <- x[
        validos
    ]
    
    y <- y[
        validos
    ]
    
    
    if (
        stats::sd(x) == 0 ||
        stats::sd(y) == 0
    ) {
        
        return(
            NA_real_
        )
    }
    
    
    stats::cor(
        x,
        y
    )
}


# ==============================================================================
# 3. CARGA Y PREPARACIÓN DE LOS DATOS
# ==============================================================================

bd_original <- readr::read_csv(
    archivo_entrada,
    show_col_types = FALSE,
    na = c(
        "",
        "NA",
        "N/A",
        ".",
        "null"
    )
)


verificar_columnas(
    bd_original,
    columnas_requeridas
)


bd <- bd_original |>
    
    dplyr::mutate(
        fecha =
            as.Date(
                fecha
            )
    ) |>
    
    dplyr::arrange(
        fecha
    ) |>
    
    dplyr::filter(
        fecha <=
            fecha_fin_muestra
    )


if (anyNA(bd$fecha)) {
    
    stop(
        "Existen fechas que no pudieron convertirse al formato Date."
    )
}


if (anyDuplicated(bd$fecha) > 0) {
    
    stop(
        "Existen fechas duplicadas en la base."
    )
}


# ------------------------------------------------------------------------------
# Selección de las variables en niveles
# ------------------------------------------------------------------------------

datos_niveles <- bd |>
    
    dplyr::transmute(
        
        fecha,
        
        inflacion =
            infl_gt,
        
        expectativa =
            exp_inf_24m
    ) |>
    
    dplyr::filter(
        is.finite(
            inflacion
        ),
        is.finite(
            expectativa
        )
    )


if (
    nrow(
        datos_niveles
    ) < 80
) {
    
    stop(
        "La muestra efectiva contiene menos de 80 observaciones."
    )
}


# ------------------------------------------------------------------------------
# Verificación de continuidad mensual
# ------------------------------------------------------------------------------

indice_mensual <-
    lubridate::year(
        datos_niveles$fecha
    ) * 12 +
    lubridate::month(
        datos_niveles$fecha
    )


if (
    any(
        diff(
            indice_mensual
        ) != 1
    )
) {
    
    stop(
        "La muestra presenta saltos mensuales. ",
        "Revise observaciones faltantes en inflación o expectativas."
    )
}


# ------------------------------------------------------------------------------
# Primeras diferencias
# ------------------------------------------------------------------------------

datos_diferencias <- datos_niveles |>
    
    dplyr::mutate(
        
        d_inflacion =
            inflacion -
            dplyr::lag(
                inflacion
            ),
        
        d_expectativa =
            expectativa -
            dplyr::lag(
                expectativa
            )
    ) |>
    
    dplyr::filter(
        is.finite(
            d_inflacion
        ),
        is.finite(
            d_expectativa
        )
    )


# ------------------------------------------------------------------------------
# VAR(1) sobre las primeras diferencias
# ------------------------------------------------------------------------------

datos_modelo <- datos_diferencias |>
    
    dplyr::mutate(
        
        d_inflacion_lag1 =
            dplyr::lag(
                d_inflacion
            ),
        
        d_expectativa_lag1 =
            dplyr::lag(
                d_expectativa
            )
    ) |>
    
    dplyr::filter(
        is.finite(
            d_inflacion_lag1
        ),
        is.finite(
            d_expectativa_lag1
        )
    )


Y_modelo <- as.matrix(
    
    datos_modelo |>
        
        dplyr::select(
            d_inflacion,
            d_expectativa
        )
)


n_obs <-
    nrow(
        datos_modelo
    )


if (
    n_obs < 60
) {
    
    stop(
        "La muestra del TVP-VAR en diferencias contiene menos de 60 observaciones."
    )
}


readr::write_csv(
    datos_modelo,
    file.path(
        carpeta_salida,
        "01_datos_diferencias.csv"
    )
)


cat(
    "\nMuestra efectiva del TVP-VAR en diferencias:\n",
    format(
        min(
            datos_modelo$fecha
        ),
        "%Y-%m"
    ),
    " a ",
    format(
        max(
            datos_modelo$fecha
        ),
        "%Y-%m"
    ),
    "\n"
)


cat(
    "Número de observaciones:",
    n_obs,
    "\n"
)


# ==============================================================================
# 4. VAR(1) OLS PARA INICIALIZACIÓN
# ==============================================================================

inicializar_desde_OLS <- function(datos) {
    
    
    # --------------------------------------------------------------------------
    # Ecuación de inflación
    # --------------------------------------------------------------------------
    
    eq_pi <- stats::lm(
        
        d_inflacion ~
            d_inflacion_lag1 +
            d_expectativa_lag1,
        
        data = datos
    )
    
    
    # --------------------------------------------------------------------------
    # Ecuación de expectativas
    # --------------------------------------------------------------------------
    
    eq_e <- stats::lm(
        
        d_expectativa ~
            d_inflacion_lag1 +
            d_expectativa_lag1,
        
        data = datos
    )
    
    
    beta_pi <-
        stats::coef(
            eq_pi
        )
    
    
    beta_e <-
        stats::coef(
            eq_e
        )
    
    
    estado_inicial <- c(
        
        a0 =
            unname(
                beta_pi[
                    "(Intercept)"
                ]
            ),
        
        a =
            unname(
                beta_pi[
                    "d_inflacion_lag1"
                ]
            ),
        
        b =
            unname(
                beta_pi[
                    "d_expectativa_lag1"
                ]
            ),
        
        c0 =
            unname(
                beta_e[
                    "(Intercept)"
                ]
            ),
        
        c =
            unname(
                beta_e[
                    "d_inflacion_lag1"
                ]
            ),
        
        d =
            unname(
                beta_e[
                    "d_expectativa_lag1"
                ]
            )
    )
    
    
    if (
        any(
            !is.finite(
                estado_inicial
            )
        )
    ) {
        
        stop(
            "La inicialización OLS produjo coeficientes no finitos. ",
            "Revise la base y posibles problemas de colinealidad."
        )
    }
    
    
    H_inicial <- pmax(
        
        c(
            
            stats::var(
                stats::residuals(
                    eq_pi
                )
            ),
            
            stats::var(
                stats::residuals(
                    eq_e
                )
            )
        ),
        
        1e-8
    )
    
    
    list(
        
        estado_inicial =
            estado_inicial,
        
        H_inicial =
            H_inicial,
        
        eq_pi =
            eq_pi,
        
        eq_e =
            eq_e
    )
}


init <-
    inicializar_desde_OLS(
        datos_modelo
    )


estado_inicial <-
    init$estado_inicial


H_inicial <-
    init$H_inicial


tabla_OLS <- dplyr::bind_rows(
    
    tibble::tibble(
        
        ecuacion =
            "Δ Inflación",
        
        parametro =
            names(
                stats::coef(
                    init$eq_pi
                )
            ),
        
        estimacion =
            as.numeric(
                stats::coef(
                    init$eq_pi
                )
            )
    ),
    
    
    tibble::tibble(
        
        ecuacion =
            "Δ Expectativa 24m",
        
        parametro =
            names(
                stats::coef(
                    init$eq_e
                )
            ),
        
        estimacion =
            as.numeric(
                stats::coef(
                    init$eq_e
                )
            )
    )
)


readr::write_csv(
    tabla_OLS,
    file.path(
        carpeta_salida,
        "02_coeficientes_OLS_inicializacion.csv"
    )
)


# ==============================================================================
# 5. CONSTRUCCIÓN DEL MODELO ESTADO-ESPACIO
# ==============================================================================
#
# Estados:
#
#   1 = a0
#   2 = a
#   3 = b
#   4 = c0_t
#   5 = c_t
#   6 = d
#
#
# M0_D:
#
#   Q = 0 para todos los estados.
#
#
# M3_D:
#
#   Q_c0 > 0
#   Q_c  > 0
#
#   Los demás elementos de Q permanecen en cero.
#
# ==============================================================================


crear_modelo <- function(
        datos,
        Y,
        estado_inicial,
        modelo = c(
            "M0_D",
            "M3_D"
        )
) {
    
    
    modelo <-
        match.arg(
            modelo
        )
    
    
    n <-
        nrow(
            datos
        )
    
    
    m <- 6L
    
    
    # --------------------------------------------------------------------------
    # Matriz Z_t
    #
    # Dimensiones:
    #
    #   2 variables observadas
    #   x 6 estados
    #   x n períodos
    # --------------------------------------------------------------------------
    
    Z_t <- array(
        0,
        dim = c(
            2,
            m,
            n
        )
    )
    
    
    # --------------------------------------------------------------------------
    # Ecuación 1:
    #
    # Δpi_t =
    #   a0
    #   + a Δpi_(t-1)
    #   + b Δe_(t-1)
    #   + eps_pi,t
    # --------------------------------------------------------------------------
    
    Z_t[
        1,
        1,
    ] <- 1
    
    
    Z_t[
        1,
        2,
    ] <-
        datos$d_inflacion_lag1
    
    
    Z_t[
        1,
        3,
    ] <-
        datos$d_expectativa_lag1
    
    
    # --------------------------------------------------------------------------
    # Ecuación 2:
    #
    # Δe_t =
    #   c0_t
    #   + c_t Δpi_(t-1)
    #   + d Δe_(t-1)
    #   + eps_e,t
    # --------------------------------------------------------------------------
    
    Z_t[
        2,
        4,
    ] <- 1
    
    
    Z_t[
        2,
        5,
    ] <-
        datos$d_inflacion_lag1
    
    
    Z_t[
        2,
        6,
    ] <-
        datos$d_expectativa_lag1
    
    
    # --------------------------------------------------------------------------
    # Matriz de varianzas de transición
    # --------------------------------------------------------------------------
    
    Q_mat <- matrix(
        0,
        nrow = m,
        ncol = m
    )
    
    
    if (
        modelo ==
        "M3_D"
    ) {
        
        Q_mat[
            4,
            4
        ] <-
            NA_real_
        
        
        Q_mat[
            5,
            5
        ] <-
            NA_real_
    }
    
    
    # --------------------------------------------------------------------------
    # Construcción KFAS
    # --------------------------------------------------------------------------
    
    KFAS::SSModel(
        
        Y ~ -1 +
            
            SSMcustom(
                
                Z =
                    Z_t,
                
                T =
                    diag(
                        m
                    ),
                
                R =
                    diag(
                        m
                    ),
                
                Q =
                    Q_mat,
                
                a1 =
                    estado_inicial,
                
                P1 =
                    matrix(
                        0,
                        m,
                        m
                    ),
                
                P1inf =
                    diag(
                        m
                    ),
                
                n =
                    n
            ),
        
        H =
            diag(
                rep(
                    NA_real_,
                    2
                )
            )
    )
}


# ==============================================================================
# 6. FUNCIONES PARA ESTIMAR H Y Q
# ==============================================================================


crear_updatefn <- function(
        modelo = c(
            "M0_D",
            "M3_D"
        )
) {
    
    
    modelo <-
        match.arg(
            modelo
        )
    
    
    function(
        par,
        objeto
    ) {
        
        
        # ------------------------------------------------------------------------
        # Varianzas del error de observación
        # ------------------------------------------------------------------------
        
        objeto$H[
            1,
            1,
            1
        ] <-
            exp(
                par[1]
            )
        
        
        objeto$H[
            2,
            2,
            1
        ] <-
            exp(
                par[2]
            )
        
        
        # Se mantiene H diagonal.
        
        objeto$H[
            1,
            2,
            1
        ] <- 0
        
        
        objeto$H[
            2,
            1,
            1
        ] <- 0
        
        
        # ------------------------------------------------------------------------
        # M3_D:
        #
        # Q_c0 y Q_c son estimadas.
        # ------------------------------------------------------------------------
        
        if (
            modelo ==
            "M3_D"
        ) {
            
            objeto$Q[
                4,
                4,
                1
            ] <-
                exp(
                    par[3]
                )
            
            
            objeto$Q[
                5,
                5,
                1
            ] <-
                exp(
                    par[4]
                )
        }
        
        
        objeto
    }
}


# ------------------------------------------------------------------------------
# Valores iniciales del optimizador
# ------------------------------------------------------------------------------

crear_inicios <- function(
        H_inicial,
        modelo = c(
            "M0_D",
            "M3_D"
        )
) {
    
    
    modelo <-
        match.arg(
            modelo
        )
    
    
    if (
        modelo ==
        "M0_D"
    ) {
        
        return(
            
            list(
                
                log(
                    H_inicial *
                        0.5
                ),
                
                log(
                    H_inicial
                ),
                
                log(
                    H_inicial *
                        2
                )
            )
        )
    }
    
    
    # Distintos órdenes de magnitud para Q_c0 y Q_c.
    
    Q_inicios <- list(
        
        c(
            1e-10,
            1e-10
        ),
        
        c(
            1e-8,
            1e-8
        ),
        
        c(
            1e-6,
            1e-6
        ),
        
        c(
            1e-5,
            1e-5
        ),
        
        c(
            1e-4,
            1e-4
        ),
        
        c(
            1e-3,
            1e-3
        ),
        
        c(
            1e-3,
            1e-6
        ),
        
        c(
            1e-6,
            1e-3
        ),
        
        c(
            1e-4,
            1e-6
        ),
        
        c(
            1e-6,
            1e-4
        )
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
    
    
    # Inicios adicionales modificando H.
    
    inicios <- c(
        
        inicios,
        
        list(
            
            log(
                c(
                    H_inicial *
                        0.5,
                    1e-5,
                    1e-5
                )
            ),
            
            log(
                c(
                    H_inicial *
                        2,
                    1e-5,
                    1e-5
                )
            )
        )
    )
    
    
    inicios
}


# ==============================================================================
# 7. ESTIMACIÓN CON MÚLTIPLES INICIOS
# ==============================================================================


estimar_modelo <- function(
        codigo,
        datos,
        Y,
        estado_inicial,
        H_inicial,
        maxit
) {
    
    
    modelo_base <- crear_modelo(
        
        datos =
            datos,
        
        Y =
            Y,
        
        estado_inicial =
            estado_inicial,
        
        modelo =
            codigo
    )
    
    
    updatefn <-
        crear_updatefn(
            codigo
        )
    
    
    inicios <-
        crear_inicios(
            H_inicial,
            codigo
        )
    
    
    numero_parametros <-
        
        if (
            codigo ==
            "M0_D"
        ) {
            
            2L
            
        } else {
            
            4L
        }
    
    
    resultados <- lapply(
        
        seq_along(
            inicios
        ),
        
        function(i) {
            
            
            fit <- tryCatch(
                
                KFAS::fitSSM(
                    
                    model =
                        modelo_base,
                    
                    inits =
                        inicios[[i]],
                    
                    updatefn =
                        updatefn,
                    
                    method =
                        "L-BFGS-B",
                    
                    lower =
                        rep(
                            log(
                                limite_varianza
                            ),
                            numero_parametros
                        ),
                    
                    upper =
                        rep(
                            log(
                                100
                            ),
                            numero_parametros
                        ),
                    
                    control =
                        list(
                            maxit =
                                maxit,
                            factr =
                                1e7
                        )
                ),
                
                error =
                    function(e) {
                        
                        NULL
                    }
            )
            
            
            if (
                is.null(
                    fit
                )
            ) {
                
                return(
                    
                    list(
                        
                        intento =
                            i,
                        
                        convergencia =
                            NA_integer_,
                        
                        logLik =
                            NA_real_,
                        
                        fit =
                            NULL
                    )
                )
            }
            
            
            ll <- tryCatch(
                
                as.numeric(
                    stats::logLik(
                        fit$model
                    )
                ),
                
                error =
                    function(e) {
                        
                        NA_real_
                    }
            )
            
            
            list(
                
                intento =
                    i,
                
                convergencia =
                    fit$optim.out$convergence,
                
                logLik =
                    ll,
                
                fit =
                    fit
            )
        }
    )
    
    
    tabla_intentos <- dplyr::bind_rows(
        
        lapply(
            
            resultados,
            
            function(x) {
                
                tibble::tibble(
                    
                    modelo =
                        codigo,
                    
                    intento =
                        x$intento,
                    
                    convergencia =
                        x$convergencia,
                    
                    log_verosimilitud =
                        x$logLik
                )
            }
        )
    )
    
    
    validos <- resultados[
        
        vapply(
            
            resultados,
            
            function(x) {
                
                !is.null(
                    x$fit
                ) &&
                    is.finite(
                        x$logLik
                    )
            },
            
            logical(1)
        )
    ]
    
    
    if (
        length(
            validos
        ) == 0
    ) {
        
        stop(
            "No fue posible estimar ",
            codigo,
            "."
        )
    }
    
    
    convergentes <- validos[
        
        vapply(
            
            validos,
            
            function(x) {
                
                isTRUE(
                    x$convergencia ==
                        0
                )
            },
            
            logical(1)
        )
    ]
    
    
    candidatos <-
        
        if (
            length(
                convergentes
            ) > 0
        ) {
            
            convergentes
            
        } else {
            
            warning(
                "No hubo convergencia formal en ",
                codigo,
                "; se conservará el intento con mayor log-verosimilitud."
            )
            
            validos
        }
    
    
    mejor <- candidatos[[
        
        which.max(
            
            vapply(
                
                candidatos,
                
                function(x) {
                    
                    x$logLik
                },
                
                numeric(1)
            )
        )
    ]]
    
    
    list(
        
        codigo =
            codigo,
        
        fit =
            mejor$fit,
        
        modelo =
            mejor$fit$model,
        
        intento =
            mejor$intento,
        
        convergencia =
            mejor$convergencia,
        
        logLik =
            mejor$logLik,
        
        tabla_intentos =
            tabla_intentos
    )
}


# ==============================================================================
# 8. ESTIMACIÓN DE M0_D Y M3_D
# ==============================================================================

cat(
    "\nEstimando M0_D: VAR constante en primeras diferencias...\n"
)


ajuste_M0 <- estimar_modelo(
    
    codigo =
        "M0_D",
    
    datos =
        datos_modelo,
    
    Y =
        Y_modelo,
    
    estado_inicial =
        estado_inicial,
    
    H_inicial =
        H_inicial,
    
    maxit =
        maxit_M0
)


cat(
    "\nEstimando M3_D: TVP-VAR restringido en primeras diferencias...\n"
)


ajuste_M3 <- estimar_modelo(
    
    codigo =
        "M3_D",
    
    datos =
        datos_modelo,
    
    Y =
        Y_modelo,
    
    estado_inicial =
        estado_inicial,
    
    H_inicial =
        H_inicial,
    
    maxit =
        maxit_M3
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


# ==============================================================================
# 9. FILTRO Y SUAVIZAMIENTO
# ==============================================================================

kfs_M0 <- KFAS::KFS(
    
    ajuste_M0$modelo,
    
    filtering =
        c(
            "state",
            "mean"
        ),
    
    smoothing =
        c(
            "state",
            "mean",
            "disturbance"
        )
)


kfs_M3 <- KFAS::KFS(
    
    ajuste_M3$modelo,
    
    filtering =
        c(
            "state",
            "mean"
        ),
    
    smoothing =
        c(
            "state",
            "mean",
            "disturbance"
        )
)


# ==============================================================================
# 10. HIPERPARÁMETROS
# ==============================================================================

extraer_hiperparametros <- function(
        codigo,
        ajuste
) {
    
    
    H <- diag(
        ajuste$modelo$H[
            ,
            ,
            1
        ]
    )
    
    
    Q <- diag(
        ajuste$modelo$Q[
            ,
            ,
            1
        ]
    )
    
    
    tibble::tibble(
        
        modelo =
            codigo,
        
        H_d_inflacion =
            H[1],
        
        H_d_expectativa =
            H[2],
        
        Q_a0 =
            Q[1],
        
        Q_a =
            Q[2],
        
        Q_b =
            Q[3],
        
        Q_c0 =
            Q[4],
        
        Q_c =
            Q[5],
        
        Q_d =
            Q[6],
        
        convergencia =
            ajuste$convergencia,
        
        intento_seleccionado =
            ajuste$intento,
        
        log_verosimilitud =
            ajuste$logLik
    )
}


tabla_hiperparametros <- dplyr::bind_rows(
    
    extraer_hiperparametros(
        "M0_D",
        ajuste_M0
    ),
    
    extraer_hiperparametros(
        "M3_D",
        ajuste_M3
    )
)


readr::write_csv(
    tabla_hiperparametros,
    file.path(
        carpeta_salida,
        "04_hiperparametros_M0_M3.csv"
    )
)


# ==============================================================================
# 11. COMPARACIÓN M0_D VS M3_D
# ==============================================================================

# Se replica el mismo conteo utilizado en el modelo central:
#
# M0_D:
#   6 coeficientes + 2 varianzas H = 8
#
# M3_D:
#   6 coeficientes + 2 H + 2 Q = 10


k_M0 <- 8
k_M3 <- 10


ll_M0 <-
    ajuste_M0$logLik


ll_M3 <-
    ajuste_M3$logLik


AIC_M0 <-
    -2 *
    ll_M0 +
    2 *
    k_M0


AIC_M3 <-
    -2 *
    ll_M3 +
    2 *
    k_M3


BIC_M0 <-
    -2 *
    ll_M0 +
    log(
        n_obs
    ) *
    k_M0


BIC_M3 <-
    -2 *
    ll_M3 +
    log(
        n_obs
    ) *
    k_M3


tabla_comparacion <- tibble::tibble(
    
    modelo =
        c(
            "VAR diferencias",
            "TVP-VAR diferencias"
        ),
    
    descripcion =
        c(
            "VAR(1) constante en primeras diferencias",
            "TVP-VAR(1) restringido en primeras diferencias"
        ),
    
    log_verosimilitud =
        c(
            ll_M0,
            ll_M3
        ),
    
    parametros_contados =
        c(
            k_M0,
            k_M3
        ),
    
    AIC =
        c(
            AIC_M0,
            AIC_M3
        ),
    
    BIC =
        c(
            BIC_M0,
            BIC_M3
        )
)


readr::write_csv(
    tabla_comparacion,
    file.path(
        carpeta_salida,
        "05_comparacion_M0_M3.csv"
    )
)


# ==============================================================================
# 12. EXTRACCIÓN DE LOS ESTADOS SUAVIZADOS
# ==============================================================================

nombres_estados <- c(
    "a0_t",
    "a_t",
    "b_t",
    "c0_t",
    "c_t",
    "d_t"
)


estados_M3 <- as.data.frame(
    kfs_M3$alphahat
)


names(
    estados_M3
) <-
    nombres_estados


trayectoria_M3 <- dplyr::bind_cols(
    
    datos_modelo |>
        
        dplyr::select(
            fecha,
            d_inflacion,
            d_expectativa,
            d_inflacion_lag1,
            d_expectativa_lag1
        ),
    
    estados_M3
)


# ------------------------------------------------------------------------------
# Errores estándar de los dos estados variantes
# ------------------------------------------------------------------------------

trayectoria_M3$c0_se <- sqrt(
    
    pmax(
        kfs_M3$V[
            4,
            4,
        ],
        0
    )
)


trayectoria_M3$c_se <- sqrt(
    
    pmax(
        kfs_M3$V[
            5,
            5,
        ],
        0
    )
)


trayectoria_M3 <- trayectoria_M3 |>
    
    dplyr::mutate(
        
        c0_li_95 =
            c0_t -
            z_critico *
            c0_se,
        
        c0_ls_95 =
            c0_t +
            z_critico *
            c0_se,
        
        c_li_95 =
            c_t -
            z_critico *
            c_se,
        
        c_ls_95 =
            c_t +
            z_critico *
            c_se,
        
        c_significativo_5pct =
            c_li_95 >
            0 |
            c_ls_95 <
            0
    )


# ==============================================================================
# 13. ESTABILIDAD LOCAL DEL TVP-VAR EN DIFERENCIAS
# ==============================================================================

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
        
        fecha =
            trayectoria_M3$fecha[t],
        
        raiz_1_real =
            Re(
                raices[1]
            ),
        
        raiz_1_imaginaria =
            Im(
                raices[1]
            ),
        
        raiz_2_real =
            Re(
                raices[2]
            ),
        
        raiz_2_imaginaria =
            Im(
                raices[2]
            ),
        
        radio_espectral =
            max(
                Mod(
                    raices
                )
            ),
        
        estable =
            max(
                Mod(
                    raices
                )
            ) <
            1
    )
}


tabla_estabilidad <- dplyr::bind_rows(
    
    lapply(
        
        seq_len(
            nrow(
                trayectoria_M3
            )
        ),
        
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
        
        by =
            "fecha"
    )


readr::write_csv(
    trayectoria_M3,
    file.path(
        carpeta_salida,
        "06_trayectoria_TVP_VAR_diferencias.csv"
    )
)


readr::write_csv(
    tabla_estabilidad,
    file.path(
        carpeta_salida,
        "10_estabilidad_local.csv"
    )
)


# ==============================================================================
# 14. RESUMEN GENERAL DE c_t
# ==============================================================================

indice_min_c <-
    which.min(
        trayectoria_M3$c_t
    )


indice_max_c <-
    which.max(
        trayectoria_M3$c_t
    )


resumen_c <- tibble::tibble(
    
    fecha_inicial =
        dplyr::first(
            trayectoria_M3$fecha
        ),
    
    fecha_final =
        dplyr::last(
            trayectoria_M3$fecha
        ),
    
    observaciones =
        nrow(
            trayectoria_M3
        ),
    
    c_promedio =
        mean(
            trayectoria_M3$c_t,
            na.rm = TRUE
        ),
    
    c_inicial =
        dplyr::first(
            trayectoria_M3$c_t
        ),
    
    c_final =
        dplyr::last(
            trayectoria_M3$c_t
        ),
    
    c_inicial_li_95 =
        dplyr::first(
            trayectoria_M3$c_li_95
        ),
    
    c_inicial_ls_95 =
        dplyr::first(
            trayectoria_M3$c_ls_95
        ),
    
    c_final_li_95 =
        dplyr::last(
            trayectoria_M3$c_li_95
        ),
    
    c_final_ls_95 =
        dplyr::last(
            trayectoria_M3$c_ls_95
        ),
    
    c_minimo =
        trayectoria_M3$c_t[
            indice_min_c
        ],
    
    fecha_c_minimo =
        trayectoria_M3$fecha[
            indice_min_c
        ],
    
    c_maximo =
        trayectoria_M3$c_t[
            indice_max_c
        ],
    
    fecha_c_maximo =
        trayectoria_M3$fecha[
            indice_max_c
        ],
    
    porcentaje_c_significativo_5pct =
        100 *
        mean(
            trayectoria_M3$c_significativo_5pct,
            na.rm = TRUE
        ),
    
    porcentaje_meses_estables =
        100 *
        mean(
            trayectoria_M3$estable,
            na.rm = TRUE
        ),
    
    radio_espectral_maximo =
        max(
            trayectoria_M3$radio_espectral,
            na.rm = TRUE
        )
)


readr::write_csv(
    resumen_c,
    file.path(
        carpeta_salida,
        "07_resumen_c_t.csv"
    )
)


# ==============================================================================
# 15. RESUMEN POR SUBPERÍODOS
# ==============================================================================

resumen_subperiodos <- trayectoria_M3 |>
    
    dplyr::mutate(
        
        subperiodo =
            dplyr::case_when(
                
                fecha <=
                    as.Date(
                        "2014-12-01"
                    ) ~
                    "2010-2014",
                
                fecha <=
                    as.Date(
                        "2019-12-01"
                    ) ~
                    "2015-2019",
                
                fecha <=
                    as.Date(
                        "2022-12-01"
                    ) ~
                    "2020-2022",
                
                TRUE ~
                    "2023-2026"
            )
    ) |>
    
    dplyr::group_by(
        subperiodo
    ) |>
    
    dplyr::summarise(
        
        fecha_inicial =
            min(
                fecha
            ),
        
        fecha_final =
            max(
                fecha
            ),
        
        observaciones =
            dplyr::n(),
        
        c_promedio =
            mean(
                c_t,
                na.rm = TRUE
            ),
        
        c_mediana =
            median(
                c_t,
                na.rm = TRUE
            ),
        
        c_minimo =
            min(
                c_t,
                na.rm = TRUE
            ),
        
        c_maximo =
            max(
                c_t,
                na.rm = TRUE
            ),
        
        c0_promedio =
            mean(
                c0_t,
                na.rm = TRUE
            ),
        
        d_promedio =
            mean(
                d_t,
                na.rm = TRUE
            ),
        
        .groups =
            "drop"
    )


readr::write_csv(
    resumen_subperiodos,
    file.path(
        carpeta_salida,
        "08_resumen_subperiodos_c_t.csv"
    )
)


# ==============================================================================
# 16. DIAGNÓSTICOS RESIDUALES
# ==============================================================================


extraer_residuos_estandarizados <- function(
        kfs,
        n_esperado
) {
    
    
    r <- tryCatch(
        
        as.matrix(
            
            stats::rstandard(
                kfs,
                type =
                    "recursive"
            )
        ),
        
        error =
            function(e) {
                
                NULL
            }
    )
    
    
    if (
        is.null(r) ||
        nrow(r) !=
        n_esperado ||
        ncol(r) !=
        2
    ) {
        
        stop(
            "No fue posible extraer correctamente los residuos ",
            "estandarizados recursivos del TVP-VAR en diferencias."
        )
    }
    
    
    colnames(r) <- c(
        "d_inflacion",
        "d_expectativa"
    )
    
    
    r
}


residuos_M3 <- extraer_residuos_estandarizados(
    kfs_M3,
    n_obs
)


tabla_residuos <- tibble::tibble(
    
    fecha =
        datos_modelo$fecha,
    
    residuo_d_inflacion =
        residuos_M3[
            ,
            "d_inflacion"
        ],
    
    residuo_d_expectativa =
        residuos_M3[
            ,
            "d_expectativa"
        ]
)


readr::write_csv(
    tabla_residuos,
    file.path(
        carpeta_salida,
        "11_residuos_estandarizados.csv"
    )
)


# ------------------------------------------------------------------------------
# Ljung-Box y Jarque-Bera por ecuación
# ------------------------------------------------------------------------------

crear_diagnosticos_ecuacion <- function(
        residuos,
        nombre_ecuacion
) {
    
    
    r <- residuos[
        is.finite(
            residuos
        )
    ]
    
    
    filas_LB <- lapply(
        
        rezagos_ljung_box,
        
        function(L) {
            
            
            if (
                length(r) <=
                L + 5
            ) {
                
                return(
                    NULL
                )
            }
            
            
            prueba <- stats::Box.test(
                
                r,
                
                lag =
                    L,
                
                type =
                    "Ljung-Box",
                
                fitdf =
                    0
            )
            
            
            tibble::tibble(
                
                ecuacion =
                    nombre_ecuacion,
                
                diagnostico =
                    "Ljung-Box",
                
                rezago =
                    L,
                
                estadistico =
                    unname(
                        prueba$statistic
                    ),
                
                p_valor =
                    prueba$p.value
            )
        }
    )
    
    
    JB <-
        jarque_bera_manual(
            r
        )
    
    
    fila_JB <- tibble::tibble(
        
        ecuacion =
            nombre_ecuacion,
        
        diagnostico =
            "Jarque-Bera",
        
        rezago =
            NA_real_,
        
        estadistico =
            JB$estadistico,
        
        p_valor =
            JB$p_valor
    )
    
    
    dplyr::bind_rows(
        filas_LB,
        fila_JB
    )
}


tabla_diagnosticos <- dplyr::bind_rows(
    
    crear_diagnosticos_ecuacion(
        
        residuos_M3[
            ,
            "d_inflacion"
        ],
        
        "Δ Inflación"
    ),
    
    crear_diagnosticos_ecuacion(
        
        residuos_M3[
            ,
            "d_expectativa"
        ],
        
        "Δ Expectativa 24m"
    )
)


# ------------------------------------------------------------------------------
# Correlación contemporánea entre innovaciones
# ------------------------------------------------------------------------------

rho_residuos <- correlacion_segura(
    
    residuos_M3[
        ,
        "d_inflacion"
    ],
    
    residuos_M3[
        ,
        "d_expectativa"
    ]
)


n_rho <- sum(
    
    is.finite(
        residuos_M3[
            ,
            "d_inflacion"
        ]
    ) &
        
        is.finite(
            residuos_M3[
                ,
                "d_expectativa"
            ]
        )
)


fisher_z <-
    
    if (
        is.finite(
            rho_residuos
        ) &&
        abs(
            rho_residuos
        ) < 1 &&
        n_rho > 3
    ) {
        
        atanh(
            rho_residuos
        ) *
            sqrt(
                n_rho - 3
            )
        
    } else {
        
        NA_real_
    }


p_fisher <-
    
    if (
        is.finite(
            fisher_z
        )
    ) {
        
        2 *
            stats::pnorm(
                abs(
                    fisher_z
                ),
                lower.tail = FALSE
            )
        
    } else {
        
        NA_real_
    }


tabla_diagnosticos <- dplyr::bind_rows(
    
    tabla_diagnosticos,
    
    tibble::tibble(
        
        ecuacion =
            "Sistema",
        
        diagnostico =
            "Correlación contemporánea Fisher",
        
        rezago =
            NA_real_,
        
        estadistico =
            fisher_z,
        
        p_valor =
            p_fisher
    )
) |>
    
    dplyr::mutate(
        
        conclusion_5pct =
            dplyr::case_when(
                
                !is.finite(
                    p_valor
                ) ~
                    "No disponible",
                
                p_valor <
                    0.05 ~
                    "Se rechaza H0",
                
                TRUE ~
                    "No se rechaza H0"
            )
    )


readr::write_csv(
    tabla_diagnosticos,
    file.path(
        carpeta_salida,
        "09_diagnosticos_residuales.csv"
    )
)


# ==============================================================================
# 17. FIGURA PRINCIPAL DE ROBUSTEZ: c_t
# ==============================================================================

p_c <- ggplot2::ggplot(
    
    trayectoria_M3,
    
    ggplot2::aes(
        x =
            fecha,
        y =
            c_t
    )
) +
    
    ggplot2::geom_ribbon(
        
        ggplot2::aes(
            ymin =
                c_li_95,
            ymax =
                c_ls_95
        ),
        
        alpha =
            0.16,
        
        linewidth =
            0
    ) +
    
    ggplot2::geom_hline(
        
        yintercept =
            0,
        
        linetype =
            "dashed",
        
        linewidth =
            0.4
    ) +
    
    ggplot2::geom_line(
        linewidth =
            0.75
    ) +
    
    ggplot2::scale_x_date(
        
        date_breaks =
            "2 years",
        
        date_labels =
            "%Y"
    ) +
    
    ggplot2::labs(
        
        x =
            NULL,
        
        y =
            expression(
                hat(c)[t]
            )
    ) +
    
    ggplot2::theme_minimal(
        base_size =
            10.5
    ) +
    
    ggplot2::theme(
        
        panel.grid.minor =
            ggplot2::element_blank(),
        
        axis.text =
            ggplot2::element_text(
                colour =
                    "grey25"
            )
    )


ggplot2::ggsave(
    
    filename =
        file.path(
            carpeta_salida,
            "figura_01_c_t_diferencias.png"
        ),
    
    plot =
        p_c,
    
    width =
        8.5,
    
    height =
        4.3,
    
    dpi =
        300,
    
    bg =
        "white"
)


print(
    p_c
)


# ==============================================================================
# 18. FIGURA DE c0_t
# ==============================================================================

p_c0 <- ggplot2::ggplot(
    
    trayectoria_M3,
    
    ggplot2::aes(
        x =
            fecha,
        y =
            c0_t
    )
) +
    
    ggplot2::geom_ribbon(
        
        ggplot2::aes(
            ymin =
                c0_li_95,
            ymax =
                c0_ls_95
        ),
        
        alpha =
            0.16,
        
        linewidth =
            0
    ) +
    
    ggplot2::geom_hline(
        
        yintercept =
            0,
        
        linetype =
            "dashed",
        
        linewidth =
            0.4
    ) +
    
    ggplot2::geom_line(
        linewidth =
            0.75
    ) +
    
    ggplot2::scale_x_date(
        
        date_breaks =
            "2 years",
        
        date_labels =
            "%Y"
    ) +
    
    ggplot2::labs(
        
        x =
            NULL,
        
        y =
            expression(
                hat(c)[0 * t]
            )
    ) +
    
    ggplot2::theme_minimal(
        base_size =
            10.5
    ) +
    
    ggplot2::theme(
        panel.grid.minor =
            ggplot2::element_blank()
    )


ggplot2::ggsave(
    
    filename =
        file.path(
            carpeta_salida,
            "figura_02_c0_t_diferencias.png"
        ),
    
    plot =
        p_c0,
    
    width =
        8.5,
    
    height =
        4.3,
    
    dpi =
        300,
    
    bg =
        "white"
)


# ==============================================================================
# 19. GUARDAR OBJETOS
# ==============================================================================

saveRDS(
    
    ajuste_M0$modelo,
    
    file.path(
        carpeta_salida,
        "modelo_M0_D.rds"
    )
)


saveRDS(
    
    ajuste_M3$modelo,
    
    file.path(
        carpeta_salida,
        "modelo_M3_D.rds"
    )
)


saveRDS(
    
    kfs_M3,
    
    file.path(
        carpeta_salida,
        "kfs_M3_D.rds"
    )
)


# ==============================================================================
# 20. RESUMEN EJECUTIVO
# ==============================================================================


Q_c0_hat <- tabla_hiperparametros |>
    
    dplyr::filter(
        modelo ==
            "M3_D"
    ) |>
    
    dplyr::pull(
        Q_c0
    )


Q_c_hat <- tabla_hiperparametros |>
    
    dplyr::filter(
        modelo ==
            "M3_D"
    ) |>
    
    dplyr::pull(
        Q_c
    )


p_LB_pi_12 <- tabla_diagnosticos |>
    
    dplyr::filter(
        ecuacion ==
            "Δ Inflación",
        diagnostico ==
            "Ljung-Box",
        rezago ==
            12
    ) |>
    
    dplyr::pull(
        p_valor
    )


p_LB_e_12 <- tabla_diagnosticos |>
    
    dplyr::filter(
        ecuacion ==
            "Δ Expectativa 24m",
        diagnostico ==
            "Ljung-Box",
        rezago ==
            12
    ) |>
    
    dplyr::pull(
        p_valor
    )


c_inicial <-
    resumen_c$c_inicial


c_final <-
    resumen_c$c_final


direccion_c <-
    
    if (
        is.finite(
            c_inicial
        ) &&
        is.finite(
            c_final
        )
    ) {
        
        if (
            c_final <
            c_inicial
        ) {
            
            "La sensibilidad c_t disminuye entre el inicio y el final de la muestra."
            
        } else if (
            c_final >
            c_inicial
        ) {
            
            "La sensibilidad c_t aumenta entre el inicio y el final de la muestra."
            
        } else {
            
            "La sensibilidad c_t no presenta cambio entre los extremos de la muestra."
        }
        
    } else {
        
        "No fue posible determinar automáticamente la dirección de c_t."
    }


lineas_resumen <- c(
    
    "ROBUSTEZ: TVP-VAR RESTRINGIDO EN PRIMERAS DIFERENCIAS",
    "=======================================================",
    "",
    
    paste0(
        "Muestra efectiva: ",
        format(
            min(
                datos_modelo$fecha
            ),
            "%Y-%m"
        ),
        " a ",
        format(
            max(
                datos_modelo$fecha
            ),
            "%Y-%m"
        )
    ),
    
    paste0(
        "Observaciones: ",
        n_obs
    ),
    
    "",
    
    "MODELO",
    "------",
    
    paste0(
        "LogLik M0_D: ",
        round(
            ll_M0,
            4
        )
    ),
    
    paste0(
        "LogLik M3_D: ",
        round(
            ll_M3,
            4
        )
    ),
    
    paste0(
        "AIC M0_D: ",
        round(
            AIC_M0,
            4
        )
    ),
    
    paste0(
        "AIC M3_D: ",
        round(
            AIC_M3,
            4
        )
    ),
    
    paste0(
        "BIC M0_D: ",
        round(
            BIC_M0,
            4
        )
    ),
    
    paste0(
        "BIC M3_D: ",
        round(
            BIC_M3,
            4
        )
    ),
    
    "",
    
    "VARIANZAS DE TRANSICIÓN",
    "-----------------------",
    
    paste0(
        "Q_c0 = ",
        format(
            Q_c0_hat,
            scientific =
                TRUE
        )
    ),
    
    paste0(
        "Q_c  = ",
        format(
            Q_c_hat,
            scientific =
                TRUE
        )
    ),
    
    "",
    
    "TRAYECTORIA DE LA SENSIBILIDAD c_t",
    "-----------------------------------",
    
    paste0(
        "c_t inicial = ",
        round(
            c_inicial,
            4
        )
    ),
    
    paste0(
        "c_t final   = ",
        round(
            c_final,
            4
        )
    ),
    
    paste0(
        "c_t promedio = ",
        round(
            resumen_c$c_promedio,
            4
        )
    ),
    
    direccion_c,
    
    "",
    
    "DIAGNÓSTICOS",
    "------------",
    
    paste0(
        "Ljung-Box 12 rezagos, Δ inflación: p = ",
        format.pval(
            p_LB_pi_12,
            digits =
                4,
            eps =
                0.0001
        )
    ),
    
    paste0(
        "Ljung-Box 12 rezagos, Δ expectativa: p = ",
        format.pval(
            p_LB_e_12,
            digits =
                4,
            eps =
                0.0001
        )
    ),
    
    paste0(
        "Porcentaje de meses localmente estables: ",
        round(
            resumen_c$porcentaje_meses_estables,
            2
        ),
        "%"
    ),
    
    paste0(
        "Radio espectral máximo: ",
        round(
            resumen_c$radio_espectral_maximo,
            4
        )
    ),
    
    "",
    
    "INTERPRETACIÓN",
    "--------------",
    
    paste0(
        "Este ejercicio debe interpretarse exclusivamente como una ",
        "comprobación de robustez del patrón de sensibilidad."
    ),
    
    paste0(
        "No deben construirse lambda_t ni pi_estrella_t a partir ",
        "de esta especificación en primeras diferencias."
    )
)


writeLines(
    
    lineas_resumen,
    
    con =
        file.path(
            carpeta_salida,
            "00_resumen_ejecutivo.txt"
        )
)


# ==============================================================================
# 21. RESULTADOS EN CONSOLA
# ==============================================================================

cat(
    "\n============================================================\n"
)

cat(
    "TVP-VAR EN PRIMERAS DIFERENCIAS: ESTIMACIÓN FINALIZADA\n"
)

cat(
    "============================================================\n\n"
)


cat(
    "Comparación M0_D vs M3_D:\n"
)

print(
    tabla_comparacion
)


cat(
    "\nHiperparámetros:\n"
)

print(
    tabla_hiperparametros
)


cat(
    "\nResumen de c_t:\n"
)

print(
    resumen_c
)


cat(
    "\nResumen por subperíodos:\n"
)

print(
    resumen_subperiodos
)


cat(
    "\nDiagnósticos residuales:\n"
)

print(
    tabla_diagnosticos
)


cat(
    "\nResultados guardados en:\n",
    carpeta_salida,
    "\n"
)


cat(
    "\n============================================================\n"
)