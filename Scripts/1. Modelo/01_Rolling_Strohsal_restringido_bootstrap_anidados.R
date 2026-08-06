# ============================================================
# 01. ROLLING REGRESSIONS — STROHSAL RESTRINGIDO
# Indicador principal con restricciones, bootstrap por bloques
# y modelos anidados de diagnóstico
# ============================================================

# MODELO COMPLETO DE STROHSAL
#
# e24_t = (1 - theta1 - theta2) * meta_t
#         + theta1 * inflacion_(t-1)
#         + theta2 * e12_(t-1)
#         + u_t
#
# Forma estimable:
#
# y_t = theta1 * x1_t + theta2 * x2_t + u_t
#
# donde:
#
# y_t  = Delta e24_t - (meta_t - e24_(t-1))
# x1_t = inflacion_(t-1) - meta_t
# x2_t = e12_(t-1) - meta_t
#
# Indicador de anclaje:
#
# A = 1 - theta1 - theta2
#
# Restricciones económicas:
#
# theta1 >= 0
# theta2 >= 0
# theta1 + theta2 <= 1
#
# Por tanto:
#
# 0 <= A <= 1
#
# ESTRATEGIA
#
# 1. Modelo restringido M3 como medición principal.
# 2. Bootstrap circular por bloques de residuos para la inferencia.
# 3. Modelos anidados irrestrictos:
#       M1: solo inflación.
#       M2: solo expectativas a 12 meses.
#       M3: inflación + expectativas a 12 meses.
# 4. Comparación entre M3 irrestricto y M3 restringido.
#
# IMPORTANTE:
# - Las restricciones garantizan una interpretación como ponderaciones.
# - El bootstrap cuantifica la incertidumbre bajo las restricciones.
# - Los modelos anidados ayudan a diagnosticar información superpuesta.
# - Este procedimiento no corrige por sí mismo una posible endogeneidad.
# ============================================================

# ------------------------------------------------------------
# Paquetes
# ------------------------------------------------------------
library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(lmtest)
library(sandwich)
library(quadprog)

# Instalar quadprog una sola vez, si fuera necesario:
# install.packages("quadprog")

# ------------------------------------------------------------
# 0) Limpieza de entorno
# ------------------------------------------------------------
rm(list = ls())
cat("\014")
if (!is.null(dev.list())) dev.off()

# ------------------------------------------------------------
# 1) Parámetros generales
# ------------------------------------------------------------

window_size  <- 36       # Ventana principal
nw_lag       <- 3        # Rezagos para Newey-West en modelos anidados
block_length <- 4        # Longitud del bloque bootstrap, en meses
B_boot       <- 1000     # Réplicas bootstrap finales
seed_boot    <- 20260801
tol_frontera <- 1e-6
ridge_eps    <- 1e-10    # Estabilización numérica mínima

# Para depurar el script con mayor rapidez:
# B_boot <- 199

# Robustez de la ventana:
# window_size <- 30
# window_size <- 48
# window_size <- 60

set.seed(seed_boot)

# ------------------------------------------------------------
# 2) Cargar base y construir variables
# ------------------------------------------------------------

df_model <- read_csv(
    "bd_sample_modelo.csv",
    show_col_types = FALSE
) |>
    mutate(
        fecha = ymd(fecha)
    ) |>
    arrange(fecha) |>
    mutate(
        # Meta central vigente
        meta_inflacion = case_when(
            fecha < ymd("2012-01-01") ~ 5.0,
            fecha < ymd("2013-01-01") ~ 4.5,
            TRUE                       ~ 4.0
        ),

        # Rezagos construidos directamente desde las series
        exp12_lag1 = lag(exp_inf_12m, 1),
        exp24_lag1_strohsal = lag(exp_inf_24m, 1),
        infl_lag1_strohsal = lag(infl_gt, 1),

        # Cambio mensual de la expectativa a 24 meses
        d_exp24 =
            exp_inf_24m - exp24_lag1_strohsal,

        # Variable dependiente transformada
        y_strohsal =
            d_exp24 -
            (
                meta_inflacion -
                exp24_lag1_strohsal
            ),

        # Fuente de desanclaje 1:
        # inflación rezagada respecto de la meta
        x_inflacion =
            infl_lag1_strohsal -
            meta_inflacion,

        # Fuente de desanclaje 2:
        # expectativas a 12 meses rezagadas respecto de la meta
        x_exp12 =
            exp12_lag1 -
            meta_inflacion
    ) |>
    filter(
        complete.cases(
            y_strohsal,
            x_inflacion,
            x_exp12
        )
    )

cat(
    "\nPrimera fecha efectiva: ",
    as.character(min(df_model$fecha)),
    "\n"
)

cat(
    "Última fecha efectiva: ",
    as.character(max(df_model$fecha)),
    "\n"
)

cat(
    "Observaciones efectivas: ",
    nrow(df_model),
    "\n"
)

if (nrow(df_model) < window_size) {
    stop(
        "La muestra efectiva tiene menos observaciones que la ventana móvil."
    )
}

if (block_length > window_size) {
    stop(
        "La longitud del bloque no puede superar el tamaño de la ventana."
    )
}

# ------------------------------------------------------------
# 3) Funciones auxiliares
# ------------------------------------------------------------

# 3.1) Extraer un coeficiente de coeftest como escalares
extraer_coeficiente <- function(coeftest_obj, nombre) {

    if (!nombre %in% rownames(coeftest_obj)) {
        return(
            c(
                estimate = NA_real_,
                std_error = NA_real_,
                statistic = NA_real_,
                p_value = NA_real_
            )
        )
    }

    c(
        estimate  = unname(coeftest_obj[nombre, "Estimate"]),
        std_error = unname(coeftest_obj[nombre, "Std. Error"]),
        statistic = unname(coeftest_obj[nombre, "t value"]),
        p_value   = unname(coeftest_obj[nombre, "Pr(>|t|)"])
    )
}

# 3.2) Estimar un modelo irrestricto con errores Newey-West HAC
estimar_ols_hac <- function(formula, data_window, nw_lag) {

    modelo <- lm(
        formula,
        data = data_window
    )

    vcov_hac <- sandwich::NeweyWest(
        modelo,
        lag = nw_lag,
        prewhite = FALSE,
        adjust = TRUE
    )

    prueba_hac <- lmtest::coeftest(
        modelo,
        vcov. = vcov_hac
    )

    residuos <- residuals(modelo)

    list(
        modelo = modelo,
        vcov_hac = vcov_hac,
        prueba_hac = prueba_hac,
        fitted = fitted(modelo),
        residuos = residuos,
        sse = sum(residuos^2),
        rmse = sqrt(mean(residuos^2)),
        n = nobs(modelo)
    )
}

# 3.3) Estimar el modelo M3 con restricciones mediante
#      mínimos cuadrados cuadráticos
estimar_restringido <- function(
    y,
    X,
    ridge_eps = 1e-10,
    tol_frontera = 1e-6
) {

    y <- as.numeric(y)
    X <- as.matrix(X)
    storage.mode(X) <- "double"

    if (ncol(X) != 2) {
        stop(
            "El modelo restringido requiere exactamente dos regresores."
        )
    }

    # solve.QP minimiza:
    # 1/2 b' D b - d' b
    #
    # Con:
    # D = 2 X'X
    # d = 2 X'y
    #
    # se minimiza la suma de residuos al cuadrado.
    Dmat <-
        2 * crossprod(X) +
        diag(ridge_eps, ncol(X))

    dvec <-
        as.numeric(
            2 * crossprod(X, y)
        )

    # Restricciones expresadas como:
    # t(Amat) %*% b >= bvec
    #
    # theta1 >= 0
    # theta2 >= 0
    # -theta1 - theta2 >= -1
    Amat <- cbind(
        c(1, 0),
        c(0, 1),
        c(-1, -1)
    )

    bvec <- c(
        0,
        0,
        -1
    )

    ajuste_qp <- quadprog::solve.QP(
        Dmat = Dmat,
        dvec = dvec,
        Amat = Amat,
        bvec = bvec,
        meq = 0
    )

    coeficientes <- as.numeric(
        ajuste_qp$solution
    )

    names(coeficientes) <- c(
        "theta1",
        "theta2"
    )

    # Limpieza de errores numéricos muy pequeños
    coeficientes[
        abs(coeficientes) < tol_frontera
    ] <- 0

    suma_theta <- sum(coeficientes)

    if (
        suma_theta > 1 &&
        suma_theta <= 1 + tol_frontera
    ) {
        coeficientes <-
            coeficientes / suma_theta
    }

    theta1 <- coeficientes["theta1"]
    theta2 <- coeficientes["theta2"]
    anclaje <- 1 - theta1 - theta2

    if (
        anclaje < 0 &&
        anclaje >= -tol_frontera
    ) {
        anclaje <- 0
    }

    ajustados <- as.numeric(
        X %*% coeficientes
    )

    residuos <- y - ajustados

    list(
        theta1 = unname(theta1),
        theta2 = unname(theta2),
        anclaje = unname(anclaje),
        fitted = ajustados,
        residuos = residuos,
        sse = sum(residuos^2),
        rmse = sqrt(mean(residuos^2)),
        theta1_frontera =
            unname(theta1 <= tol_frontera),
        theta2_frontera =
            unname(theta2 <= tol_frontera),
        anclaje_frontera =
            unname(anclaje <= tol_frontera)
    )
}

# 3.4) Remuestreo circular por bloques
remuestrear_bloques_circulares <- function(
    residuos,
    block_length
) {

    residuos <- as.numeric(residuos)
    n <- length(residuos)

    # Los residuos se centran antes del remuestreo.
    residuos <- residuos - mean(residuos)

    n_bloques <- ceiling(
        n / block_length
    )

    inicios <- sample.int(
        n,
        size = n_bloques,
        replace = TRUE
    )

    residuos_boot <- unlist(
        lapply(
            inicios,
            function(inicio) {

                indices <-
                    (
                        inicio - 1 +
                        0:(block_length - 1)
                    ) %% n + 1

                residuos[indices]
            }
        ),
        use.names = FALSE
    )

    residuos_boot[
        seq_len(n)
    ]
}

# 3.5) Bootstrap por bloques para el modelo restringido
bootstrap_restringido <- function(
    y,
    X,
    ajuste_restringido,
    B_boot,
    block_length,
    ridge_eps,
    tol_frontera
) {

    n <- length(y)

    resultados_boot <- matrix(
        NA_real_,
        nrow = B_boot,
        ncol = 3,
        dimnames = list(
            NULL,
            c(
                "theta1",
                "theta2",
                "anclaje"
            )
        )
    )

    for (b in seq_len(B_boot)) {

        residuos_b <-
            remuestrear_bloques_circulares(
                ajuste_restringido$residuos,
                block_length
            )

        y_b <-
            ajuste_restringido$fitted +
            residuos_b

        ajuste_b <- tryCatch(
            estimar_restringido(
                y = y_b,
                X = X,
                ridge_eps = ridge_eps,
                tol_frontera = tol_frontera
            ),
            error = function(e) NULL
        )

        if (!is.null(ajuste_b)) {

            resultados_boot[b, ] <- c(
                ajuste_b$theta1,
                ajuste_b$theta2,
                ajuste_b$anclaje
            )
        }
    }

    resultados_boot <-
        resultados_boot[
            complete.cases(resultados_boot),
            ,
            drop = FALSE
        ]

    n_validas <- nrow(
        resultados_boot
    )

    if (n_validas == 0) {
        return(
            tibble(
                n_boot_validas = 0,
                theta1_se_boot = NA_real_,
                theta1_ic_inf_boot = NA_real_,
                theta1_ic_sup_boot = NA_real_,
                theta2_se_boot = NA_real_,
                theta2_ic_inf_boot = NA_real_,
                theta2_ic_sup_boot = NA_real_,
                anclaje_se_boot = NA_real_,
                anclaje_ic_inf_boot = NA_real_,
                anclaje_ic_sup_boot = NA_real_,
                anclaje_sesgo_boot = NA_real_,
                freq_theta1_frontera_boot = NA_real_,
                freq_theta2_frontera_boot = NA_real_,
                freq_anclaje_frontera_boot = NA_real_
            )
        )
    }

    cuantiles <- apply(
        resultados_boot,
        2,
        quantile,
        probs = c(0.025, 0.975),
        na.rm = TRUE,
        names = FALSE
    )

    tibble(
        n_boot_validas = n_validas,

        theta1_se_boot =
            sd(resultados_boot[, "theta1"]),
        theta1_ic_inf_boot =
            cuantiles[1, "theta1"],
        theta1_ic_sup_boot =
            cuantiles[2, "theta1"],

        theta2_se_boot =
            sd(resultados_boot[, "theta2"]),
        theta2_ic_inf_boot =
            cuantiles[1, "theta2"],
        theta2_ic_sup_boot =
            cuantiles[2, "theta2"],

        anclaje_se_boot =
            sd(resultados_boot[, "anclaje"]),
        anclaje_ic_inf_boot =
            cuantiles[1, "anclaje"],
        anclaje_ic_sup_boot =
            cuantiles[2, "anclaje"],

        anclaje_sesgo_boot =
            mean(
                resultados_boot[, "anclaje"]
            ) -
            ajuste_restringido$anclaje,

        freq_theta1_frontera_boot =
            mean(
                resultados_boot[, "theta1"] <=
                tol_frontera
            ),

        freq_theta2_frontera_boot =
            mean(
                resultados_boot[, "theta2"] <=
                tol_frontera
            ),

        freq_anclaje_frontera_boot =
            mean(
                resultados_boot[, "anclaje"] <=
                tol_frontera
            )
    )
}

# ------------------------------------------------------------
# 4) Función integral para cada ventana
# ------------------------------------------------------------

estimar_ventana <- function(
    data_window,
    nw_lag,
    B_boot,
    block_length,
    ridge_eps,
    tol_frontera
) {

    # --------------------------------------------------------
    # 4.1) Modelos anidados irrestrictos
    # --------------------------------------------------------

    # M1: solo inflación
    m1 <- estimar_ols_hac(
        y_strohsal ~ 0 + x_inflacion,
        data_window,
        nw_lag
    )

    # M2: solo expectativas a 12 meses
    m2 <- estimar_ols_hac(
        y_strohsal ~ 0 + x_exp12,
        data_window,
        nw_lag
    )

    # M3: modelo completo irrestricto
    m3_ols <- estimar_ols_hac(
        y_strohsal ~
            0 +
            x_inflacion +
            x_exp12,
        data_window,
        nw_lag
    )

    coef_theta1_m1 <-
        extraer_coeficiente(
            m1$prueba_hac,
            "x_inflacion"
        )

    coef_theta2_m2 <-
        extraer_coeficiente(
            m2$prueba_hac,
            "x_exp12"
        )

    coef_theta1_m3 <-
        extraer_coeficiente(
            m3_ols$prueba_hac,
            "x_inflacion"
        )

    coef_theta2_m3 <-
        extraer_coeficiente(
            m3_ols$prueba_hac,
            "x_exp12"
        )

    anclaje_m3_ols <-
        1 -
        coef_theta1_m3["estimate"] -
        coef_theta2_m3["estimate"]

    # Intervalo HAC de A en el modelo irrestricto
    g_anclaje <- c(
        -1,
        -1
    )

    var_anclaje_m3_ols <-
        as.numeric(
            t(g_anclaje) %*%
                m3_ols$vcov_hac[
                    c(
                        "x_inflacion",
                        "x_exp12"
                    ),
                    c(
                        "x_inflacion",
                        "x_exp12"
                    )
                ] %*%
                g_anclaje
        )

    se_anclaje_m3_ols <-
        sqrt(
            max(
                var_anclaje_m3_ols,
                0
            )
        )

    # --------------------------------------------------------
    # 4.2) Modelo completo restringido
    # --------------------------------------------------------

    y <- data_window$y_strohsal

    X <- as.matrix(
        data_window |>
            select(
                x_inflacion,
                x_exp12
            )
    )

    m3_restringido <-
        estimar_restringido(
            y = y,
            X = X,
            ridge_eps = ridge_eps,
            tol_frontera = tol_frontera
        )

    # --------------------------------------------------------
    # 4.3) Bootstrap por bloques del modelo restringido
    # --------------------------------------------------------

    boot <-
        bootstrap_restringido(
            y = y,
            X = X,
            ajuste_restringido =
                m3_restringido,
            B_boot = B_boot,
            block_length =
                block_length,
            ridge_eps = ridge_eps,
            tol_frontera =
                tol_frontera
        )

    # --------------------------------------------------------
    # 4.4) Diagnósticos de información superpuesta
    # --------------------------------------------------------

    correlacion_x1_x2 <-
        cor(
            data_window$x_inflacion,
            data_window$x_exp12,
            use = "complete.obs"
        )

    vif_aproximado <-
        ifelse(
            abs(correlacion_x1_x2) < 1,
            1 /
                (
                    1 -
                    correlacion_x1_x2^2
                ),
            Inf
        )

    numero_condicion_X <-
        kappa(X)

    # Coste de imponer las restricciones
    penalizacion_sse <-
        m3_restringido$sse -
        m3_ols$sse

    penalizacion_sse_pct <-
        ifelse(
            m3_ols$sse > 0,
            100 *
                penalizacion_sse /
                m3_ols$sse,
            NA_real_
        )

    # Aporte marginal en términos de RMSE:
    # positivo = el modelo completo mejora frente al modelo reducido
    mejora_rmse_al_agregar_inflacion <-
        m2$rmse -
        m3_ols$rmse

    mejora_rmse_al_agregar_exp12 <-
        m1$rmse -
        m3_ols$rmse

    tibble(
        # ----------------------------------------------------
        # M1: solo inflación
        # ----------------------------------------------------
        theta1_m1 =
            coef_theta1_m1["estimate"],
        se_theta1_m1_hac =
            coef_theta1_m1["std_error"],
        p_theta1_m1_hac =
            coef_theta1_m1["p_value"],
        rmse_m1 =
            m1$rmse,
        sse_m1 =
            m1$sse,

        # ----------------------------------------------------
        # M2: solo expectativas a 12 meses
        # ----------------------------------------------------
        theta2_m2 =
            coef_theta2_m2["estimate"],
        se_theta2_m2_hac =
            coef_theta2_m2["std_error"],
        p_theta2_m2_hac =
            coef_theta2_m2["p_value"],
        rmse_m2 =
            m2$rmse,
        sse_m2 =
            m2$sse,

        # ----------------------------------------------------
        # M3 irrestricto
        # ----------------------------------------------------
        theta1_m3_ols =
            coef_theta1_m3["estimate"],
        se_theta1_m3_hac =
            coef_theta1_m3["std_error"],
        p_theta1_m3_hac =
            coef_theta1_m3["p_value"],

        theta2_m3_ols =
            coef_theta2_m3["estimate"],
        se_theta2_m3_hac =
            coef_theta2_m3["std_error"],
        p_theta2_m3_hac =
            coef_theta2_m3["p_value"],

        anclaje_m3_ols =
            unname(anclaje_m3_ols),
        se_anclaje_m3_hac =
            se_anclaje_m3_ols,
        anclaje_m3_ols_ic_inf =
            unname(anclaje_m3_ols) -
            1.96 * se_anclaje_m3_ols,
        anclaje_m3_ols_ic_sup =
            unname(anclaje_m3_ols) +
            1.96 * se_anclaje_m3_ols,

        rmse_m3_ols =
            m3_ols$rmse,
        sse_m3_ols =
            m3_ols$sse,

        # ----------------------------------------------------
        # M3 restringido: indicador principal
        # ----------------------------------------------------
        theta1_restringido =
            m3_restringido$theta1,
        theta2_restringido =
            m3_restringido$theta2,
        anclaje_restringido =
            m3_restringido$anclaje,

        theta1_frontera =
            m3_restringido$theta1_frontera,
        theta2_frontera =
            m3_restringido$theta2_frontera,
        anclaje_frontera =
            m3_restringido$anclaje_frontera,

        rmse_m3_restringido =
            m3_restringido$rmse,
        sse_m3_restringido =
            m3_restringido$sse,

        # ----------------------------------------------------
        # Diagnósticos
        # ----------------------------------------------------
        correlacion_x1_x2 =
            correlacion_x1_x2,
        vif_aproximado =
            vif_aproximado,
        numero_condicion_X =
            numero_condicion_X,

        mejora_rmse_al_agregar_inflacion =
            mejora_rmse_al_agregar_inflacion,
        mejora_rmse_al_agregar_exp12 =
            mejora_rmse_al_agregar_exp12,

        penalizacion_sse =
            penalizacion_sse,
        penalizacion_sse_pct =
            penalizacion_sse_pct,

        n = nrow(data_window)
    ) |>
        bind_cols(boot)
}

# ------------------------------------------------------------
# 5) Ciclo móvil
# ------------------------------------------------------------

indices_fin <-
    window_size:nrow(df_model)

rolling_list <-
    vector(
        "list",
        length(indices_fin)
    )

for (j in seq_along(indices_fin)) {

    i <- indices_fin[j]

    data_window <- df_model[
        (i - window_size + 1):i,
    ]

    if (
        j == 1 ||
        j %% 10 == 0 ||
        j == length(indices_fin)
    ) {
        message(
            "Estimando ventana ",
            j,
            " de ",
            length(indices_fin),
            " | Fecha final: ",
            max(data_window$fecha)
        )
    }

    resultado_j <- estimar_ventana(
        data_window =
            data_window,
        nw_lag =
            nw_lag,
        B_boot =
            B_boot,
        block_length =
            block_length,
        ridge_eps =
            ridge_eps,
        tol_frontera =
            tol_frontera
    )

    rolling_list[[j]] <-
        tibble(
            fecha_inicio =
                min(data_window$fecha),
            fecha_fin =
                max(data_window$fecha)
        ) |>
        bind_cols(
            resultado_j
        )
}

rolling_results <-
    bind_rows(
        rolling_list
    )

# ------------------------------------------------------------
# 6) Guardar resultados
# ------------------------------------------------------------

if (!dir.exists("output")) {
    dir.create("output")
}

nombre_base <- paste0(
    "01_rolling_strohsal_restringido_",
    window_size,
    "m_B",
    B_boot,
    "_L",
    block_length
)

write_csv(
    rolling_results,
    paste0(
        "output/",
        nombre_base,
        ".csv"
    )
)

# Resumen de las restricciones activas
resumen_fronteras <- tibble(
    indicador = c(
        "Theta 1 en cero",
        "Theta 2 en cero",
        "Anclaje en cero"
    ),
    numero_ventanas = c(
        sum(
            rolling_results$theta1_frontera,
            na.rm = TRUE
        ),
        sum(
            rolling_results$theta2_frontera,
            na.rm = TRUE
        ),
        sum(
            rolling_results$anclaje_frontera,
            na.rm = TRUE
        )
    ),
    porcentaje_ventanas = 100 * numero_ventanas /
        nrow(rolling_results)
)

write_csv(
    resumen_fronteras,
    paste0(
        "output/",
        nombre_base,
        "_resumen_fronteras.csv"
    )
)

# ------------------------------------------------------------
# 7) Gráfico principal: indicador restringido de anclaje
# ------------------------------------------------------------

p_anclaje <- ggplot(
    rolling_results,
    aes(
        x = fecha_fin,
        y = anclaje_restringido
    )
) +
    geom_ribbon(
        aes(
            ymin =
                anclaje_ic_inf_boot,
            ymax =
                anclaje_ic_sup_boot
        ),
        alpha = 0.2
    ) +
    geom_line(
        linewidth = 1
    ) +
    geom_hline(
        yintercept = 1,
        linetype = "dotted"
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    coord_cartesian(
        ylim = c(0, 1)
    ) +
    labs(
        title = paste0(
            "Indicador restringido de anclaje ",
            "(ventana de ",
            window_size,
            " meses)"
        ),
        subtitle = paste0(
            "Intervalos percentiles al 95% mediante bootstrap ",
            "circular por bloques"
        ),
        x = "Fecha final de la ventana",
        y = expression(
            A[w] == 1 - theta[1*w] - theta[2*w]
        )
    ) +
    theme_minimal()

ggsave(
    paste0(
        "output/",
        nombre_base,
        "_anclaje.png"
    ),
    plot = p_anclaje,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_anclaje)

# ------------------------------------------------------------
# 8) Gráfico de los componentes restringidos
# ------------------------------------------------------------

componentes_restringidos <- bind_rows(
    rolling_results |>
        transmute(
            fecha_fin,
            componente =
                "Inflación pasada",
            estimacion =
                theta1_restringido,
            ic_inf =
                theta1_ic_inf_boot,
            ic_sup =
                theta1_ic_sup_boot
        ),

    rolling_results |>
        transmute(
            fecha_fin,
            componente =
                "Expectativas a 12 meses",
            estimacion =
                theta2_restringido,
            ic_inf =
                theta2_ic_inf_boot,
            ic_sup =
                theta2_ic_sup_boot
        )
)

p_componentes <- ggplot(
    componentes_restringidos,
    aes(
        x = fecha_fin,
        y = estimacion
    )
) +
    geom_ribbon(
        aes(
            ymin = ic_inf,
            ymax = ic_sup
        ),
        alpha = 0.2
    ) +
    geom_line(
        linewidth = 0.9
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    facet_wrap(
        ~ componente,
        ncol = 1
    ) +
    coord_cartesian(
        ylim = c(0, 1)
    ) +
    labs(
        title = paste0(
            "Fuentes restringidas de desanclaje ",
            "(ventana de ",
            window_size,
            " meses)"
        ),
        subtitle = paste0(
            "Intervalos percentiles al 95% mediante bootstrap ",
            "circular por bloques"
        ),
        x = "Fecha final de la ventana",
        y = "Ponderación estimada"
    ) +
    theme_minimal()

ggsave(
    paste0(
        "output/",
        nombre_base,
        "_componentes.png"
    ),
    plot = p_componentes,
    width = 10,
    height = 8,
    dpi = 300
)

print(p_componentes)

# ------------------------------------------------------------
# 9) Comparación: anclaje irrestricto y restringido
# ------------------------------------------------------------

comparacion_anclaje <- rolling_results |>
    select(
        fecha_fin,
        anclaje_restringido,
        anclaje_m3_ols
    ) |>
    pivot_longer(
        cols = c(
            anclaje_restringido,
            anclaje_m3_ols
        ),
        names_to = "modelo",
        values_to = "anclaje"
    ) |>
    mutate(
        modelo = recode(
            modelo,
            anclaje_restringido =
                "M3 restringido",
            anclaje_m3_ols =
                "M3 irrestricto"
        )
    )

p_comparacion_anclaje <- ggplot(
    comparacion_anclaje,
    aes(
        x = fecha_fin,
        y = anclaje,
        linetype = modelo
    )
) +
    geom_line(
        linewidth = 0.9
    ) +
    geom_hline(
        yintercept = 1,
        linetype = "dotted"
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    labs(
        title = paste0(
            "Comparación del indicador de anclaje ",
            "(ventana de ",
            window_size,
            " meses)"
        ),
        subtitle = paste0(
            "Modelo completo irrestricto frente al modelo ",
            "económicamente restringido"
        ),
        x = "Fecha final de la ventana",
        y = "Indicador de anclaje",
        linetype = "Modelo"
    ) +
    theme_minimal()

ggsave(
    paste0(
        "output/",
        nombre_base,
        "_comparacion_anclaje.png"
    ),
    plot = p_comparacion_anclaje,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_comparacion_anclaje)

# ------------------------------------------------------------
# 10) Modelos anidados: efecto de la inflación
# ------------------------------------------------------------

comparacion_theta1 <- rolling_results |>
    select(
        fecha_fin,
        theta1_m1,
        theta1_m3_ols
    ) |>
    pivot_longer(
        cols = c(
            theta1_m1,
            theta1_m3_ols
        ),
        names_to = "modelo",
        values_to = "theta1"
    ) |>
    mutate(
        modelo = recode(
            modelo,
            theta1_m1 =
                "M1: solo inflación",
            theta1_m3_ols =
                "M3: modelo completo"
        )
    )

p_anidados_theta1 <- ggplot(
    comparacion_theta1,
    aes(
        x = fecha_fin,
        y = theta1,
        linetype = modelo
    )
) +
    geom_line(
        linewidth = 0.9
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    labs(
        title = paste0(
            "Modelos anidados: coeficiente de inflación ",
            "(ventana de ",
            window_size,
            " meses)"
        ),
        subtitle = paste0(
            "Comparación del efecto bruto con el efecto parcial ",
            "al incorporar expectativas a 12 meses"
        ),
        x = "Fecha final de la ventana",
        y = expression(theta[1*w]),
        linetype = "Modelo"
    ) +
    theme_minimal()

ggsave(
    paste0(
        "output/",
        nombre_base,
        "_anidados_theta1.png"
    ),
    plot = p_anidados_theta1,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_anidados_theta1)

# ------------------------------------------------------------
# 11) Modelos anidados: expectativas a 12 meses
# ------------------------------------------------------------

comparacion_theta2 <- rolling_results |>
    select(
        fecha_fin,
        theta2_m2,
        theta2_m3_ols
    ) |>
    pivot_longer(
        cols = c(
            theta2_m2,
            theta2_m3_ols
        ),
        names_to = "modelo",
        values_to = "theta2"
    ) |>
    mutate(
        modelo = recode(
            modelo,
            theta2_m2 =
                "M2: solo expectativas a 12 meses",
            theta2_m3_ols =
                "M3: modelo completo"
        )
    )

p_anidados_theta2 <- ggplot(
    comparacion_theta2,
    aes(
        x = fecha_fin,
        y = theta2,
        linetype = modelo
    )
) +
    geom_line(
        linewidth = 0.9
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    labs(
        title = paste0(
            "Modelos anidados: expectativas a 12 meses ",
            "(ventana de ",
            window_size,
            " meses)"
        ),
        subtitle = paste0(
            "Comparación del efecto bruto con el efecto parcial ",
            "al incorporar la inflación"
        ),
        x = "Fecha final de la ventana",
        y = expression(theta[2*w]),
        linetype = "Modelo"
    ) +
    theme_minimal()

ggsave(
    paste0(
        "output/",
        nombre_base,
        "_anidados_theta2.png"
    ),
    plot = p_anidados_theta2,
    width = 10,
    height = 5,
    dpi = 300
)

print(p_anidados_theta2)

# ------------------------------------------------------------
# 12) Resumen en consola
# ------------------------------------------------------------

cat(
    "\n============================================================\n"
)

cat(
    "RESULTADOS DEL MODELO STROHSAL RESTRINGIDO\n"
)

cat(
    "============================================================\n"
)

cat(
    "Ventana móvil: ",
    window_size,
    " meses\n"
)

cat(
    "Bootstrap: ",
    B_boot,
    " réplicas por ventana\n"
)

cat(
    "Longitud del bloque: ",
    block_length,
    " meses\n"
)

cat(
    "Número de ventanas: ",
    nrow(rolling_results),
    "\n"
)

cat(
    "Promedio del indicador restringido de anclaje: ",
    round(
        mean(
            rolling_results$anclaje_restringido,
            na.rm = TRUE
        ),
        3
    ),
    "\n"
)

cat(
    "Theta 1 en la frontera cero: ",
    sum(
        rolling_results$theta1_frontera,
        na.rm = TRUE
    ),
    " ventanas\n"
)

cat(
    "Theta 2 en la frontera cero: ",
    sum(
        rolling_results$theta2_frontera,
        na.rm = TRUE
    ),
    " ventanas\n"
)

cat(
    "Anclaje en la frontera cero: ",
    sum(
        rolling_results$anclaje_frontera,
        na.rm = TRUE
    ),
    " ventanas\n"
)

cat(
    "Penalización media del SSE por las restricciones (%): ",
    round(
        mean(
            rolling_results$penalizacion_sse_pct,
            na.rm = TRUE
        ),
        3
    ),
    "\n"
)

cat(
    "Correlación media entre los regresores: ",
    round(
        mean(
            rolling_results$correlacion_x1_x2,
            na.rm = TRUE
        ),
        3
    ),
    "\n"
)

cat(
    "\nArchivos generados en la carpeta output:\n",
    "- Base completa de resultados\n",
    "- Resumen de restricciones activas\n",
    "- Indicador principal con intervalo bootstrap\n",
    "- Componentes restringidos con intervalos bootstrap\n",
    "- Comparación restringido vs. irrestricto\n",
    "- Comparación de modelos anidados para theta 1\n",
    "- Comparación de modelos anidados para theta 2\n"
)

cat(
    "============================================================\n"
)
