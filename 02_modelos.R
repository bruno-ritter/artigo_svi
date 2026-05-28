# ==============================================================================
# 02_modelos.R
# TCC: Previsao de Volatilidade com SVI (Google Trends)
# M0 (Random Walk), M1 (GARCH(1,1)), M2 (GARCH-X com svi_log_dev)
# Janela movel de 156 semanas; painel balanceado por ticker
# ==============================================================================

# ── 0. Bibliotecas e constantes ───────────────────────────────────────────────
.libPaths("~/R/library")
library(rugarch)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(parallel)
library(tseries)

CAMINHO         <- "./data/"
JANELA          <- 52L * 3L                          # 156 semanas de treino
T_START         <- JANELA + 2L
USE_PARALLEL    <- TRUE
N_CORES         <- max(1L, parallel::detectCores() - 2L)
TIMEOUT_TICKER  <- 30L * 60L                         # 30 min/ticker em segundos
DIR_CHECKPOINT  <- paste0(CAMINHO, "checkpoints/")   # 1 .rds por ticker
FORCE_RECOMPUTE <- FALSE                             # se TRUE, ignora checkpoints
set.seed(42)

dir.create(DIR_CHECKPOINT, showWarnings = FALSE, recursive = TRUE)

cat("=== 02_modelos.R ===\n")
cat(sprintf("Janela de treino : %d semanas\n", JANELA))
cat(sprintf("Primeira prev.   : indice t = %d\n", T_START))
cat(sprintf("Nucleos usados   : %d\n", if (USE_PARALLEL) N_CORES else 1L))
cat(sprintf("Timeout/ticker   : %d min\n", TIMEOUT_TICKER %/% 60L))
cat(sprintf("Checkpoints      : %s\n", DIR_CHECKPOINT))

# ── 1. Carrega dados ──────────────────────────────────────────────────────────
df_final <- read_csv(
  paste0(CAMINHO, "base_final_tcc.csv"),
  col_types = cols(semana = col_date(format = "%Y-%m-%d")),
  show_col_types = FALSE
)

# Filtra apenas os tickers de interesse
# TICKERS_SELECIONADOS <- c("PETR4", "GFSA3", "CPLE3", "CVCB3", "BEES3", "TECN3", "PATI4")  # <- edite conforme necessário
# df_final <- df_final |> filter(ticker_b3 %in% TICKERS_SELECIONADOS)

cat(sprintf("\nbase_final_tcc   : %d linhas, %d tickers\n",
            nrow(df_final), n_distinct(df_final$ticker_b3)))

# ── 2. Funcoes auxiliares ─────────────────────────────────────────────────────

# Verifica se a janela de treino tem dados completos para a semana t
janela_ok <- function(df_t, t, janela = JANELA) {
  w <- (t - janela):(t - 1L)
  w_svi_lag <- (t - janela - 1L):(t - 2L)
  all(!is.na(df_t$ret_semanal[w])) &&
    all(!is.na(df_t$svi_log_dev[w_svi_lag])) &&
    !is.na(df_t$svi_log_dev[t - 1L]) &&
    !is.na(df_t$rv[t - 1L])
}

# Diagnosticos dos residuos padronizados z_t = eps_t / sigma_t
# Ljung-Box em z (autocorrelacao na media),
# Ljung-Box em z^2 (efeitos ARCH remanescentes),
# Jarque-Bera (normalidade)
diag_residuos <- function(fit, modelo) {
  if (is.null(fit)) return(list(resumo = NULL, z = NULL))

  z <- tryCatch(as.numeric(residuals(fit, standardize = TRUE)),
                error = function(e) NULL)
  if (is.null(z)) return(list(resumo = NULL, z = NULL))

  z <- z[is.finite(z)]
  n <- length(z)
  if (n < 30L) return(list(resumo = NULL, z = NULL))

  lb_z  <- tryCatch(Box.test(z,   lag = 10L, type = "Ljung-Box"),
                    error = function(e) NULL)
  lb_z2 <- tryCatch(Box.test(z^2, lag = 10L, type = "Ljung-Box"),
                    error = function(e) NULL)
  jb    <- tryCatch(tseries::jarque.bera.test(z),
                    error = function(e) NULL)
  if (is.null(lb_z) || is.null(lb_z2) || is.null(jb))
    return(list(resumo = NULL, z = NULL))

  m <- mean(z)
  v <- mean((z - m)^2)

  resumo <- tibble(
    modelo       = modelo,
    n_obs        = n,
    media_z      = m,
    sd_z         = sd(z),
    assimetria   = mean((z - m)^3) / v^1.5,
    curtose      = mean((z - m)^4) / v^2,
    lb_z_stat    = unname(lb_z$statistic),
    lb_z_pvalor  = unname(lb_z$p.value),
    lb_z2_stat   = unname(lb_z2$statistic),
    lb_z2_pvalor = unname(lb_z2$p.value),
    jb_stat      = unname(jb$statistic),
    jb_pvalor    = unname(jb$p.value)
  )

  list(resumo = resumo, z = z)
}

# M0: passeio aleatorio (previsao = rv da semana anterior)
calc_m0 <- function(df_t, idx) {
  tibble(
    semana     = df_t$semana[idx],
    rv_pred_m0 = df_t$rv[idx - 1L]
  )
}

# M1: GARCH(1,1) com janela movel
# Retornos escalados por 100 para estabilidade numerica; previsao reescalada
calc_m1 <- function(df_t, idx, janela = JANELA) {
  spec <- ugarchspec(
    variance.model     = list(model = "eGARCH", garchOrder = c(1L, 1L)),
    mean.model         = list(armaOrder = c(0L, 0L), include.mean = TRUE),
    distribution.model = "norm"
  )

  resultados <- vector("list", length(idx))
  ultimo_fit <- NULL

  for (k in seq_along(idx)) {
    t         <- idx[k]
    train_ret <- df_t$ret_semanal[(t - janela):(t - 1L)] * 100

    linha <- list(
      semana     = df_t$semana[t],
      rv         = df_t$rv[t],
      rv_pred_m1 = NA_real_,
      omega      = NA_real_,
      alpha      = NA_real_,
      beta       = NA_real_,
      status_m1  = NA_character_
    )

    fit_res <- tryCatch({
      fit <- ugarchfit(spec, data = train_ret, solver = "hybrid",
                       fit.control = list(stationarity = 1),
                       solver.control = list(trace = 0))
      fc  <- ugarchforecast(fit, n.ahead = 1L)
      list(ok = TRUE, fit = fit, fc = fc)
    }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))

    if (fit_res$ok) {
      mc               <- fit_res$fit@fit$matcoef
      linha$rv_pred_m1 <- as.numeric(sigma(fit_res$fc))^2 / 1e4
      linha$omega      <- mc["omega",  1L]
      linha$alpha      <- mc["alpha1", 1L]
      linha$beta       <- mc["beta1",  1L]
      linha$status_m1  <- "ok"
      ultimo_fit       <- fit_res$fit
    } else {
      linha$status_m1 <- fit_res$msg
    }

    resultados[[k]] <- linha
  }

  out <- bind_rows(resultados)
  attr(out, "ultimo_fit") <- ultimo_fit
  out
}

# M2: GARCH-X(1,1) com svi_log_dev como regressor externo na eq. da variancia
# Usa svi_log_dev defasado em uma semana: retorno da semana s explicado por SVI de s - 1.
calc_m2 <- function(df_t, idx, janela = JANELA) {
  resultados <- vector("list", length(idx))
  ultimo_fit <- NULL

  for (k in seq_along(idx)) {
    t         <- idx[k]
    train_ret <- df_t$ret_semanal[(t - janela):(t - 1L)] * 100
    train_x   <- matrix(df_t$svi_log_dev[(t - janela - 1L):(t - 2L)], ncol = 1L)
    fc_x      <- matrix(df_t$svi_log_dev[t - 1L], nrow = 1L, ncol = 1L)

    linha <- list(
      semana     = df_t$semana[t],
      rv         = df_t$rv[t],
      rv_pred_m2 = NA_real_,
      alpha_m2   = NA_real_,
      beta_m2    = NA_real_,
      gamma      = NA_real_,
      gamma_pval = NA_real_,
      status_m2  = NA_character_
    )

    warnings_m2 <- character(0)
    fit_res <- tryCatch(
      withCallingHandlers({
        spec <- ugarchspec(
          variance.model     = list(model = "eGARCH", garchOrder = c(1L, 1L),
                                    external.regressors = train_x),
          mean.model         = list(armaOrder = c(0L, 0L), include.mean = TRUE),
          distribution.model = "norm"
        )
        fit <- ugarchfit(spec, data = train_ret, solver = "hybrid",
                         fit.control = list(stationarity = 1),
                         solver.control = list(trace = 0))
        fc  <- ugarchforecast(fit, n.ahead = 1L,
                              external.forecasts = list(vregfor = fc_x))
        list(ok = TRUE, fit = fit, fc = fc, warnings = warnings_m2)
      }, warning = function(w) {
        msg <- gsub("[\r\n]+", " ", conditionMessage(w))
        warnings_m2 <<- c(warnings_m2, trimws(msg))
        invokeRestart("muffleWarning")
      }),
      error = function(e) list(ok = FALSE, msg = conditionMessage(e), warnings = warnings_m2)
    )

    if (fit_res$ok) {
      mc               <- fit_res$fit@fit$matcoef
      linha$rv_pred_m2 <- as.numeric(sigma(fit_res$fc))^2 / 1e4
      linha$alpha_m2   <- mc["alpha1", 1L]
      linha$beta_m2    <- mc["beta1",  1L]
      linha$gamma      <- mc["vxreg1", 1L]
      linha$gamma_pval <- mc["vxreg1", 4L]
      linha$status_m2  <- if (length(fit_res$warnings) == 0L) {
        "ok"
      } else {
        paste("ok_warning:", paste(unique(fit_res$warnings), collapse = " | "))
      }
      ultimo_fit       <- fit_res$fit
    } else {
      linha$status_m2 <- fit_res$msg
    }

    resultados[[k]] <- linha
  }

  out <- bind_rows(resultados)
  attr(out, "ultimo_fit") <- ultimo_fit
  out
}

# ── 3. Pre-voo: calendario global e tickers elegiveis ────────────────────────

tickers_all <- sort(unique(df_final$ticker_b3))
grupo_ref <- df_final %>%
  filter(ticker_b3 == tickers_all[1L], !is.na(svi_log_dev)) %>%
  arrange(semana)
n_ref <- nrow(grupo_ref)

if (n_ref <= T_START)
  stop(sprintf("Serie de referencia muito curta (n=%d, T_START=%d)", n_ref, T_START))

idx_global          <- seq(T_START, n_ref)
semanas_previsao    <- grupo_ref$semana[idx_global]
n_semanas_previsao  <- length(idx_global)

# Filtra tickers com painel completo e sem NAs nas janelas de previsao
tickers_elegiveis <- character(0)
for (tk in tickers_all) {
  g <- df_final %>%
    filter(ticker_b3 == tk, !is.na(svi_log_dev)) %>%
    arrange(semana)
  if (nrow(g) != n_ref || !identical(g$semana, grupo_ref$semana)) next
  if (all(vapply(idx_global, janela_ok, logical(1L), df_t = g)))
    tickers_elegiveis <- c(tickers_elegiveis, tk)
}

cat("\nPre-voo:\n")
cat(sprintf("  Semanas de previsao : %d (%s a %s)\n",
            n_semanas_previsao, min(semanas_previsao), max(semanas_previsao)))
cat(sprintf("  Tickers elegiveis   : %d / %d\n",
            length(tickers_elegiveis), length(tickers_all)))

if (length(tickers_elegiveis) == 0L) stop("Nenhum ticker passou no pre-voo.")

# ── 4. Loop principal por ticker ──────────────────────────────────────────────

processar_ticker <- function(i) {
  ticker      <- tickers_elegiveis[i]
  arquivo_ckp <- file.path(DIR_CHECKPOINT, paste0(ticker, ".rds"))

  # Reaproveita resultado salvo (retomada incremental)
  if (!FORCE_RECOMPUTE && file.exists(arquivo_ckp)) {
    return(readRDS(arquivo_ckp))
  }

  t_inicio <- Sys.time()
  setTimeLimit(elapsed = TIMEOUT_TICKER, transient = TRUE)

  resultado <- tryCatch({
    grupo <- df_final %>%
      filter(ticker_b3 == ticker, !is.na(svi_log_dev)) %>%
      arrange(semana)

    df_m0 <- calc_m0(grupo, idx_global)
    df_m1 <- calc_m1(grupo, idx_global)
    df_m2 <- calc_m2(grupo, idx_global)

    stopifnot(identical(df_m1$semana, df_m2$semana))

    df_tick <- df_m1 %>%
      mutate(rv_pred_m0 = df_m0$rv_pred_m0) %>%
      left_join(
        df_m2 %>% select(semana, rv_pred_m2, alpha_m2, beta_m2,
                          gamma, gamma_pval, status_m2),
        by = "semana"
      ) %>%
      mutate(ticker_b3 = ticker)

    resumo <- tibble(
      ticker_b3        = ticker,
      n_previsoes      = nrow(df_m1),
      gamma_medio      = mean(df_m2$gamma, na.rm = TRUE),
      gamma_pval_medio = mean(df_m2$gamma_pval, na.rm = TRUE),
      pct_signif_5pct  = mean(df_m2$gamma_pval < 0.05, na.rm = TRUE) * 100,
      persistencia_m1  = mean(df_m1$alpha + df_m1$beta, na.rm = TRUE),
      persistencia_m2  = mean(df_m2$alpha_m2 + df_m2$beta_m2, na.rm = TRUE)
    )

    diag_m1 <- diag_residuos(attr(df_m1, "ultimo_fit"), "m1")
    diag_m2 <- diag_residuos(attr(df_m2, "ultimo_fit"), "m2")

    diagnosticos <- bind_rows(diag_m1$resumo, diag_m2$resumo)
    if (nrow(diagnosticos) > 0L)
      diagnosticos <- mutate(diagnosticos, ticker_b3 = ticker, .before = 1L)

    z_lista <- list()
    if (!is.null(diag_m1$z))
      z_lista[["m1"]] <- tibble(modelo = "m1", i = seq_along(diag_m1$z), z = diag_m1$z)
    if (!is.null(diag_m2$z))
      z_lista[["m2"]] <- tibble(modelo = "m2", i = seq_along(diag_m2$z), z = diag_m2$z)

    residuos_z <- bind_rows(z_lista)
    if (nrow(residuos_z) > 0L)
      residuos_z <- mutate(residuos_z, ticker_b3 = ticker, .before = 1L)

    list(
      ticker_b3    = ticker,
      status       = "ok",
      tempo_s      = as.numeric(difftime(Sys.time(), t_inicio, units = "secs")),
      previsoes    = df_tick,
      resumo       = resumo,
      diagnosticos = diagnosticos,
      residuos_z   = residuos_z
    )
  },
  error = function(e) {
    msg <- conditionMessage(e)
    eh_timeout <- grepl("reached elapsed time limit|reached CPU time limit", msg)
    list(
      ticker_b3    = ticker,
      status       = if (eh_timeout) "timeout" else paste0("erro: ", msg),
      tempo_s      = as.numeric(difftime(Sys.time(), t_inicio, units = "secs")),
      previsoes    = NULL,
      resumo       = NULL,
      diagnosticos = NULL,
      residuos_z   = NULL
    )
  })

  setTimeLimit()  # remove o limite antes de devolver controle ao master

  # Grava checkpoint mesmo em falha
  saveRDS(resultado, arquivo_ckp)
  resultado
}

cat(sprintf("\nProcessando %d tickers elegiveis...\n", length(tickers_elegiveis)))

if (USE_PARALLEL && N_CORES > 1L) {
  cl <- parallel::makeCluster(N_CORES, type = "PSOCK")
  parallel::clusterExport(cl, varlist = c(
    "df_final", "JANELA", "T_START", "idx_global",
    "tickers_elegiveis", "processar_ticker",
    "calc_m0", "calc_m1", "calc_m2", "janela_ok", "diag_residuos",
    "DIR_CHECKPOINT", "TIMEOUT_TICKER", "FORCE_RECOMPUTE"
  ), envir = environment())
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      .libPaths("~/R/library")
      library(rugarch); library(dplyr); library(tseries)
    })
  })
  # chunk.size = 1L garante load balancing real: 1 ticker por vez, sem pre-bloco
  results_list <- tryCatch(
    parallel::parLapplyLB(cl, seq_along(tickers_elegiveis),
                          processar_ticker, chunk.size = 1L),
    finally = parallel::stopCluster(cl)
  )
} else {
  results_list <- lapply(seq_along(tickers_elegiveis), processar_ticker)
}

# ── 5. Consolidacao, validacao e exportacao ───────────────────────────────────

results_ok  <- Filter(Negate(is.null), results_list)
status_runs <- vapply(results_ok, `[[`, character(1L), "status")
tempo_runs  <- vapply(results_ok, `[[`, numeric(1L),   "tempo_s")
ticker_runs <- vapply(results_ok, `[[`, character(1L), "ticker_b3")

falhas <- which(status_runs != "ok")
if (length(falhas) > 0L) {
  cat(sprintf("\n%d ticker(s) com falha/timeout (excluidos do painel):\n",
              length(falhas)))
  for (k in falhas)
    cat(sprintf("  %-8s  %s  (%.1f s)\n",
                ticker_runs[k], status_runs[k], tempo_runs[k]))
}

# Painel apenas com tickers que completaram com sucesso
tem_previsoes    <- vapply(results_ok, function(x) !is.null(x$previsoes),
                           logical(1L))
results_panel    <- results_ok[tem_previsoes]
df_previsoes     <- bind_rows(lapply(results_panel, `[[`, "previsoes"))
df_resumo_params <- bind_rows(lapply(results_panel, `[[`, "resumo"))
df_diagnosticos  <- bind_rows(lapply(results_panel, `[[`, "diagnosticos"))
df_residuos_z    <- bind_rows(lapply(results_panel, `[[`, "residuos_z"))

# Valida painel balanceado entre tickers sobreviventes
n_por_ticker <- df_previsoes %>% count(ticker_b3)
if (length(unique(n_por_ticker$n)) != 1L)
  stop("Painel desbalanceado: tickers com numero diferente de semanas.")
if (unique(n_por_ticker$n) != n_semanas_previsao)
  stop(sprintf("Esperado %d semanas/ticker, obtido %d.",
               n_semanas_previsao, unique(n_por_ticker$n)))

# Ordena colunas e salva
df_previsoes <- df_previsoes %>%
  select(ticker_b3, semana, rv,
         rv_pred_m0,
         rv_pred_m1, omega, alpha, beta, status_m1,
         rv_pred_m2, alpha_m2, beta_m2, gamma, gamma_pval, status_m2)

diagnostico_previsoes <- df_previsoes %>%
  filter(!is.na(rv_pred_m1), !is.na(rv_pred_m2)) %>%
  summarise(
    n_validas        = n(),
    iguais_exatas    = sum(rv_pred_m1 == rv_pred_m2),
    iguais_6_casas   = sum(round(rv_pred_m1, 6) == round(rv_pred_m2, 6)),
    iguais_8_casas   = sum(round(rv_pred_m1, 8) == round(rv_pred_m2, 8)),
    iguais_10_casas  = sum(round(rv_pred_m1, 10) == round(rv_pred_m2, 10)),
    gamma_abs_lt_1e_6 = sum(abs(gamma) < 1e-6, na.rm = TRUE)
  )

cat(sprintf(
  paste0(
    "Diagnostico M1/M2: %d validas, %d iguais exatas, ",
    "%d iguais em 6 casas, %d em 8 casas, %d em 10 casas, ",
    "%d com |gamma| < 1e-6\n"
  ),
  diagnostico_previsoes$n_validas,
  diagnostico_previsoes$iguais_exatas,
  diagnostico_previsoes$iguais_6_casas,
  diagnostico_previsoes$iguais_8_casas,
  diagnostico_previsoes$iguais_10_casas,
  diagnostico_previsoes$gamma_abs_lt_1e_6
))

write_csv(df_previsoes,     paste0(CAMINHO, "previsoes_consolidadas.csv"))
write_csv(df_resumo_params, paste0(CAMINHO, "resumo_parametros.csv"))
write_csv(df_diagnosticos,  paste0(CAMINHO, "diagnosticos_residuos.csv"))
write_csv(df_residuos_z,    paste0(CAMINHO, "residuos_padronizados.csv"))

cat(sprintf("\nprevisoes_consolidadas.csv : %d linhas, %d tickers, %d semanas/ticker\n",
            nrow(df_previsoes), n_distinct(df_previsoes$ticker_b3),
            unique(n_por_ticker$n)))
cat(sprintf("resumo_parametros.csv      : %d tickers\n", nrow(df_resumo_params)))
cat(sprintf("diagnosticos_residuos.csv  : %d linhas\n", nrow(df_diagnosticos)))
cat(sprintf("residuos_padronizados.csv  : %d linhas\n", nrow(df_residuos_z)))
