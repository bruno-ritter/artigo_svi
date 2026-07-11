# ==============================================================================
# 02_modelos.R
# Previsao de Volatilidade com SVI (Google Trends)
# Modelos:
#   M0  Random Walk (previsao = RV da semana anterior)
#   M1  eGARCH(1,1)            (assimetrico, puro)
#   M2  eGARCH-X(1,1)          (assimetrico + svi_log_dev)
#   M3  GARCH(1,1)  / sGARCH   (simetrico,  puro)
#   M4  GARCH-X(1,1)/ sGARCH-X (simetrico  + svi_log_dev)
# Objetivo da extensao (artigo): comparar o p-valor do coeficiente do SVI (delta)
# no GARCH-X (M4) vs eGARCH-X (M2). A hipotese e que delta perde relevancia quando
# o modelo captura a assimetria nativamente (efeito alavancagem, gamma1 do eGARCH),
# indicando que SVI e alavancagem disputam a mesma fonte de variacao na variancia.
# Janela movel de 156 semanas; painel balanceado por ticker.
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

CAMINHO <- "./data/"
JANELA <- 52L * 3L # 156 semanas de treino
T_START <- JANELA + 2L
USE_PARALLEL <- TRUE
N_CORES <- max(1L, parallel::detectCores() - 1L)
TIMEOUT_TICKER <- 30L * 60L # 30 min/ticker em segundos
# checkpoints_v2: estrutura de 4 modelos (M1-M4). Os checkpoints antigos (3 modelos)
# ficam preservados em checkpoints/ e nao sao reutilizados aqui.
DIR_CHECKPOINT <- paste0(CAMINHO, "checkpoints_v2/") # 1 .rds por ticker
FORCE_RECOMPUTE <- FALSE # se TRUE, ignora checkpoints
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
# TICKERS_SELECIONADOS <- c("PETR4", "GFSA3", "CPLE3", "CVCB3", "BEES3", "TECN3", "PATI4")
# df_final <- df_final |> filter(ticker_b3 %in% TICKERS_SELECIONADOS)

cat(sprintf(
  "\nbase_final_tcc   : %d linhas, %d tickers\n",
  nrow(df_final), n_distinct(df_final$ticker_b3)
))

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
  if (is.null(fit)) {
    return(list(resumo = NULL, z = NULL))
  }

  z <- tryCatch(as.numeric(residuals(fit, standardize = TRUE)),
    error = function(e) NULL
  )
  if (is.null(z)) {
    return(list(resumo = NULL, z = NULL))
  }

  z <- z[is.finite(z)]
  n <- length(z)
  if (n < 30L) {
    return(list(resumo = NULL, z = NULL))
  }

  lb_z <- tryCatch(Box.test(z, lag = 10L, type = "Ljung-Box"),
    error = function(e) NULL
  )
  lb_z2 <- tryCatch(Box.test(z^2, lag = 10L, type = "Ljung-Box"),
    error = function(e) NULL
  )
  jb <- tryCatch(tseries::jarque.bera.test(z),
    error = function(e) NULL
  )
  if (is.null(lb_z) || is.null(lb_z2) || is.null(jb)) {
    return(list(resumo = NULL, z = NULL))
  }

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

# ── 3. Funcao geral dos modelos GARCH ─────────────────────────────────────────
# Cobre M1-M4 conforme (model, use_svi):
#   ("eGARCH", FALSE) -> M1 eGARCH puro
#   ("eGARCH", TRUE ) -> M2 eGARCH-X
#   ("sGARCH", FALSE) -> M3 GARCH puro
#   ("sGARCH", TRUE ) -> M4 GARCH-X
# Retornos escalados por 100 para estabilidade numerica; previsao reescalada (/1e4).
# Colunas de coeficientes padronizadas (rugarch):
#   omega, alpha (alpha1), beta (beta1),
#   leverage = gamma1 (so no eGARCH; NA no sGARCH),
#   delta = vxreg1 = coeficiente do SVI (so quando use_svi = TRUE) e seu p-valor.
# use_svi usa svi_log_dev defasado 1 semana: o retorno da semana s e explicado pelo
# SVI de s-1 (evita look-ahead), coerente entre M2 e M4.
calc_garch <- function(df_t, idx, model = c("eGARCH", "sGARCH"),
                       use_svi = FALSE, janela = JANELA) {
  model <- match.arg(model)
  resultados <- vector("list", length(idx))
  ultimo_fit <- NULL

  for (k in seq_along(idx)) {
    t <- idx[k]
    train_ret <- df_t$ret_semanal[(t - janela):(t - 1L)] * 100

    if (use_svi) {
      train_x <- matrix(df_t$svi_log_dev[(t - janela - 1L):(t - 2L)], ncol = 1L)
      fc_x <- matrix(df_t$svi_log_dev[t - 1L], nrow = 1L, ncol = 1L)
    } else {
      train_x <- NULL
      fc_x <- NULL
    }

    linha <- list(
      semana     = df_t$semana[t],
      rv         = df_t$rv[t],
      rv_pred    = NA_real_,
      omega      = NA_real_,
      alpha      = NA_real_,
      beta       = NA_real_,
      leverage   = NA_real_, # gamma1 (assimetria/alavancagem), so no eGARCH
      delta      = NA_real_, # vxreg1 (coeficiente do SVI), so com use_svi
      delta_pval = NA_real_,
      status     = NA_character_
    )

    warns <- character(0)
    fit_res <- tryCatch(
      withCallingHandlers(
        {
          spec <- ugarchspec(
            variance.model = list(
              model = model, garchOrder = c(1L, 1L),
              external.regressors = train_x
            ),
            mean.model = list(armaOrder = c(0L, 0L), include.mean = TRUE),
            distribution.model = "norm"
          )
          fit <- ugarchfit(spec,
            data = train_ret, solver = "hybrid",
            fit.control = list(stationarity = 1),
            solver.control = list(trace = 0)
          )
          fc <- if (use_svi) {
            ugarchforecast(fit,
              n.ahead = 1L,
              external.forecasts = list(vregfor = fc_x)
            )
          } else {
            ugarchforecast(fit, n.ahead = 1L)
          }
          list(ok = TRUE, fit = fit, fc = fc, warnings = warns)
        },
        warning = function(w) {
          msg <- gsub("[\r\n]+", " ", conditionMessage(w))
          warns <<- c(warns, trimws(msg))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) list(ok = FALSE, msg = conditionMessage(e), warnings = warns)
    )

    if (fit_res$ok) {
      mc <- fit_res$fit@fit$matcoef
      rn <- rownames(mc)
      linha$rv_pred <- as.numeric(sigma(fit_res$fc))^2 / 1e4
      linha$omega <- mc["omega", 1L]
      linha$alpha <- mc["alpha1", 1L]
      linha$beta <- mc["beta1", 1L]
      if ("gamma1" %in% rn) linha$leverage <- mc["gamma1", 1L]
      if (use_svi && "vxreg1" %in% rn) {
        linha$delta <- mc["vxreg1", 1L]
        linha$delta_pval <- mc["vxreg1", 4L]
      }
      linha$status <- if (length(fit_res$warnings) == 0L) {
        "ok"
      } else {
        paste("ok_warning:", paste(unique(fit_res$warnings), collapse = " | "))
      }
      ultimo_fit <- fit_res$fit
    } else {
      linha$status <- fit_res$msg
    }

    resultados[[k]] <- linha
  }

  out <- bind_rows(resultados)
  attr(out, "ultimo_fit") <- ultimo_fit
  out
}

# ── 3b. Pre-voo: calendario global e tickers elegiveis ───────────────────────

tickers_all <- sort(unique(df_final$ticker_b3))
grupo_ref <- df_final %>%
  filter(ticker_b3 == tickers_all[1L], !is.na(svi_log_dev)) %>%
  arrange(semana)
n_ref <- nrow(grupo_ref)

if (n_ref <= T_START) {
  stop(sprintf("Serie de referencia muito curta (n=%d, T_START=%d)", n_ref, T_START))
}

idx_global <- seq(T_START, n_ref)
semanas_previsao <- grupo_ref$semana[idx_global]
n_semanas_previsao <- length(idx_global)

# Filtra tickers com painel completo e sem NAs nas janelas de previsao
tickers_elegiveis <- character(0)
for (tk in tickers_all) {
  g <- df_final %>%
    filter(ticker_b3 == tk, !is.na(svi_log_dev)) %>%
    arrange(semana)
  if (nrow(g) != n_ref || !identical(g$semana, grupo_ref$semana)) next
  if (all(vapply(idx_global, janela_ok, logical(1L), df_t = g))) {
    tickers_elegiveis <- c(tickers_elegiveis, tk)
  }
}

cat("\nPre-voo:\n")
cat(sprintf(
  "  Semanas de previsao : %d (%s a %s)\n",
  n_semanas_previsao, min(semanas_previsao), max(semanas_previsao)
))
cat(sprintf(
  "  Tickers elegiveis   : %d / %d\n",
  length(tickers_elegiveis), length(tickers_all)
))

if (length(tickers_elegiveis) == 0L) stop("Nenhum ticker passou no pre-voo.")

# ── 4. Loop principal por ticker ──────────────────────────────────────────────

processar_ticker <- function(i) {
  ticker <- tickers_elegiveis[i]
  arquivo_ckp <- file.path(DIR_CHECKPOINT, paste0(ticker, ".rds"))

  # Reaproveita resultado salvo (retomada incremental)
  if (!FORCE_RECOMPUTE && file.exists(arquivo_ckp)) {
    return(readRDS(arquivo_ckp))
  }

  t_inicio <- Sys.time()
  setTimeLimit(elapsed = TIMEOUT_TICKER, transient = TRUE)

  resultado <- tryCatch(
    {
      grupo <- df_final %>%
        filter(ticker_b3 == ticker, !is.na(svi_log_dev)) %>%
        arrange(semana)

      df_m0 <- calc_m0(grupo, idx_global)
      df_m1 <- calc_garch(grupo, idx_global, model = "eGARCH", use_svi = FALSE)
      df_m2 <- calc_garch(grupo, idx_global, model = "eGARCH", use_svi = TRUE)
      df_m3 <- calc_garch(grupo, idx_global, model = "sGARCH", use_svi = FALSE)
      df_m4 <- calc_garch(grupo, idx_global, model = "sGARCH", use_svi = TRUE)

      stopifnot(
        identical(df_m1$semana, df_m2$semana),
        identical(df_m1$semana, df_m3$semana),
        identical(df_m1$semana, df_m4$semana)
      )

      # Painel por semana. Mantem o esquema legado do M0/M1/M2 (colunas usadas
      # pelo notebook 03: gamma == delta do M2 = vxreg1) e acrescenta M3/M4 e
      # a alavancagem (leverage = gamma1) do eGARCH.
      df_tick <- tibble(
        ticker_b3   = ticker,
        semana      = df_m1$semana,
        rv          = df_m1$rv,
        rv_pred_m0  = df_m0$rv_pred_m0,
        # M1 eGARCH puro
        rv_pred_m1  = df_m1$rv_pred,
        omega       = df_m1$omega,
        alpha       = df_m1$alpha,
        beta        = df_m1$beta,
        leverage_m1 = df_m1$leverage,
        status_m1   = df_m1$status,
        # M2 eGARCH-X
        rv_pred_m2  = df_m2$rv_pred,
        alpha_m2    = df_m2$alpha,
        beta_m2     = df_m2$beta,
        leverage_m2 = df_m2$leverage,
        gamma       = df_m2$delta, # legado: gamma = coeficiente do SVI (vxreg1)
        gamma_pval  = df_m2$delta_pval,
        status_m2   = df_m2$status,
        # M3 GARCH puro (simetrico)
        rv_pred_m3  = df_m3$rv_pred,
        omega_m3    = df_m3$omega,
        alpha_m3    = df_m3$alpha,
        beta_m3     = df_m3$beta,
        status_m3   = df_m3$status,
        # M4 GARCH-X (simetrico)
        rv_pred_m4    = df_m4$rv_pred,
        alpha_m4      = df_m4$alpha,
        beta_m4       = df_m4$beta,
        delta_m4      = df_m4$delta, # coeficiente do SVI no GARCH-X simetrico
        delta_pval_m4 = df_m4$delta_pval,
        status_m4     = df_m4$status
      )

      resumo <- tibble(
        ticker_b3        = ticker,
        n_previsoes      = nrow(df_m1),
        # --- legado (notebook 03): SVI no eGARCH-X (M2) ---
        gamma_medio      = mean(df_m2$delta, na.rm = TRUE),
        gamma_pval_medio = mean(df_m2$delta_pval, na.rm = TRUE),
        pct_signif_5pct  = mean(df_m2$delta_pval < 0.05, na.rm = TRUE) * 100,
        persistencia_m1  = mean(df_m1$alpha + df_m1$beta, na.rm = TRUE),
        persistencia_m2  = mean(df_m2$alpha + df_m2$beta, na.rm = TRUE),
        # --- delta do SVI: eGARCH-X (M2) vs GARCH-X (M4) ---
        delta_egarch_medio      = mean(df_m2$delta, na.rm = TRUE),
        delta_egarch_pval_medio = mean(df_m2$delta_pval, na.rm = TRUE),
        pct_signif_egarch_5pct  = mean(df_m2$delta_pval < 0.05, na.rm = TRUE) * 100,
        delta_garch_medio       = mean(df_m4$delta, na.rm = TRUE),
        delta_garch_pval_medio  = mean(df_m4$delta_pval, na.rm = TRUE),
        pct_signif_garch_5pct   = mean(df_m4$delta_pval < 0.05, na.rm = TRUE) * 100,
        # --- alavancagem (gamma1) do eGARCH ---
        leverage_m1_medio = mean(df_m1$leverage, na.rm = TRUE),
        leverage_m2_medio = mean(df_m2$leverage, na.rm = TRUE),
        # --- persistencias dos simetricos ---
        persistencia_m3  = mean(df_m3$alpha + df_m3$beta, na.rm = TRUE),
        persistencia_m4  = mean(df_m4$alpha + df_m4$beta, na.rm = TRUE)
      )

      diag_lista <- list(
        diag_residuos(attr(df_m1, "ultimo_fit"), "m1"),
        diag_residuos(attr(df_m2, "ultimo_fit"), "m2"),
        diag_residuos(attr(df_m3, "ultimo_fit"), "m3"),
        diag_residuos(attr(df_m4, "ultimo_fit"), "m4")
      )

      diagnosticos <- bind_rows(lapply(diag_lista, `[[`, "resumo"))
      if (nrow(diagnosticos) > 0L) {
        diagnosticos <- mutate(diagnosticos, ticker_b3 = ticker, .before = 1L)
      }

      z_lista <- list()
      for (d in diag_lista) {
        if (!is.null(d$z) && !is.null(d$resumo)) {
          mdl <- d$resumo$modelo[1L]
          z_lista[[mdl]] <- tibble(modelo = mdl, i = seq_along(d$z), z = d$z)
        }
      }

      residuos_z <- bind_rows(z_lista)
      if (nrow(residuos_z) > 0L) {
        residuos_z <- mutate(residuos_z, ticker_b3 = ticker, .before = 1L)
      }

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
    }
  )

  setTimeLimit() # remove o limite antes de devolver controle ao master

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
    "calc_m0", "calc_garch", "janela_ok", "diag_residuos",
    "DIR_CHECKPOINT", "TIMEOUT_TICKER", "FORCE_RECOMPUTE"
  ), envir = environment())
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      .libPaths("~/R/library")
      library(rugarch)
      library(dplyr)
      library(tseries)
    })
  })
  # chunk.size = 1L garante load balancing real: 1 ticker por vez, sem pre-bloco
  results_list <- tryCatch(
    parallel::parLapplyLB(cl, seq_along(tickers_elegiveis),
      processar_ticker,
      chunk.size = 1L
    ),
    finally = parallel::stopCluster(cl)
  )
} else {
  results_list <- lapply(seq_along(tickers_elegiveis), processar_ticker)
}

# ── 5. Consolidacao, validacao e exportacao ───────────────────────────────────

results_ok <- Filter(Negate(is.null), results_list)
status_runs <- vapply(results_ok, `[[`, character(1L), "status")
tempo_runs <- vapply(results_ok, `[[`, numeric(1L), "tempo_s")
ticker_runs <- vapply(results_ok, `[[`, character(1L), "ticker_b3")

falhas <- which(status_runs != "ok")
if (length(falhas) > 0L) {
  cat(sprintf(
    "\n%d ticker(s) com falha/timeout (excluidos do painel):\n",
    length(falhas)
  ))
  for (k in falhas) {
    cat(sprintf(
      "  %-8s  %s  (%.1f s)\n",
      ticker_runs[k], status_runs[k], tempo_runs[k]
    ))
  }
}

# Painel apenas com tickers que completaram com sucesso
tem_previsoes <- vapply(
  results_ok, function(x) !is.null(x$previsoes),
  logical(1L)
)
results_panel <- results_ok[tem_previsoes]
df_previsoes <- bind_rows(lapply(results_panel, `[[`, "previsoes"))
df_resumo_params <- bind_rows(lapply(results_panel, `[[`, "resumo"))
df_diagnosticos <- bind_rows(lapply(results_panel, `[[`, "diagnosticos"))
df_residuos_z <- bind_rows(lapply(results_panel, `[[`, "residuos_z"))

# Valida painel balanceado entre tickers sobreviventes
n_por_ticker <- df_previsoes %>% count(ticker_b3)
if (length(unique(n_por_ticker$n)) != 1L) {
  stop("Painel desbalanceado: tickers com numero diferente de semanas.")
}
if (unique(n_por_ticker$n) != n_semanas_previsao) {
  stop(sprintf(
    "Esperado %d semanas/ticker, obtido %d.",
    n_semanas_previsao, unique(n_por_ticker$n)
  ))
}

# Ordena colunas e salva
df_previsoes <- df_previsoes %>%
  select(
    ticker_b3, semana, rv,
    rv_pred_m0,
    rv_pred_m1, omega, alpha, beta, leverage_m1, status_m1,
    rv_pred_m2, alpha_m2, beta_m2, leverage_m2, gamma, gamma_pval, status_m2,
    rv_pred_m3, omega_m3, alpha_m3, beta_m3, status_m3,
    rv_pred_m4, alpha_m4, beta_m4, delta_m4, delta_pval_m4, status_m4
  )

diagnostico_previsoes <- df_previsoes %>%
  filter(!is.na(rv_pred_m1), !is.na(rv_pred_m2)) %>%
  summarise(
    n_validas = n(),
    iguais_exatas = sum(rv_pred_m1 == rv_pred_m2),
    iguais_6_casas = sum(round(rv_pred_m1, 6) == round(rv_pred_m2, 6)),
    iguais_8_casas = sum(round(rv_pred_m1, 8) == round(rv_pred_m2, 8)),
    iguais_10_casas = sum(round(rv_pred_m1, 10) == round(rv_pred_m2, 10)),
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

# ── 5b. Tabela pooled: p-valor do SVI (delta) em GARCH-X vs eGARCH-X ──────────
# Agrega TODAS as janelas OOS (todos os tickers). E o resultado-chave do artigo:
# se o p-valor medio/mediano do delta sobe (e o % de janelas significativas cai)
# ao passar de GARCH-X (M4) para eGARCH-X (M2), o SVI perde relevancia quando o
# modelo passa a capturar a assimetria nativamente (efeito alavancagem).
pooled_delta <- tibble(
  modelo = c("GARCH-X (simetrico, M4)", "eGARCH-X (assimetrico, M2)"),
  n_janelas = c(
    sum(!is.na(df_previsoes$delta_pval_m4)),
    sum(!is.na(df_previsoes$gamma_pval))
  ),
  delta_medio = c(
    mean(df_previsoes$delta_m4, na.rm = TRUE),
    mean(df_previsoes$gamma, na.rm = TRUE)
  ),
  delta_mediano = c(
    median(df_previsoes$delta_m4, na.rm = TRUE),
    median(df_previsoes$gamma, na.rm = TRUE)
  ),
  pval_medio = c(
    mean(df_previsoes$delta_pval_m4, na.rm = TRUE),
    mean(df_previsoes$gamma_pval, na.rm = TRUE)
  ),
  pval_mediano = c(
    median(df_previsoes$delta_pval_m4, na.rm = TRUE),
    median(df_previsoes$gamma_pval, na.rm = TRUE)
  ),
  pct_janelas_signif_5pct = c(
    mean(df_previsoes$delta_pval_m4 < 0.05, na.rm = TRUE) * 100,
    mean(df_previsoes$gamma_pval < 0.05, na.rm = TRUE) * 100
  ),
  pct_janelas_signif_10pct = c(
    mean(df_previsoes$delta_pval_m4 < 0.10, na.rm = TRUE) * 100,
    mean(df_previsoes$gamma_pval < 0.10, na.rm = TRUE) * 100
  )
)

# % de tickers em que o SVI e significativo (5%) na MAIORIA das janelas
pct_tickers_signif_maioria <- df_resumo_params %>%
  summarise(
    garch_x = mean(pct_signif_garch_5pct >= 50, na.rm = TRUE) * 100,
    egarch_x = mean(pct_signif_egarch_5pct >= 50, na.rm = TRUE) * 100
  )

cat("\n=== SVI (delta): GARCH-X (M4) vs eGARCH-X (M2) ===\n")
print(as.data.frame(pooled_delta), digits = 4L)
cat(sprintf(
  "\n%% de tickers com SVI signif. (5%%) na maioria das janelas: GARCH-X = %.1f%%, eGARCH-X = %.1f%%\n",
  pct_tickers_signif_maioria$garch_x, pct_tickers_signif_maioria$egarch_x
))

write_csv(df_previsoes, paste0(CAMINHO, "previsoes_consolidadas.csv"))
write_csv(df_resumo_params, paste0(CAMINHO, "resumo_parametros.csv"))
write_csv(df_diagnosticos, paste0(CAMINHO, "diagnosticos_residuos.csv"))
write_csv(df_residuos_z, paste0(CAMINHO, "residuos_padronizados.csv"))
write_csv(pooled_delta, paste0(CAMINHO, "delta_svi_garchx_vs_egarchx.csv"))

cat(sprintf(
  "\nprevisoes_consolidadas.csv : %d linhas, %d tickers, %d semanas/ticker\n",
  nrow(df_previsoes), n_distinct(df_previsoes$ticker_b3),
  unique(n_por_ticker$n)
))
cat(sprintf("resumo_parametros.csv      : %d tickers\n", nrow(df_resumo_params)))
cat(sprintf("diagnosticos_residuos.csv  : %d linhas\n", nrow(df_diagnosticos)))
cat(sprintf("residuos_padronizados.csv  : %d linhas\n", nrow(df_residuos_z)))
cat(sprintf("delta_svi_garchx_vs_egarchx.csv : %d linhas\n", nrow(pooled_delta)))
