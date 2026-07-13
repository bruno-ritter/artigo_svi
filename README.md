# TCC — SVI e Volatilidade de Ações da B3

Trabalho de Conclusão de Curso em Economia.

**Questão de pesquisa:** O Search Volume Index (SVI) do Google Trends melhora a previsão de volatilidade de ações listadas na B3?

**Status:** pipeline completo (coleta, estimação e avaliação). 300 tickers estimados; 240 na comparação final do cenário principal (demais tickers/janelas removidos por erro de estimação no eGARCH-X). Além do trio de previsão (M0 Random Walk, M1 eGARCH, M2 eGARCH-X), o pipeline agora estima também os **simétricos M3 GARCH e M4 GARCH-X**, usados apenas para comparar o p-valor do coeficiente do ASVI ($\delta$) entre a especificação simétrica e a assimétrica. No MCS (alpha = 10%) do cenário principal o conjunto de confiança é {M1, M2} — o M2 (eGARCH-X) permanece (p-valor = 0.284) e o M0 (Random Walk) é excluído (p-valor = 0.016). Sob remoção de janelas com erro de Hessiana ou filtros de zeros do SVI, o M2 sai do conjunto (p-valor ≈ 0.08).

---

## Metodologia

### Variáveis

**Volatilidade Realizada (RV)**
Proxy semanal de volatilidade calculada como a soma dos quadrados dos retornos logarítmicos diários dentro de cada semana:

$$RV_t = \sum_{d \in t} r_d^2, \quad r_d = \ln\left(\frac{P_d}{P_{d-1}}\right)$$

**Search Volume Index (SVI)**
O SVI é o índice semanal normalizado (0-100) do Google Trends para o termo de busca correspondente ao ticker, com geo=BR. A partir dele constrói-se o **desvio logarítmico** (ASVI, na linha de DA, ENGELBERG e GAO, 2011) em relação à mediana móvel de 8 semanas passadas:

$$\text{svi-log-dev}_t = \ln(1 + \text{SVI}_t) - \ln\left(1 + \underset{k \in \{t-8,\ldots,t-1\}}{\mathrm{median}}(\text{SVI}_k)\right)$$

O regressor externo efetivamente usado nos modelos (M2 e M4) é a **forma exponencial** desse desvio — uma razão de atenção, sempre positiva e centrada em 1:

$$\text{svi-exp-dev}_t = \exp(\text{svi-log-dev}_t) = \frac{1 + \text{SVI}_t}{1 + \mathrm{median}_{8}(\text{SVI})}$$

O uso de $\ln(1+\cdot)$ (`log1p`) evita problemas quando SVI = 0. Apenas os tickers cujo `svi_exp_dev` passa no teste de estacionariedade ADF entram na estimação.

---

### Modelos Competidores

**Trio de previsão** (entram no MCS):

| # | Modelo | Equação da variância condicional |
|---|--------|----------------------------------|
| M0 | Random Walk | $\hat{\sigma}^2_t = RV_{t-1}$ |
| M1 | eGARCH(1,1) | $\ln\sigma^2_t = \omega + \alpha\,z_{t-1} + \xi\,(\lvert z_{t-1}\rvert - \mathbb{E}\lvert z_{t-1}\rvert) + \beta\,\ln\sigma^2_{t-1}$ |
| M2 | eGARCH-X(1,1) | $\ln\sigma^2_t = \omega + \alpha\,z_{t-1} + \xi\,(\lvert z_{t-1}\rvert - \mathbb{E}\lvert z_{t-1}\rvert) + \beta\,\ln\sigma^2_{t-1} + \delta\,\text{svi-exp-dev}_{t-1}$ |

onde $z_{t-1} = \varepsilon_{t-1}/\sigma_{t-1}$ são os resíduos padronizados.

M0 é o benchmark ingênuo. M1 é o eGARCH de Nelson (1991), que modela o **logaritmo** da variância: garante positividade sem restrições nos parâmetros e captura o efeito alavancagem (assimetria entre choques positivos e negativos) pelo coeficiente $\alpha$ (termo de sinal), enquanto $\xi$ mede o efeito de magnitude (sobre $\lvert z\rvert$) e $\beta$ a persistência. M2 estende M1 com o `svi_exp_dev` como regressor externo na equação de log-variância; o coeficiente $\delta$ indica se o interesse de busca ajuda a explicar a volatilidade.

---

### Extensão simétrica (M3, M4) e o $\delta$ do ASVI

Para investigar se o poder do SVI se confunde com o efeito de assimetria (alavancagem) que o eGARCH modela nativamente, o pipeline estima também dois modelos **simétricos** — o GARCH(1,1) de Bollerslev (M3) e sua versão com o ASVI (M4) — e compara o p-valor do coeficiente $\delta$ do ASVI entre **M4 (GARCH-X, simétrico)** e **M2 (eGARCH-X, assimétrico)**. Esses dois modelos **não** entram no MCS; servem apenas à inferência sobre $\delta$.

| # | Modelo | Equação da variância condicional |
|---|--------|----------------------------------|
| M3 | GARCH(1,1) | $\sigma^2_t = \omega + \alpha\,\varepsilon^2_{t-1} + \beta\,\sigma^2_{t-1}$ |
| M4 | GARCH-X(1,1) | $\sigma^2_t = \omega + \alpha\,\varepsilon^2_{t-1} + \beta\,\sigma^2_{t-1} + \delta\,\text{svi-exp-dev}_{t-1}$ |

**Estimativas médias** (janelas OOS com `status = ok`):

| Modelo | $\omega$ | $\alpha$ | $\beta$ | $\alpha+\beta$ | $\xi$ (magnitude) | $\delta$ (ASVI) | p-valor $\delta$ |
|--------|---------:|---------:|--------:|---------------:|------------------:|----------------:|-----------------:|
| M3 GARCH    | 2.780 | 0.070 | 0.877 | 0.947 | —     | —      | —     |
| M4 GARCH-X  | —     | 0.071 | 0.876 | 0.947 | —     | 0.423  | 0.953 |
| M1 eGARCH   | 1.624 | 0.039 | 0.505 | —     | 0.122 | —      | —     |
| M2 eGARCH-X | —     | 0.039 | 0.365 | —     | 0.138 | −0.264 | 0.267 |

> No eGARCH (M1/M2), $\alpha$ é o coeficiente do termo de **sinal** (efeito alavancagem, sobre $z_{t-1}$) e $\xi$ o do termo de **magnitude** (sobre $\lvert z_{t-1}\rvert$); a persistência é $\beta$. No GARCH simétrico (M3/M4), $\alpha$ é o coeficiente **ARCH** (sobre $\varepsilon^2$) e a persistência é $\alpha+\beta$.

**Comparação do $\delta$ (ASVI)** — *pooled* sobre todas as janelas OOS (`data/delta_svi_garchx_vs_egarchx.csv`):

| Especificação | $\delta$ médio | $\delta$ mediano | p-valor médio | p-valor mediano | % janelas signif. (5%) |
|---------------|---------------:|-----------------:|--------------:|----------------:|-----------------------:|
| GARCH-X (M4, simétrico)    | 0.423  | ≈ 0    | 0.953 | ≈ 1.000 | 1.5%  |
| eGARCH-X (M2, assimétrico) | −0.284 | −0.012 | 0.272 | 0.103   | 43.1% |

> **Nota metodológica — o p-valor do M4 não é comparável ao do M2.** No GARCH-X (variância em **nível**), a positividade da variância impõe uma restrição de não-negatividade sobre $\delta$; como o `svi_exp_dev` é sempre positivo e centrado em 1 (colinear com $\omega$), o $\delta$ do M4 empilha na fronteira (mediana ≈ 0, p-valor ≈ 1) na maioria dos tickers. Assim, a comparação direta dos p-valores **não constitui um teste limpo** da hipótese de redundância — e, nos números brutos, o $\delta$ é *mais* (não menos) significativo no eGARCH-X. A evidência robusta de que o SVI é redundante com a memória de volatilidade do eGARCH vem de outra via — co-movimento contemporâneo, precedência da volatilidade passada sobre a atenção e queda de ~97% do coeficiente preditivo do SVI ao controlar pela RV defasada —, documentada em `03_analise_resultados.ipynb`.

---

### Estimação em Janela Móvel (*Rolling Window*)

Os modelos são estimados e avaliados **fora da amostra** (*out-of-sample*) por janela móvel de tamanho fixo:

- **Janela de treino:** 156 semanas (3 anos)
- **Horizonte de previsão:** 1 passo à frente (semana seguinte)
- **Procedimento:** a cada semana *t*, estima-se o modelo nas 156 observações anteriores e gera-se uma previsão para *t*; a janela avança uma semana e repete-se o processo
- **Primeira previsão:** indice t = 158 na série já filtrada (`T_START = 156 + 2`, pois o regressor externo do M2 precisa de uma defasagem adicional; as 8 primeiras semanas com `svi_log_dev` = NaN por falta de histórico para a mediana são descartadas antes da estimação na coleta de dados)
- **Total de previsões:** 103 previsões *out-of-sample* por ticker
- **Painel alinhado:** M0, M1 e M2 são estimados nas **mesmas** semanas calendário por ticker, bloco contíguo, sem pular semanas. Tickers que falham em qualquer semana do bloco (dados incompletos) são excluídos no pré-voo; os restantes têm o **mesmo** número de previsões

Isso evita *look-ahead bias*, reproduz o ambiente real de previsão e garante comparabilidade justa entre modelos e entre ativos.

A estimação dos quatro modelos GARCH-family (GARCH, GARCH-X, eGARCH, eGARCH-X) é feita em R (`rugarch`) via `02_modelos.R`.

---

### Função de Perda — QLIKE

A avaliação pontual utiliza a função de perda **QLIKE** (Quasi-Likelihood):

$$\text{QLIKE}(RV_t, \hat{\sigma}^2_t) = \frac{RV_t}{\hat{\sigma}^2_t} - \ln\frac{RV_t}{\hat{\sigma}^2_t} - 1$$

QLIKE é assimétrica e penaliza mais fortemente previsões subestimadas, sendo considerada robusta para comparação de modelos de volatilidade (Patton, 2011). O valor médio de QLIKE por ticker e por modelo é armazenado em `data/qlike_por_ticker.csv`.

---

### Inferência Formal — Model Confidence Set (MCS)

A comparação entre os três modelos é feita pelo **Model Confidence Set** de Hansen, Lunde & Nason (2011), implementado via `arch.bootstrap.MCS`. O MCS é um teste de hipótese sequencial que, a partir de um conjunto inicial de modelos, elimina iterativamente o pior modelo enquanto ele for significativamente inferior aos demais (alpha = 10%). O conjunto final contém todos os modelos que *não podem ser rejeitados* como igualmente bons ao melhor.

Resultado salvo em `data/mcs_resultado.csv`.

---

### Filtro de Tickers

Tickers com alta proporção de semanas com SVI = 0 indicam que o Google Trends não registrou buscas suficientes para o papel. Isso causa dois problemas: erro de estimação (matriz singular / falha de convergência) no eGARCH-X e viés por ruído puro no regressor externo.

O repositório avalia **três estratégias de filtro** comparativamente (diagnóstico salvo em `data/diagnostico_filtros.txt`):

| Estratégia | Descrição | Tickers na comparação final |
|------------|-----------|-----------------------------|
| Sem filtro de zeros | Não remove nenhum ticker | 240 |
| Zeros >= 70% | Remove tickers com 70% ou mais das semanas com SVI = 0, além dos com erro | 162 |
| Zeros >= 50% | Remove tickers com 50% ou mais das semanas com SVI = 0, além dos com erro | 147 |

A análise principal usa a estratégia **sem filtro de zeros** (apenas remoção das janelas com erro de estimação no eGARCH-X), totalizando 240 tickers e 23.443 observações *out-of-sample*. O M1 (eGARCH puro) está no MCS em **todas** as estratégias; o M2 (eGARCH-X) só permanece no conjunto no cenário sem filtro (sai sob remoção de janelas com erro de Hessiana ou sob os filtros de zeros).

---

## Resultados Preliminares

| Modelo | QLIKE médio | QLIKE mediano | Melhor em N tickers | MCS (p-valor) |
|--------|-------------|---------------|---------------------|---------------|
| M0_RW | 66.765 | 1.677 | 5 | 0.016 (fora) |
| M1_GARCH | 14.527 | 0.703 | 165 | 1.000 (dentro) |
| M2_GARCHX | 11.447 | 0.767 | 70 | 0.284 (dentro) |

> Médias e medianas das médias de QLIKE por ticker (240 tickers, `data/qlike_por_ticker.csv`). O QLIKE **médio** é fortemente influenciado por alguns tickers com perdas extremas (por isso a média do M2 fica abaixo da do M1); a **mediana** por ticker — M1 (0.703) < M2 (0.767) < M0 (1.677) — e a contagem de vitórias por ticker confirmam de forma robusta a vantagem do M1.

O M1 (eGARCH puro) domina o ticker típico (melhor em 165 dos 240) e ancora o MCS (p-valor = 1.000). O M2 (eGARCH-X) permanece no conjunto de confiança no cenário principal (p-valor = 0.284), mas vence individualmente em apenas 70 dos 240 tickers e não supera o M1 — a inclusão do SVI não gera ganho preditivo robusto. O coeficiente $\delta$ do SVI é significativo (5%) em pelo menos metade das janelas para 98 dos 240 tickers, concentrando-se nos ativos mais voláteis. Os diagnósticos de resíduos padronizados (Ljung-Box em $z_t$ e $z_t^2$, Jarque-Bera) estão em `data/diagnosticos_residuos.csv`.

---

## Estrutura do repositório

```
.
├── 01_coleta_dados.ipynb        # Coleta e preparo dos dados
├── 02_modelos.R                 # Estimação M0-M4 (rugarch), painel alinhado
├── 03_analise_resultados.ipynb  # QLIKE, MCS, diagnósticos e análise exploratória
├── requirements.txt
└── data/
    ├── acoes-listadas-b3.csv           # Lista de tickers da B3
    ├── acoes_elegiveis.csv             # Tickers que passaram nos filtros de liquidez/histórico
    ├── precos_semanais.csv             # Preços e volatilidade realizada semanal por ticker
    ├── google_trends_svi.csv           # SVI semanal por ticker (Google Trends)
    ├── log_coleta_svi.csv              # Log da coleta do Google Trends
    ├── log_adf_svi.csv                 # Resultados do teste ADF por ticker
    ├── base_final_tcc.csv              # Base combinada (preços + SVI + tratamentos)
    ├── previsoes_consolidadas.csv      # Previsões OOS M0-M4 + coeficientes por janela
    ├── resumo_parametros.csv           # Parâmetros médios por ticker (M1-M4: alpha, beta, delta, leverage)
    ├── apendice_resumo_params.csv      # Resumo arredondado para apêndice
    ├── qlike_por_ticker.csv            # QLIKE médio por ticker e modelo (240 tickers)
    ├── mcs_resultado.csv               # Resultado do Model Confidence Set (M0, M1, M2)
    ├── delta_svi_garchx_vs_egarchx.csv # p-valor do delta (ASVI): GARCH-X (M4) vs eGARCH-X (M2)
    ├── diagnosticos_residuos.csv       # Ljung-Box, Jarque-Bera dos resíduos padronizados
    ├── residuos_padronizados.csv       # Série de resíduos padronizados por ticker/modelo
    ├── diagnostico_filtros.txt         # Comparação das estratégias de filtro + estimativas M1-M4
    └── checkpoints_v2/                 # 1 .rds por ticker (retomada incremental, 4 modelos)
```

---

## Como executar

### 1. Instalar dependências

```bash
pip install -r requirements.txt
```

### 2. Rodar o pipeline em ordem

**Coleta e análise (Python):**

```bash
jupyter notebook
```

| Etapa | Arquivo | O que faz | Tempo estimado |
|-------|---------|-----------|----------------|
| 1 | `01_coleta_dados.ipynb` | Coleta preços (Yahoo Finance), calcula RV, coleta SVI (Google Trends), testa ADF, salva CSVs | ~1h (rate limit do Google Trends) |
| 3 | `03_analise_resultados.ipynb` | QLIKE, MCS (3 modelos), diagnósticos de resíduos, análise exploratória | ~5min (MCS) |

**Modelagem (R):**

```bash
Rscript 02_modelos.R
```

| Etapa | Arquivo | O que faz | Tempo estimado |
|-------|---------|-----------|----------------|
| 2 | `02_modelos.R` | Pré-voo do painel; estima M0-M4 em janela móvel; salva previsões, parâmetros e diagnósticos | ~30-60min (paralelo) |

Pacotes R necessários: `rugarch`, `dplyr`, `tidyr`, `readr`, `lubridate`, `tseries`.

> **Atalho:** com `data/base_final_tcc.csv` já salvo, pule o passo 1 e execute `02_modelos.R` seguido do notebook de análise.

---

## Dados

| Fonte | Biblioteca | Detalhes |
|-------|-----------|----------|
| Yahoo Finance | `yfinance` | Preços diários de fechamento ajustado; tickers B3 com sufixo `.SA` |
| Google Trends | `pytrends` | SVI semanal, geo=BR, termo de busca = ticker sem `.SA`; pausa de 12s entre requisições |

**Filtros de elegibilidade aplicados na coleta:**
- Volume financeiro médio diário >= R$ 1 milhão
- Pelo menos 5 anos de histórico de preços disponíveis
- Série de SVI coletada com sucesso

---

## Dependências

**Python** (`requirements.txt`):

```
arch          # Model Confidence Set (avaliação)
statsmodels   # teste ADF de estacionariedade
pytrends      # coleta do Google Trends
yfinance      # dados financeiros do Yahoo Finance
pandas
numpy
matplotlib
seaborn
tqdm
```

**R** (`02_modelos.R`):

```
rugarch       # estimação GARCH/GARCH-X e eGARCH/eGARCH-X (1,1)
dplyr
tidyr
readr
lubridate
tseries       # teste Jarque-Bera nos resíduos padronizados
parallel      # estimação paralela por ticker
```

---

## Referências

- Hansen, P. R., Lunde, A., & Nason, J. M. (2011). The model confidence set. *Econometrica*, 79(2), 453-497.
- Patton, A. J. (2011). Volatility forecast comparison using imperfect volatility proxies. *Journal of Econometrics*, 160(1), 246-256.
- Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity. *Journal of Econometrics*, 31(3), 307-327.
- Nelson, D. B. (1991). Conditional heteroskedasticity in asset returns: A new approach. *Econometrica*, 59(2), 347-370.
- DA, ZHI; ENGELBERG, JOSEPH; GAO, PENGJIE. In Search of Attention. The Journal of Finance, v. 66, n. 5, p. 1461–1499, 2011. DOI: https://doi.org/10.1111/j.1540-6261.2011.01679.x.