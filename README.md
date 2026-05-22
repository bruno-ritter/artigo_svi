# TCC — SVI e Volatilidade de Ações da B3

Trabalho de Conclusão de Curso em Economia.

**Questão de pesquisa:** O Search Volume Index (SVI) do Google Trends melhora a previsão de volatilidade de ações listadas na B3?

**Status:** pipeline completo (coleta, estimação e avaliação). 302 tickers estimados, 291 utilizados na comparação final (11 removidos por erro de matriz singular no M2). O MCS (alpha = 10%) inclui apenas o M1 (GARCH puro); o M2 (GARCH-X com SVI log-dev) foi excluído com p-valor = 0.074.

---

## Metodologia

### Variáveis

**Volatilidade Realizada (RV)**
Proxy semanal de volatilidade calculada como a soma dos quadrados dos retornos logarítmicos diários dentro de cada semana:

$$RV_t = \sum_{d \in t} r_d^2, \quad r_d = \ln\left(\frac{P_d}{P_{d-1}}\right)$$

**Search Volume Index (SVI)**
O SVI é o índice semanal normalizado (0-100) do Google Trends para o termo de busca correspondente ao ticker, com geo=BR. O tratamento do SVI que entra como regressor externo no modelo GARCH-X é o desvio logarítmico em relação à mediana móvel de 8 semanas passadas:

$$\text{SVI-log-dev}_t = \ln(1 + \text{SVI}_t) - \ln\left(1 + \underset{k \in \{t-8,\ldots,t-1\}}{\mathrm{median}}(\text{SVI}_k)\right)$$

O uso de $\ln(1+\cdot)$ (`log1p`) evita problemas quando SVI = 0. A série passa pelo teste de estacionariedade ADF na etapa de coleta.

O `svi_log_dev` é calculado conforme DA, ENGELBERG e GAO (2011).

---

### Modelos Competidores

| # | Modelo | Equação da variância condicional |
|---|--------|----------------------------------|
| M0 | Random Walk | $\hat{\sigma}^2_t = RV_{t-1}$ |
| M1 | GARCH(1,1) | $\sigma^2_t = \omega + \alpha\,\varepsilon^2_{t-1} + \beta\,\sigma^2_{t-1}$ |
| M2 | GARCH-X(1,1) | $\sigma^2_t = \omega + \alpha\,\varepsilon^2_{t-1} + \beta\,\sigma^2_{t-1} + \gamma\,\text{SVI-log-dev}_{t-1}$ |

M0 é o benchmark ingênuo. M1 é o modelo padrão de volatilidade condicional. M2 estende M1 com o SVI log-dev como regressor externo; o coeficiente $\gamma$ indica se o interesse de busca antecipa volatilidade futura.

---

### Estimação em Janela Móvel (*Rolling Window*)

Os modelos são estimados e avaliados **fora da amostra** (*out-of-sample*) por janela móvel de tamanho fixo:

- **Janela de treino:** 156 semanas (3 anos)
- **Horizonte de previsão:** 1 passo à frente (semana seguinte)
- **Procedimento:** a cada semana *t*, estima-se o modelo nas 156 observações anteriores e gera-se uma previsão para *t*; a janela avança uma semana e repete-se o processo
- **Primeira previsão:** indice t = 157 na série já filtrada (as 8 primeiras semanas com `svi_log_dev` = NaN por falta de histórico para a mediana são descartadas antes da estimação na coleta de dados)
- **Painel alinhado:** M0, M1 e M2 são estimados nas **mesmas** semanas calendário por ticker, bloco contíguo, sem pular semanas. Tickers que falham em qualquer semana do bloco (dados incompletos) são excluídos no pré-voo; os restantes têm o **mesmo** número de previsões

Isso evita *look-ahead bias*, reproduz o ambiente real de previsão e garante comparabilidade justa entre modelos e entre ativos.

A estimação GARCH/GARCH-X é feita em R (`rugarch`) via `02_modelos.R`.

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

Tickers com alta proporção de semanas com SVI = 0 indicam que o Google Trends não registrou buscas suficientes para o papel. Isso causa dois problemas: erro de matriz singular na estimação do GARCH-X e viés por ruído puro no regressor externo.

O repositório avalia **três estratégias de filtro** comparativamente (diagnóstico salvo em `data/diagnostico_filtros.txt`):

| Estratégia | Descrição |
|------------|-----------|
| Sem filtro de zeros | Não remove nenhum ticker |
| Zeros >= 70% | Remove tickers com 70% ou mais das semanas com SVI = 0, além dos com erro |
| Zeros >= 50% | Remove tickers com 50% ou mais das semanas com SVI = 0, além dos com erro |

A análise principal usa a estratégia **sem filtro de zeros** (apenas remoção de erros de singular matrix), totalizando 291 tickers.

---

## Resultados Preliminares

| Modelo | QLIKE médio | Melhor em N tickers |
|--------|-------------|---------------------|
| M0_RW | 63.107 | 2 |
| M1_GARCH | 0.928 | 160 |
| M2_GARCHX | 0.966 | 129 |

O M1 (GARCH puro) domina no agregado e é o único modelo no MCS. Apesar disso, o M2 vence individualmente em 129 dos 291 tickers, sugerindo que o SVI tem poder preditivo em um subconjunto de ativos. Os diagnósticos de resíduos padronizados (Ljung-Box em $z_t$ e $z_t^2$, Jarque-Bera) estão em `data/diagnosticos_residuos.csv`.

---

## Estrutura do repositório

```
.
├── 01_coleta_dados.ipynb        # Coleta e preparo dos dados
├── 02_modelos.R                 # Estimação M0-M2 (rugarch), painel alinhado
├── 03_analise_resultados.ipynb  # QLIKE, MCS, diagnósticos e análise exploratória
├── requirements.txt
└── data/
    ├── acoes-listadas-b3.csv        # Lista de tickers da B3
    ├── acoes_elegiveis.csv          # Tickers que passaram nos filtros de liquidez/histórico
    ├── precos_semanais.csv          # Preços e volatilidade realizada semanal por ticker
    ├── google_trends_svi.csv        # SVI semanal por ticker (Google Trends)
    ├── log_coleta_svi.csv           # Log da coleta do Google Trends
    ├── log_adf_svi.csv              # Resultados do teste ADF por ticker
    ├── base_final_tcc.csv           # Base combinada (preços + SVI + tratamentos)
    ├── previsoes_consolidadas.csv   # Previsões OOS M0-M2 (mesmas semanas por ticker)
    ├── resumo_parametros.csv        # Parâmetros médios por ticker (omega, alpha, beta, gamma)
    ├── apendice_resumo_params.csv   # Resumo arredondado para apêndice
    ├── qlike_por_ticker.csv         # QLIKE médio por ticker e modelo (291 tickers)
    ├── mcs_resultado.csv            # Resultado do Model Confidence Set (3 modelos)
    ├── diagnosticos_residuos.csv    # Ljung-Box, Jarque-Bera dos resíduos padronizados
    ├── residuos_padronizados.csv    # Série de resíduos padronizados por ticker/modelo
    └── diagnostico_filtros.txt      # Comparação das estratégias de filtro de zeros
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
| 2 | `02_modelos.R` | Pré-voo do painel; estima M0-M2 em janela móvel; salva previsões, parâmetros e diagnósticos | ~30-60min (paralelo) |

Pacotes R necessários: `rugarch`, `dplyr`, `readr`, `lubridate`.

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
rugarch       # estimação GARCH(1,1) e GARCH-X(1,1)
dplyr
readr
lubridate
parallel      # estimação paralela por ticker
```

---

## Referências

- Hansen, P. R., Lunde, A., & Nason, J. M. (2011). The model confidence set. *Econometrica*, 79(2), 453-497.
- Patton, A. J. (2011). Volatility forecast comparison using imperfect volatility proxies. *Journal of Econometrics*, 160(1), 246-256.
- Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity. *Journal of Econometrics*, 31(3), 307-327.
- DA, ZHI; ENGELBERG, JOSEPH; GAO, PENGJIE. In Search of Attention. The Journal of Finance, v. 66, n. 5, p. 1461–1499, 2011. DOI: https://doi.org/10.1111/j.1540-6261.2011.01679.x.