# sensitivity.m — Otimização de sensibilidade angular (MATLAB + COMSOL)

Busca hierárquica 4D que encontra, via modelo COMSOL, a geometria que maximiza a sensibilidade angular

```
S = d(alpha_peak)/dn   [deg/RIU]
```

onde `alpha_peak` é o ângulo de incidência onde `|TMOKE(alpha)|` é máximo, e `n` é o índice de refração do meio externo. O script é single-objective: não há otimização ou seleção por `|TMOKE|`, apenas por `|S|`. O cálculo de TMOKE(alpha) é feito internamente só porque `alpha_peak` é definido como o pico dessa curva — é o meio de localizar a ressonância rastreada, não um objetivo concorrente.

## Repositório

```
ic/
├── README.md
└── sensitivity.m
```

Os demais arquivos do projeto (outros scripts MATLAB, docs, manuscrito) existem localmente mas estão fora do controle de versão (`.gitignore`: `docs/`, `manuscript/`, `src/`).

## Pré-requisitos

- MATLAB com o LiveLink for MATLAB do COMSOL configurado (`com.comsol.model.*`).
- Um modelo `.mph` com:
  - parâmetros `h_au`, `L_domain`, `l_dente`, `h_si`, `n`;
  - um sweep paramétrico no `std1` com variáveis `alpha` [deg] e `m` (magnetização, `1`/`-1`);
  - Derived Values numéricos redirecionáveis para tabelas `tblTplus` / `tblTminus`, contendo colunas de ângulo e transmissão total.
- Ajustar `projectRootDir` (linha ~55) e `comsolModelFile` para o caminho local do projeto/`.mph`.

## Como a sensibilidade é calculada

Para um conjunto de índices de refração `n` (ex.: `[1.30, 1.33, 1.36]`), com `n = 1.33` como referência:

1. Resolve-se a curva de referência (`n = 1.33`) e localiza-se o pico global de `|TMOKE(alpha)|` → define `alpha_ref` e o sinal do pico.
2. Para cada `n` não-referência, resolve-se `alpha` apenas dentro da janela `alpha_ref ± trackingHalfWindowDeg` (20°, recortada aos limites globais de alpha) e escolhe-se o pico dentro dessa janela — preferindo o mesmo sinal do pico de referência.
3. Ajusta-se uma reta `alpha_peak(n)` por regressão linear (`polyfit`); o coeficiente angular é `S`.

Esse "rastreamento" evita que a métrica salte entre ressonâncias diferentes quando `n` muda — o que aconteceria se cada `n` fosse resolvido independentemente e o pico global de cada curva fosse usado.

TMOKE em si é calculado como:

```
TMOKE(alpha) = 2 * (T+(alpha) - T-(alpha)) / (T+(alpha) + T-(alpha))
```

com `T+`/`T-` sendo a transmissão total para `m = +1` e `m = -1`, lidas das tabelas COMSOL `tblTplus`/`tblTminus`.

## Geometria (busca 4D)

| Parâmetro | Significado | Grade COARSE |
|---|---|---|
| `L_domain` | período do domínio [nm] | `800:50:850` (2 pts) |
| `l_dente` | largura do dente [nm] | `500:50:600` (3 pts) |
| `h_si` | altura do silício [nm] | `[220, 240, 260]` (3 pts) |
| `h_au` | altura do ouro [nm] | `20:10:60` (5 pts) |

Total COARSE: `2×3×3×5 = 90` pontos de geometria. Cada ponto custa `2 × numel(fastRefractiveIndexSamples)` runs COMSOL (2 magnetizações × 3 índices de refração = 6 runs/ponto), logo `540` runs só no COARSE.

## Índices de refração

- `fastRefractiveIndexSamples = [1.30, 1.33, 1.36]` — usados em COARSE/FINE/SUPER.
- `trackingReferenceRefractiveIndex = 1.33` — curva de referência para travar a ressonância.
- `trackingHalfWindowDeg = 20` — meia-largura da janela de busca ao redor do pico de referência.
- `baselineRefractiveIndex = 1.33` — usado no snapshot final.
- `validationRefractiveIndexList = [1.30, 1.33, 1.36]` — usado no estágio SENSITIVITY FULL (mesma lista, alpha bem mais denso).

## Pipeline hierárquico

```
COARSE → seleciona TOP-K (topKCoarse=1) por |S|
  → FINE → seleciona TOP-K (topKFine=1) por |S|
    → SUPER → bestSensitivityCandidate (maior |S|)
      → SENSITIVITY FULL (curvas densas, sensitivityDense)
        → Snapshot .mph (opcional)
          → Export CSV/XLSX + figuras
```

| Estágio | Passo angular | Janela geométrica ao redor da semente | Janela angular |
|---|---|---|---|
| COARSE | 1.0° (`0`–`89`) | grade completa acima | `0`–`89` |
| FINE | 0.1° | `±10 nm` (L_domain/l_dente), `±5 nm` (h_si), `±2 nm` (h_au) | `alpha_peak_base ± 6°` |
| SUPER | 0.01° | `±4 nm` (L_domain/l_dente), `±2 nm` (h_si), `±1 nm` (h_au) | `alpha_peak_base ± 4°` |
| SENSITIVITY FULL | 0.01° | geometria fixa (`bestSensitivityCandidate`) | `0`–`89` completo, todos os `n` |

Seleção de sementes em cada estágio: `selectTopK_single_abs` — ordena por `abs(S_est_deg_per_RIU)` descendente e pega os `K` primeiros. Não existe score combinado nem seleção por TMOKE.

## Checkpoint e retomada

- `checkpointFilePath`: `<projectRootDir>/checkpoints/sensitivity_4d_checkpoint.mat`
- `progressWorkbookPath`: `<projectRootDir>/checkpoints/sensitivity_4d_progress.xlsx`
- `checkpointSchemaTag = 'tracked_same_peak_sensitivity_v1'` — checkpoints de outra versão/métrica são descartados (com warning) em vez de reutilizados.
- Salva a cada `checkpointEveryPoints = 10` pontos avaliados.
- Ao reiniciar, o script detecta o `stage` salvo (`COARSE`/`FINE`/`SUPER`/`SENSITIVITY FULL`/`FINAL`) e pula estágios já concluídos, restaurando as tabelas de resultado e sementes correspondentes.
- Ao final de uma execução completa, o arquivo de checkpoint é apagado.
- `MAX_RUNS = 20000`: o script aborta (`error`) se o total planejado de runs (exato a partir do planejamento do FINE) ultrapassar esse teto.

## Saídas

- **CSV** (em `projectRootDir`): `sensitivity_4d_coarse.csv`, `sensitivity_4d_fine.csv`, `sensitivity_4d_super.csv`, `sensitivity_4d_bestSensitivity_dense_ALLn.csv`.
- **XLSX de progresso**: abas `coarse`, `fine`, `super`, `dense`.
- **Snapshot `.mph`** (se `SAVE_SNAPSHOT = true`): geometria de `bestSensitivityCandidate` em `n = baselineRefractiveIndex`, `alpha = alphaBestDeg`.
- **Figuras** (se `MAKE_PLOTS`/`SAVE_FIGS = true`) em `<projectRootDir>/plots/sensitivity_4d_<timestamp>/`:
  - `phase_outputs/`: métrica `|S|` por candidato em cada estágio, curvas densas do SENSITIVITY FULL.
  - `iteration_tmoke/`: um PNG por ponto avaliado, se `PLOT_LIVE` e `SAVE_ITER_PLOTS` estiverem ativos.
  - Figuras finais: `TMOKE(alpha)` por `n`, `alpha_peak vs n` com ajuste linear, `|TMOKE|` no pico rastreado vs `n`.

## Funções principais

- `evaluateTrackedSensitivityAndCurves(...)` — orquestra o rastreamento de ressonância descrito acima; retorna `alphaPeakDegreesByN`, `trackedTmokeAbsByN`, `linearFit` e `sensitivitySlope` (= `S`).
- `solveAndGetTplusTminus(...)` — roda o COMSOL para `m=+1` e `m=-1` no sweep de `alpha` pedido, lê `tblTplus`/`tblTminus` e calcula `TMOKE(alpha)`.
- `selectTopK_single_abs(T, K, col)` — único critério de seleção de sementes/candidato final, por `abs(T.(col))` descendente.
- `save_checkpoint` / `maybe_checkpoint` / `write_progress_xlsx` — persistência de estado e progresso.
