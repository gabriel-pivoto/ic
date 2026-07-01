# Fluxograma textual detalhado do script MATLAB + COMSOL (5D Hierarchical Search)

## Visão geral

Este script executa uma **busca hierárquica em 5 dimensões geométricas** para encontrar a geometria que melhor equilibra:

* **máximo valor de |TMOKE|**
* **máxima sensibilidade angular**
  `S = d(alpha_peak)/dn [deg/RIU]`

A lógica é dividida em estágios:

1. **COARSE**: busca ampla e barata
2. **FINE**: refinamento local ao redor das melhores sementes
3. **SUPER**: refinamento final mais fino
4. **VALID**: validação densa com múltiplos índices de refração
5. **FULL**: curva final completa na condição baseline
6. **EXPORT / PLOTS**: salvar CSV, XLSX, snapshot e figuras

---

## Observação importante sobre a grade COARSE

No comentário do código aparece:

* `COARSE points = 270`

Mas, **pelo código real**, o total é **180 pontos**, porque:

* `domainPeriodGridNm = 800:50:850` → 2 valores (`800`, `850`)
* `toothWidthGridNm = 500:50:600` → 3 valores
* `siliconHeightGridNm = [220, 240, 260]` → 3 valores
* `ceyigHeightGridNm = [100, 140]` → 2 valores
* `goldHeightGridNm = 20:10:60` → 5 valores

Cálculo:

* `2 × 3 × 3 × 2 × 5 = 180`

Logo:

* **COARSE real = 180 pontos**
* cada ponto custa **4 runs**
* então **COARSE = 720 runs**

---

# 1. Fluxo principal do script

## 1.1 Início

**INÍCIO**
→ limpar ambiente MATLAB

* `clear`
* `clc`
* `close all`
* `format long`
* `tic`

→ importar bibliotecas do COMSOL

* `import com.comsol.model.*`
* `import com.comsol.model.util.*`

---

## 1.2 Configuração inicial

→ definir caminhos do projeto

* `projectRootDir`
* `comsolModelFile`

→ adicionar subpastas ao path

* `addpath(genpath(projectRootDir))`

→ definir orçamento máximo de runs

* `MAX_RUNS = 20000`

→ configurar checkpoint

* pasta de checkpoint
* arquivo `.mat` de checkpoint
* arquivo `.xlsx` de progresso
* salvar checkpoint a cada `10` pontos

---

## 1.3 Lógica de resume (retomada)

→ inicializar:

* `resumeFromCheckpoint = false`
* `checkpointData = struct`
* `resumeStageTag = ""`

→ verificar se existe arquivo de checkpoint

### Se existir

→ tentar carregar `checkpointData`

* se carregar:

  * `resumeFromCheckpoint = true`
  * recuperar:

    * estágio salvo
    * número de pontos já feitos
    * número de runs já executados
* se falhar:

  * emitir warning
  * começar do zero

### Se não existir

→ seguir como execução nova

→ abrir `diary(...)` para log em arquivo

---

## 1.4 Definição de parâmetros e tags

→ definir tag do estudo COMSOL

* `STUDY_TAG = 'std1'`

→ definir nomes dos parâmetros geométricos

* `h_au`
* `h_ceyig`
* `L_domain`
* `l_dente`
* `h_si`
* `n`

→ definir parâmetros de sweep

* `alpha`
* `m`

→ definir estados magnéticos

* `m = +1`
* `m = -1`

→ definir tags das tabelas de reflectância

* `tblRplus`
* `tblRminus`

---

## 1.5 Índices de refração usados

### Busca rápida (COARSE/FINE/SUPER)

* `fastRefractiveIndexSamples = [1.33, 1.39]`

### Índice baseline para TMOKE

* `baselineRefractiveIndex = 1.33`

### Validação final densa

* `validationRefractiveIndexList = [1.33, 1.36, 1.39]`

---

## 1.6 Grades geométricas COARSE

→ definir a grade 5D inicial

* `L_domain`: `800:50:850`
* `l_dente`: `500:50:600`
* `h_si`: `[220, 240, 260]`
* `h_ceyig`: `[100, 140]`
* `h_au`: `20:10:60`

---

## 1.7 Configuração angular

### COARSE

* `alphaCoarseRange = [0, 1.0, 89]`

### Refinamentos

* `alphaFineStep = 0.1`
* `alphaSuperStep = 0.01`
* `alphaDenseStep = 0.01`
* `alphaFullStep = 0.01`

---

## 1.8 Estratégia de promoção (TOP-K)

* `topKCoarse = 1`
* `topKFine = 1`

Isso significa:

* do COARSE sai **1** semente para o FINE
* do FINE sai **1** semente para o SUPER

---

## 1.9 Refinamentos geométricos

### FINE

* ouro: `±2`, passo `1`
* ceyig: `±5`, passo `1`
* domínio: `±10`, passo `5`
* dente: `±10`, passo `5`
* silício: `±5`, passo `5`

### SUPER

* ouro: `±1`, passo `0.5`
* ceyig: `±2`, passo `1`
* domínio: `±4`, passo `2`
* dente: `±4`, passo `2`
* silício: `±2`, passo `1`

---

## 1.10 Janelas angulares locais

### Janela nominal

* FINE: `±5°`
* SUPER: `±2°`

### Janela usada para sensibilidade (mais larga)

* FINE: `±6°`
* SUPER: `±4°`

Objetivo:

* evitar perder o pico quando o índice `n` muda

---

## 1.11 Configuração de outputs

* `SAVE_SNAPSHOT = true`
* `MAKE_PLOTS = true`
* `PLOT_LIVE = true`

→ salvar figuras:

* `SAVE_FIGS = true`
* formatos: `png` e `pdf`

→ criar pasta de plots com timestamp

---

## 1.12 Carregar modelo COMSOL

→ limpar modelos carregados

* `ModelUtil.clear`

→ carregar `.mph`

* `model = mphload(comsolModelFile)`

→ habilitar progresso do COMSOL no console

* `ModelUtil.showProgress(true)`

---

# 2. Planejamento global de runs

## 2.1 Inicialização global

→ iniciar cronômetro geral

* `globalTimerStart = tic`

→ inicializar:

* `runsCompletedGlobal = 0`
* `isTotalRunEstimateExact = false`

---

## 2.2 Custo por ponto

Cada ponto usa **2 índices de refração** (`1.33` e `1.39`).

Para cada índice:

* 1 run com `m = +1`
* 1 run com `m = -1`

Logo:

* `2 runs por n`
* `2 n`
* **4 runs por ponto**

---

## 2.3 Custo COARSE

* `coarseTotalPoints = 180`
* `coarseTotalRuns = 180 × 4 = 720`

---

## 2.4 Estimativa inicial do FINE

→ calcular número estimado de pontos ao redor de cada semente
→ multiplicar pelo `topKCoarse`
→ multiplicar por `4 runs por ponto`

Resultado:

* `fineTotalRunsEstimate`

---

## 2.5 Estimativa inicial do SUPER

→ calcular quantos pontos por semente existem no refinamento SUPER
→ multiplicar por `topKFine`
→ multiplicar por `4 runs por ponto`

Resultado:

* `superTotalRuns`

---

## 2.6 Runs extras fixos

Incluem:

* VALID: `2 × número de n de validação`
* SNAPSHOT: `2 runs` (se habilitado)
* FULL: `2 runs`

---

## 2.7 Total global inicial

→ somar:

* COARSE
* FINE estimado
* SUPER
* extras

→ imprimir resumo inicial no console

---

# 3. Estágio COARSE

## 3.1 Decisão: executar ou restaurar

### Se estiver retomando de um estágio posterior

Se `resumeStageTag` for:

* `FINE`
* `SUPER`
* `VALID`
* `FULL`

então:

→ restaurar do checkpoint:

* `coarseResultsTable`
* `coarseSeedCandidates`

→ **pular COARSE**

### Caso contrário

→ executar COARSE normalmente

---

## 3.2 Preparação do COARSE

→ criar acumuladores:

* `coarseRows = []`
* `coarsePointIndex = 0`

→ guardar referência para ETA:

* `stageRunsStart`
* `stageTotalRuns`
* `stageTimerStart`

### Se estiver retomando exatamente no COARSE

→ restaurar:

* `runsCompletedGlobal`
* `coarseRows`
* `coarsePointIndex`

---

## 3.3 Laço 5D completo do COARSE

Estrutura:

* para cada `L_domain`

  * setar no COMSOL
* para cada `l_dente`

  * setar no COMSOL
* para cada `h_si`

  * setar no COMSOL
* para cada `h_ceyig`

  * setar no COMSOL
* para cada `h_au`

  * setar no COMSOL

Para cada combinação geométrica:

### 3.3.1 Controle de retomada

→ incrementar `coarsePointIndex`

→ se esse ponto já foi feito em execução anterior:

* `continue`

---

### 3.3.2 Resolver para `n = 1.33`

→ setar `n = 1.33`

→ chamar `solveAndGetRplusRminus(...)` com sweep:

* `alpha = 0:1:89`

A função retorna:

* grade angular
* `Rplus`
* `Rminus`
* curva TMOKE

→ calcular:

* pico de `|TMOKE|`
* índice do pico
* ângulo do pico
* valor de TMOKE no pico

→ se o pico estiver na borda da faixa angular:

* emitir warning

→ se `PLOT_LIVE = true`:

* atualizar figura ao vivo

---

### 3.3.3 Resolver para `n = 1.39`

→ setar `n = 1.39`

→ repetir o solver no mesmo sweep angular

→ encontrar o ângulo do pico em `n2`

---

### 3.3.4 Calcular sensibilidade rápida

→ usar diferença finita entre os dois picos:

`sensitivityEstimateFast = (alpha_peak_n2 - alpha_peak_n1) / (n2 - n1)`

---

### 3.3.5 Armazenar o ponto

→ adicionar linha em `coarseRows` com:

* geometria
* `maxAbsTMOKE_base`
* `alpha_peak_base_deg`
* `TMOKE_at_peak_base`
* `alpha_peak_n2_deg`
* `S_est_deg_per_RIU`

---

### 3.3.6 Atualizar métricas de progresso

→ somar `4` em `runsCompletedGlobal`

→ calcular:

* fração concluída do estágio
* ETA do estágio
* fração global
* ETA global

→ escrever linha detalhada de log

---

### 3.3.7 Checkpoint periódico

→ incrementar `pointsSinceCheckpoint`

→ chamar `maybe_checkpoint(...)`

Se o número de pontos desde o último checkpoint atingir `10`:

* salvar checkpoint parcial
* tentar atualizar Excel
* resetar contador

---

## 3.4 Final do COARSE

→ converter `coarseRows` para `coarseResultsTable`

→ selecionar TOP-K por trade-off:

* rank por `|TMOKE|`
* rank por `|S|`
* somar ranks
* menor soma = melhor equilíbrio

→ guardar em `coarseSeedCandidates`

→ salvar checkpoint promovendo para `FINE`

→ escrever aba `coarse` no Excel

---

# 4. Planejamento exato do FINE

## 4.1 Garantir sementes válidas

Se `coarseSeedCandidates` estiver vazio, mas a execução retomou em `FINE`:

→ tentar restaurar do checkpoint
→ se não existir, abortar com erro

---

## 4.2 Calcular FINE exato

Para cada semente do COARSE:

→ criar listas refinadas com **clamping** aos limites das grades originais

* ouro
* ceyig
* domínio
* dente
* silício

→ multiplicar os tamanhos das listas
→ somar em `fineTotalPoints`

Depois:

* `fineTotalRuns = fineTotalPoints × 4`

---

## 4.3 Tornar o total global exato

→ atualizar:

* `globalRunTargetEstimate = coarse + fine(exato) + super + extras`

→ marcar:

* `isTotalRunEstimateExact = true`

→ imprimir:

* runs exatos do FINE
* total global exato

→ se passar de `MAX_RUNS`:

* abortar com erro

---

# 5. Estágio FINE

## 5.1 Decisão: executar ou restaurar

Se `resumeStageTag` for:

* `SUPER`
* `VALID`
* `FULL`

então:

→ restaurar:

* `fineResultsTable`
* `superSeedCandidates`

→ **pular FINE**

Caso contrário:

→ executar FINE

---

## 5.2 Preparação do FINE

→ criar:

* `fineRows = []`
* `finePointIndex = 0`

→ preparar ETA do estágio

### Se retomando exatamente em `FINE`

→ restaurar:

* `runsCompletedGlobal`
* `finePointIndex`
* `fineRows`
* `coarseSeedCandidates`

---

## 5.3 Laço do FINE

Para cada semente do COARSE:

### 5.3.1 Definir centro angular

→ usar o pico do COARSE em baseline:

* `alphaWindowCenterDeg = coarseSeedCandidates.alpha_peak_base_deg(s)`

---

### 5.3.2 Construir listas geométricas refinadas

→ gerar listas ao redor da semente, com clamping

---

### 5.3.3 Construir faixa angular FINE

→ usar:

* `alphaStartDeg = max(0, center - 6)`
* `alphaStopDeg = min(89, center + 6)`

---

### 5.3.4 Laço 5D refinado

Para cada combinação refinada:

→ setar geometria no COMSOL
→ incrementar `finePointIndex`

Se o ponto já tiver sido feito em execução anterior:

* `continue`

---

### 5.3.5 Resolver para `n = 1.33`

→ chamar `solveAndGetRplusRminus(...)` com:

* `alpha = alphaStart : 0.1 : alphaStop`

→ calcular pico de `|TMOKE|` e ângulo do pico

→ se `PLOT_LIVE = true`:

* atualizar figura

---

### 5.3.6 Resolver para `n = 1.39`

→ repetir na mesma janela angular

→ obter o ângulo do pico em `n2`

---

### 5.3.7 Calcular sensibilidade rápida

→ mesma fórmula de diferença finita

---

### 5.3.8 Armazenar e atualizar

→ adicionar linha em `fineRows`

→ somar `4` em `runsCompletedGlobal`

→ atualizar ETAs

→ escrever log

→ checkpoint periódico com `maybe_checkpoint(...)`

---

## 5.4 Final do FINE

→ converter `fineRows` para `fineResultsTable`

→ selecionar melhor(s) semente(s) por trade-off

→ salvar em `superSeedCandidates`

→ salvar checkpoint promovendo para `SUPER`

→ escrever aba `fine` no Excel

---

# 6. Planejamento exato do SUPER

## 6.1 Garantir sementes do FINE

Se `superSeedCandidates` estiver vazio e o resume estiver em `SUPER`:

→ tentar restaurar
→ se não existir, abortar com erro

---

## 6.2 Calcular SUPER exato

Para cada semente:

→ criar listas SUPER ao redor da semente

**Observação:** aqui o código não faz clamping explícito aos limites da grade COARSE.

→ multiplicar os tamanhos das listas
→ somar em `superTotalPoints`

Depois:

* `superTotalRuns = superTotalPoints × 4`

→ imprimir runs exatos do SUPER

---

# 7. Estágio SUPER

## 7.1 Decisão: executar ou restaurar

Se `resumeStageTag` for:

* `VALID`
* `FULL`

então:

→ restaurar:

* `superResultsTable`
* `bestTradeoffCandidate`
* `bestTmokeCandidate`
* `bestSensitivityCandidate`

→ **pular SUPER**

Caso contrário:

→ executar SUPER

---

## 7.2 Preparação do SUPER

→ criar:

* `superRows = []`
* `superPointIndex = 0`

→ preparar ETA

### Se retomando exatamente em `SUPER`

→ restaurar:

* `runsCompletedGlobal`
* `superPointIndex`
* `superRows`
* `superSeedCandidates`

---

## 7.3 Laço do SUPER

Para cada semente do FINE:

### 7.3.1 Definir centro angular

→ `alphaWindowCenterDeg = superSeedCandidates.alpha_peak_base_deg(s)`

---

### 7.3.2 Construir listas geométricas SUPER

→ listas mais finas em torno da semente

---

### 7.3.3 Construir faixa angular SUPER

→ usar:

* `alphaStartDeg = max(0, center - 4)`
* `alphaStopDeg = min(89, center + 4)`

---

### 7.3.4 Laço 5D superfino

Para cada combinação:

→ setar geometria
→ incrementar `superPointIndex`

Se já foi feito:

* `continue`

---

### 7.3.5 Resolver para `n = 1.33`

→ chamar solver com:

* `alpha = alphaStart : 0.01 : alphaStop`

→ calcular pico e ângulo do pico

→ se `PLOT_LIVE = true`:

* atualizar gráfico

---

### 7.3.6 Resolver para `n = 1.39`

→ repetir na mesma janela angular

→ obter ângulo do pico em `n2`

---

### 7.3.7 Calcular sensibilidade rápida

→ mesma diferença finita entre picos

---

### 7.3.8 Armazenar e atualizar

→ adicionar linha em `superRows`

→ somar `4` em `runsCompletedGlobal`

→ atualizar ETAs

→ logar

→ checkpoint periódico

---

## 7.4 Final do SUPER

→ converter `superRows` em `superResultsTable`

→ selecionar três resultados finais:

### Melhor trade-off

* `bestTradeoffCandidate`

### Melhor TMOKE puro

* `bestTmokeCandidate`

### Melhor sensibilidade pura

* `bestSensitivityCandidate`

→ imprimir os três no console

→ salvar checkpoint promovendo para `VALID`

→ escrever aba `super` no Excel

---

# 8. Estágio VALID (validação densa)

## 8.1 Garantir candidato final de trade-off

Se `bestTradeoffCandidate` estiver vazio e a execução tiver retomado em `VALID` ou `FULL`:

→ restaurar do checkpoint

---

## 8.2 Fixar geometria final

→ extrair da melhor linha:

* `Ldom_best`
* `Lden_best`
* `hsi_best`
* `hcey_best`
* `hau_best`

→ setar todos esses parâmetros no COMSOL

---

## 8.3 Preparar armazenamento denso

→ criar vetores e células para:

* `alphaPeakDegreesByN`
* `tmokeMaxAbsByN`
* curvas TMOKE
* grades angulares
* `Rplus`
* `Rminus`

---

## 8.4 Laço nos `n` de validação

Para cada `n` em `[1.33, 1.36, 1.39]`:

→ setar `n`

→ chamar `solveAndGetRplusRminus(...)` com sweep denso:

* `alpha = 0 : 0.01 : 89`

→ obter:

* `a`
* `Rp`
* `Rm`
* `TM`

→ calcular:

* `|TMOKE|_max`
* índice do pico
* `alpha_peak(n)`

→ armazenar:

* curva TMOKE
* reflectâncias
* grade angular

→ somar `2` em `runsCompletedGlobal`

---

## 8.5 Ajuste linear de sensibilidade

→ ajustar reta:

* `alpha_peak` em função de `n`

→ obter coeficiente angular:

* `sensitivityDense`

Esse é o valor final validado de:

* `d(alpha_peak)/dn`

---

## 8.6 Definir referência baseline

→ localizar o `n` da validação mais próximo do `baselineRefractiveIndex`

→ definir:

* `alphaBestDeg = alpha_peak` nesse `n`

→ usando a curva já armazenada da baseline:

* calcular `tmokeBestValue`

---

## 8.7 Construir tabela densa consolidada

→ iniciar `bestDenseTable = table()`

Para cada `n`:

→ criar uma tabela temporária com:

* geometria repetida
* `n`
* `alpha`
* `Rplus`
* `Rminus`
* `TMOKE`

→ concatenar tudo em `bestDenseTable`

---

## 8.8 Final do VALID

→ salvar checkpoint promovendo para `FULL`

→ escrever aba `dense` no Excel

---

# 9. Snapshot opcional

## 9.1 Verificar flag

Se `SAVE_SNAPSHOT = true`:

→ executar snapshot

Caso contrário:

→ pular

---

## 9.2 Gerar snapshot

→ setar `n = baselineRefractiveIndex`

→ configurar sweep com:

* `alpha = alphaBestDeg` (ponto único)
* `m = +1` e `m = -1`

→ rodar estudo COMSOL
→ atualizar derived values
→ somar `2` em `runsCompletedGlobal`

→ gerar nome do arquivo com:

* geometria
* `n`
* `alpha`
* timestamp

→ salvar snapshot `.mph`

→ registrar no log

---

# 10. Estágio FULL (curva baseline final)

## 10.1 Fixar baseline

→ setar `n = baselineRefractiveIndex`

---

## 10.2 Rodar sweep completo final

→ chamar `solveAndGetRplusRminus(...)` com:

* `alpha = 0 : 0.01 : 89`

→ obter:

* `alphaFull`
* `reflectancePlusFull`
* `reflectanceMinusFull`
* `tmokeFull`

→ somar `2` em `runsCompletedGlobal`

---

## 10.3 Plot ao vivo (opcional)

Se `PLOT_LIVE = true`:

→ atualizar figura com a curva FULL

---

## 10.4 Criar tabela final baseline

→ construir `bestFullTable` com:

* geometria final
* `n` baseline
* `alpha`
* `Rplus`
* `Rminus`
* `TMOKE`

---

# 11. Exportação de arquivos

## 11.1 Exportar tabelas dos estágios

Se existirem:

* `coarseResultsTable` → `tmoke_sens_5D_coarse.csv`
* `fineResultsTable` → `tmoke_sens_5D_fine.csv`
* `superResultsTable` → `tmoke_sens_5D_super.csv`

---

## 11.2 Exportar resultados finais

→ salvar:

* `tmoke_sens_bestTradeoff_dense_ALLn.csv`
* `tmoke_sens_bestTradeoff_full_baseline.csv`

→ atualizar aba `full` no Excel

---

## 11.3 Limpar checkpoint

Se o arquivo de checkpoint ainda existir:

→ deletar

---

# 12. Geração de plots finais

## 12.1 Verificar flag

Se `MAKE_PLOTS = true`:

→ gerar figuras finais

Caso contrário:

→ pular toda a seção

---

## 12.2 Figura 1 — TMOKE(alpha) para cada n

→ abrir figura
→ para cada `n` da validação:

* plotar `TMOKE(alpha)`

→ adicionar:

* labels
* título
* legenda

---

## 12.3 Figura 2 — alpha_peak vs n + ajuste linear

→ abrir figura
→ plotar pontos experimentais

→ gerar reta ajustada

→ adicionar:

* labels
* título
* legenda

---

## 12.4 Figura 3 — |TMOKE| máximo vs n

→ abrir figura
→ plotar:

* `validationRefractiveIndexList` vs `tmokeMaxAbsByN`

→ adicionar labels e título

---

## 12.5 Figura 4 — mapa de trade-off

→ concatenar todos os candidatos:

* COARSE
* FINE
* SUPER

→ abrir figura
→ plotar scatter:

* eixo x = `|TMOKE|_max`
* eixo y = `|S_est|`

→ destacar o `bestTradeoffCandidate`

---

# 13. Salvar todas as figuras

## 13.1 Verificar flag

Se `SAVE_FIGS = true`:

→ chamar `saveAllOpenFigures(...)`

→ exportar todas as figuras abertas em:

* `png`
* `pdf`

→ imprimir pasta de saída

Caso contrário:

→ não salvar figuras

---

# 14. Resumo final

→ calcular tempo total

→ imprimir no console:

* melhor trade-off
* melhor TMOKE puro
* melhor sensibilidade pura
* `sensitivityDense`
* número total de runs
* tempo total formatado

→ fechar `diary`

**FIM**

---

# 15. Subfluxo interno: `solveAndGetRplusRminus(...)`

Essa é a função central do script.

## Entrada

Recebe:

* modelo COMSOL
* tag do estudo
* nomes de `alpha` e `m`
* faixa angular
* strings de `m=+1` e `m=-1`
* tags das tabelas de saída

---

## Fluxo

### 15.1 Preparação

→ garantir existência das tabelas:

* `tblRplus`
* `tblRminus`

→ calcular quantos pontos angulares são esperados (`Npts`)

---

### 15.2 Resolver com `m = +1`

→ redirecionar todos os Derived Values para `tblRplus`
→ limpar `tblRplus`
→ configurar sweep paramétrico com:

* `alpha = range(...)`
* `m = +1`

→ rodar o estudo COMSOL
→ atualizar valores numéricos
→ ler a tabela
→ extrair:

* coluna de `alpha`
* coluna de reflectância

→ armazenar como:

* `alpha1_deg`
* `R1`

---

### 15.3 Resolver com `m = -1`

→ redirecionar numericals para `tblRminus`
→ limpar `tblRminus`
→ configurar sweep com:

* `alpha = range(...)`
* `m = -1`

→ rodar o estudo
→ atualizar numericals
→ ler a tabela
→ extrair:

* `alpha2_deg`
* `R2`

---

### 15.4 Validações

→ garantir que ambos os sweeps têm o número esperado de pontos
→ garantir que as grades de `alpha` são idênticas

---

### 15.5 Cálculo do TMOKE

→ definir:

* `Rplus = R1`
* `Rminus = R2`

→ calcular denominador:

* `denom = Rplus + Rminus`

→ evitar divisão por zero:

* onde `|denom| < 1e-9`, substituir por `1e-9`

→ calcular:

`TMOKE = (Rplus - Rminus) ./ denom`

---

### 15.6 Saída

Retornar:

* `alpha_deg`
* `Rplus`
* `Rminus`
* `TMOKE`

---

# 16. Subfluxo interno: `readAlphaAndRFromNamedTable(...)`

## Objetivo

Ler a tabela do COMSOL e identificar automaticamente:

* qual coluna é `alpha`
* qual coluna é reflectância

---

## Fluxo

→ chamar `mphtable(mdl, ttag)`

→ tentar obter cabeçalhos em diferentes campos:

* `colhead`
* `head`
* `header`
* `colnames`

### Se não houver cabeçalhos

→ assumir:

* coluna 1 = `alpha`
* coluna 2 = reflectância

### Se houver cabeçalhos

→ procurar uma coluna contendo `"alpha"`
→ procurar uma coluna contendo:

* `"reflectance"`
* `"total reflectance"`
* `"total r"`

### Se não encontrar nenhuma

→ usar fallback:

* coluna 1 = `alpha`
* coluna 2 = reflectância

---

## Validações finais

→ garantir pelo menos `Npts` linhas
→ usar apenas as últimas `Npts` linhas

### Se a unidade parecer estar em radianos

→ converter `alpha` para graus

---

## Saída

Retornar:

* `alpha_deg`
* `Rcol`

---

# 17. Subfluxo interno: `selectTopK_tradeoff(...)`

## Objetivo

Selecionar os melhores candidatos equilibrando:

* `|TMOKE|`
* `|S|`

---

## Fluxo

→ calcular:

* `tm = abs(TMOKE)`
* `ss = abs(S)`

→ ranquear `tm` em ordem decrescente
→ ranquear `ss` em ordem decrescente

→ calcular score:

* `score = rank_TM + rank_S`

→ menor score = melhor compromisso

→ ordenar por `score`

→ retornar as `K` primeiras linhas

---

# 18. Subfluxo interno: checkpoint

## 18.1 `save_checkpoint(...)`

### Objetivo

Salvar o estado parcial da execução.

### Fluxo

→ montar `checkpointData` com:

* estágio atual
* runs concluídos
* pontos concluídos
* payload parcial

→ salvar em arquivo temporário `.tmp`
→ mover para o arquivo final

Isso reduz risco de corromper o checkpoint.

---

## 18.2 `maybe_checkpoint(...)`

### Objetivo

Salvar automaticamente a cada N pontos.

### Fluxo

→ verificar se `pointsSinceCheckpoint >= checkpointEveryPoints`

### Se sim

→ chamar `save_checkpoint(...)`
→ tentar atualizar o workbook Excel
→ resetar contador para `0`

### Se não

→ não salvar nada

---

# 19. Subfluxo interno: `updateLivePlot(...)`

## Objetivo

Manter uma figura única e atualizada em tempo real.

---

## Fluxo

→ usar handles persistentes para:

* figura
* eixo
* curva TMOKE
* ponto de pico
* curva `R+`
* curva `R-`

### Se a figura ainda não existir

→ criar:

* figura
* eixo
* lado esquerdo para TMOKE
* lado direito para reflectância

→ criar linhas:

* `TMOKE(alpha)`
* pico TMOKE
* `R+`
* `R-`

---

## Atualização

→ atualizar dados das curvas
→ atualizar posição do pico
→ atualizar título com:

* estágio
* geometria
* valor de `|TMOKE|`
* ângulo do pico

→ `drawnow('limitrate')`

---

# 20. Subfluxo interno: `saveAllOpenFigures(...)`

## Objetivo

Salvar todas as figuras abertas.

---

## Fluxo

→ garantir que a pasta de saída exista
→ localizar todas as figuras abertas
→ ordenar por número da figura

Para cada figura:

→ obter `Name`
→ se não houver, usar `Tag`
→ se também não houver, usar `Figure_N`

→ sanitizar o nome do arquivo

Para cada formato solicitado:

* se `png` → `exportgraphics(..., Resolution=300)`
* se `pdf` → `exportgraphics(..., ContentType='vector')`
* senão → `saveas(...)`

→ ao final, imprimir quantas figuras foram salvas

---

# 21. Resumo conceitual em uma linha

O algoritmo executa:

**carregamento do modelo → estimativa de custo → varredura ampla (COARSE) → refinamento local (FINE) → refinamento superfino (SUPER) → escolha do melhor compromisso → validação densa em múltiplos índices (VALID) → curva final baseline (FULL) → exportação de dados, snapshot e figuras**

---

# 22. Fluxograma resumido em estilo seta

**INÍCIO**
→ limpar ambiente MATLAB
→ importar COMSOL
→ definir paths, parâmetros e flags
→ configurar checkpoint / resume
→ carregar modelo `.mph`
→ estimar runs e ETA
→ **COARSE**
→ selecionar melhor(es) por trade-off
→ recalcular orçamento exato do **FINE**
→ **FINE**
→ selecionar melhor(es) por trade-off
→ recalcular orçamento exato do **SUPER**
→ **SUPER**
→ selecionar:

* melhor trade-off
* melhor TMOKE
* melhor sensibilidade
  → **VALID** com `n = [1.33, 1.36, 1.39]`
  → ajustar `alpha_peak vs n`
  → obter `sensitivityDense`
  → **SNAPSHOT** (opcional)
  → **FULL** em `n = baseline`
  → exportar CSV / XLSX
  → gerar plots
  → salvar figuras
  → imprimir resumo final
  → **FIM**

