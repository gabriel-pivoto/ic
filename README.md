# Fluxograma detalhado do script MATLAB + COMSOL  
## 5D Hierarchical Search para otimização de TMOKE e sensibilidade angular

---

## 1. Objetivo do script

Este script executa uma **busca hierárquica em 5 dimensões geométricas** para encontrar a geometria que melhor equilibra dois objetivos:

1. **Maximizar o valor absoluto de TMOKE**
2. **Maximizar a sensibilidade angular**
   - `S = d(alpha_peak)/dn [deg/RIU]`

A busca é feita em múltiplos estágios, reduzindo custo computacional ao começar com uma varredura ampla e depois refinar localmente apenas os melhores candidatos.

---

## 2. Métricas otimizadas

### 2.1 Métrica de TMOKE

Para cada geometria, o script calcula:

- `maxAbsTMOKE_base = max(|TMOKE(alpha)|)` no índice de refração baseline

Isso representa a intensidade máxima do efeito magneto-óptico transversal para aquela geometria.

---

### 2.2 Métrica de sensibilidade rápida

Durante as etapas de busca (`COARSE`, `FINE`, `SUPER`), a sensibilidade é estimada rapidamente usando apenas **dois índices de refração**:

- `n1 = 1.33`
- `n2 = 1.39`

A estimativa é:

- `S_est = (alpha_peak(n2) - alpha_peak(n1)) / (n2 - n1)`

Essa é uma aproximação rápida da derivada angular por índice de refração.

---

### 2.3 Métrica de sensibilidade densa (validação final)

Na etapa `VALID`, a geometria final escolhida é reavaliada com uma lista mais densa de índices:

- `n = [1.33, 1.36, 1.39]`

O script então:

1. calcula `alpha_peak` para cada `n`
2. ajusta uma reta `alpha_peak vs n`
3. usa o coeficiente angular dessa reta como:

- `sensitivityDense`

Esse valor é a estimativa final validada de `d(alpha_peak)/dn`.

---

## 3. Variáveis geométricas da busca 5D

A geometria é varrida em 5 dimensões:

- `L_domain` → período do domínio
- `l_dente` → largura do dente
- `h_si` → altura do silício
- `h_ceyig` → altura da camada de Ce:YIG
- `h_au` → altura da camada de ouro

---

## 4. Estrutura global do algoritmo

O fluxo geral do algoritmo é:

1. carregar ambiente MATLAB e modelo COMSOL
2. configurar parâmetros, checkpoints e orçamento de runs
3. estimar custo computacional total
4. executar busca `COARSE`
5. promover os melhores candidatos para `FINE`
6. executar refinamento `FINE`
7. promover os melhores candidatos para `SUPER`
8. executar refinamento `SUPER`
9. escolher:
   - melhor trade-off
   - melhor TMOKE puro
   - melhor sensibilidade pura
10. validar o melhor trade-off com curvas densas (`VALID`)
11. gerar curva final baseline (`FULL`)
12. exportar resultados (CSV/XLSX)
13. gerar e salvar figuras
14. imprimir resumo final

---

## 5. Configuração inicial

### 5.1 Inicialização do MATLAB

O script começa limpando completamente o ambiente:

- limpa variáveis
- limpa console
- fecha figuras
- define exibição numérica longa
- inicia cronômetro global

---

### 5.2 Importação da API do COMSOL

São importadas as classes:

- `com.comsol.model.*`
- `com.comsol.model.util.*`

Essas bibliotecas permitem controlar o modelo `.mph` diretamente via MATLAB.

---

### 5.3 Configuração de caminhos

O script define:

- diretório raiz do projeto
- caminho do arquivo `.mph`
- inclusão de subpastas via `addpath(genpath(...))`

---

### 5.4 Orçamento máximo de execução

É definido um teto de segurança:

- `MAX_RUNS = 20000`

Se o número total planejado de runs ultrapassar esse limite, o script aborta.

---

## 6. Sistema de checkpoint e retomada

### 6.1 Objetivo

Como o processo pode ser longo, o script salva progresso periodicamente para permitir:

- retomada após interrupção
- preservação dos resultados parciais
- continuidade sem perder runs já executados

---

### 6.2 Arquivos usados

O sistema de checkpoint usa:

- um arquivo `.mat` com estado interno
- um arquivo `.xlsx` com tabelas de progresso legíveis

---

### 6.3 Comportamento de retomada

Na inicialização, o script verifica se existe um checkpoint.

### Se existir:

- tenta carregar o conteúdo salvo
- recupera:
  - estágio salvo
  - quantidade de pontos já processados
  - número total de runs já executados

### Se o checkpoint estiver corrompido:

- emite warning
- reinicia do zero

### Se não existir:

- inicia normalmente do começo

---

## 7. Índices de refração utilizados

### 7.1 Busca rápida

Durante `COARSE`, `FINE` e `SUPER`, o script usa:

- `fastRefractiveIndexSamples = [1.33, 1.39]`

Esses dois valores são suficientes para estimar a sensibilidade com baixo custo.

---

### 7.2 Baseline

O índice baseline é:

- `baselineRefractiveIndex = 1.33`

Ele é usado como referência principal para a métrica de TMOKE.

---

### 7.3 Validação final

Na etapa `VALID`, são usados:

- `validationRefractiveIndexList = [1.33, 1.36, 1.39]`

Esses três pontos permitem calcular uma inclinação linear mais confiável.

---

## 8. Grade geométrica COARSE

A varredura inicial usa as seguintes listas:

- `L_domain = 800:50:850`
- `l_dente = 500:50:600`
- `h_si = [220, 240, 260]`
- `h_ceyig = [100, 140]`
- `h_au = 20:10:60`

### Total real de pontos

Apesar do comentário do código dizer `270`, o total real é:

- `2 × 3 × 3 × 2 × 5 = 180 pontos`

Cada ponto custa:

- `4 runs`

Logo:

- `COARSE = 720 runs`

---

## 9. Estratégia hierárquica

### 9.1 Estágio COARSE

Objetivo:

- explorar a região ampla do espaço 5D
- medir TMOKE e sensibilidade rápida em toda a grade inicial

Características:

- sweep angular grosso: `0:1:89`
- custo menor por ponto
- alta cobertura, menor precisão angular

---

### 9.2 Estágio FINE

Objetivo:

- refinar a busca ao redor dos melhores candidatos do COARSE

Características:

- usa janelas geométricas menores
- usa sweep angular mais fino
- calcula novamente TMOKE e sensibilidade rápida

---

### 9.3 Estágio SUPER

Objetivo:

- refinamento final mais preciso

Características:

- passos geométricos ainda menores
- sweep angular superfino
- define os melhores candidatos finais

---

### 9.4 Estágio VALID

Objetivo:

- validar o melhor candidato de trade-off com curvas densas em múltiplos `n`

Características:

- sweep angular global denso
- ajuste linear de `alpha_peak vs n`
- cálculo final de `sensitivityDense`

---

### 9.5 Estágio FULL

Objetivo:

- gerar a curva final completa da melhor geometria na condição baseline

Características:

- sweep completo de `alpha`
- dados limpos para exportação e comparação

---

## 10. Seleção dos melhores candidatos

### 10.1 Critério de trade-off

Em cada etapa, os candidatos são ranqueados por dois critérios:

- rank em `|TMOKE|`
- rank em `|S|`

Depois o script calcula:

- `score_tradeoff = rank_TM + rank_S`

### Interpretação

- menor `score_tradeoff` = melhor compromisso entre TMOKE e sensibilidade

---

### 10.2 Seleções feitas ao final do SUPER

O script extrai três candidatos finais:

1. **bestTradeoffCandidate**
   - melhor equilíbrio entre TMOKE e sensibilidade

2. **bestTmokeCandidate**
   - melhor apenas em `|TMOKE|`

3. **bestSensitivityCandidate**
   - melhor apenas em `|S|`

---

## 11. Fluxo detalhado por estágio

### 11.1 COARSE

Para cada ponto da grade 5D:

1. setar parâmetros geométricos no COMSOL
2. resolver com `n = 1.33`
3. encontrar:
   - curva TMOKE
   - pico de `|TMOKE|`
   - ângulo do pico
4. resolver com `n = 1.39`
5. encontrar o novo ângulo de pico
6. calcular:
   - `S_est`
7. armazenar resultados
8. atualizar:
   - runs globais
   - ETA do estágio
   - ETA global
9. salvar checkpoint periodicamente

No final:

- converter resultados em tabela
- selecionar TOP-K por trade-off
- promover para `FINE`

---

### 11.2 FINE

Para cada semente do COARSE:

1. usar `alpha_peak_base` como centro angular
2. construir janelas geométricas refinadas
3. limitar `alpha` a uma faixa local ao redor do pico
4. repetir, para cada ponto:
   - resolver em `n = 1.33`
   - resolver em `n = 1.39`
   - recalcular `S_est`
   - salvar resultados
   - atualizar ETAs
   - salvar checkpoint

No final:

- formar `fineResultsTable`
- selecionar sementes para `SUPER`

---

### 11.3 SUPER

Para cada semente do FINE:

1. usar o pico anterior como centro angular
2. construir janelas geométricas ainda mais finas
3. usar step angular superfino
4. repetir:
   - cálculo em `n = 1.33`
   - cálculo em `n = 1.39`
   - cálculo de `S_est`
   - atualização de progresso
   - checkpoint

No final:

- formar `superResultsTable`
- selecionar:
  - melhor trade-off
  - melhor TMOKE
  - melhor sensibilidade

---

### 11.4 VALID

Com a geometria do `bestTradeoffCandidate` fixada:

1. setar a geometria ótima
2. para cada `n` em `[1.33, 1.36, 1.39]`:
   - varrer `alpha = 0:0.01:89`
   - calcular curva TMOKE
   - localizar `alpha_peak`
   - registrar `|TMOKE|_max`
3. ajustar `alpha_peak vs n`
4. extrair a inclinação:
   - `sensitivityDense`

No final:

- consolidar todas as curvas em `bestDenseTable`

---

### 11.5 FULL

Com a mesma geometria ótima:

1. fixar `n = baseline`
2. varrer `alpha = 0:0.01:89`
3. obter:
   - `Rplus`
   - `Rminus`
   - `TMOKE`
4. montar `bestFullTable`

---

## 12. Solver central: `solveAndGetRplusRminus(...)`

Essa função é o coração do script.

### Ela faz:

1. garantir existência das tabelas no COMSOL
2. configurar sweep em `alpha`
3. rodar com `m = +1`
4. ler tabela e extrair:
   - `alpha`
   - `Rplus`
5. configurar sweep em `alpha`
6. rodar com `m = -1`
7. ler tabela e extrair:
   - `alpha`
   - `Rminus`
8. verificar consistência entre as grades angulares
9. calcular:

- `TMOKE = (Rplus - Rminus) / (Rplus + Rminus)`

com proteção contra divisão por zero

### Saída da função

Retorna:

- `alpha_deg`
- `Rplus`
- `Rminus`
- `TMOKE`

---

## 13. Outputs gerados

Ao final, o script pode gerar:

### 13.1 CSVs

- `tmoke_sens_5D_coarse.csv`
- `tmoke_sens_5D_fine.csv`
- `tmoke_sens_5D_super.csv`
- `tmoke_sens_bestTradeoff_dense_ALLn.csv`
- `tmoke_sens_bestTradeoff_full_baseline.csv`

---

### 13.2 XLSX de progresso

Com abas como:

- `coarse`
- `fine`
- `super`
- `dense`
- `full`

---

### 13.3 Snapshot `.mph`

Se habilitado, salva um snapshot COMSOL da melhor geometria.

---

### 13.4 Figuras

Se habilitado, salva:

- PNG
- PDF

das figuras abertas ao final da execução.

---

## 14. Gráficos finais gerados

Se `MAKE_PLOTS = true`, o script gera:

1. **TMOKE(alpha) para cada n**
2. **alpha_peak vs n + ajuste linear**
3. **|TMOKE| máximo vs n**
4. **mapa de trade-off**
   - todos os candidatos
   - destaque do melhor trade-off

---

## 15. Custos computacionais e lógica de eficiência

O script foi estruturado para economizar tempo computacional:

- usa apenas **2 valores de `n`** nas etapas de busca
- faz sweep angular grosso no início
- refina apenas em torno de candidatos promissores
- usa checkpoint para evitar perder progresso
- só faz varredura densa no final, em uma única geometria

Isso evita varrer todo o espaço 5D com alta resolução desde o começo.

---

## 16. Pontos de atenção / gargalos

### 16.1 Tempo de execução

Os maiores custos estão em:

- chamadas repetidas ao COMSOL
- sweeps angulares densos
- etapas FINE e SUPER se as janelas forem grandes

---

### 16.2 Robustez

O script protege contra:

- perda de progresso com checkpoint
- erros de leitura de tabela com tentativas de fallback
- divisão por zero no cálculo de TMOKE

---

### 16.3 Possível inconsistência comentada

O comentário no código diz:

- `COARSE points = 270`

Mas a conta real da grade é:

- `180`

Essa diferença é documental, não algorítmica.

---

## 17. Resumo conceitual final

Em termos simples, o script faz:

**carregar o modelo → explorar o espaço geométrico de forma ampla → refinar localmente os melhores candidatos → escolher o melhor compromisso entre TMOKE e sensibilidade → validar com múltiplos índices de refração → exportar curvas, tabelas e figuras**

---

# Fluxograma ASCII (para documentação)

```text
+--------------------------------------------------------------+
|                           INÍCIO                             |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| Limpa ambiente MATLAB                                         |
| - clear / clc / close all                                     |
| - format long                                                 |
| - tic                                                         |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| Importa bibliotecas do COMSOL                                 |
| - com.comsol.model.*                                          |
| - com.comsol.model.util.*                                     |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| Configura projeto                                             |
| - paths                                                       |
| - arquivo .mph                                                |
| - MAX_RUNS                                                    |
| - pasta de checkpoint                                         |
| - workbook de progresso                                       |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| Verifica checkpoint existente?                                |
+--------------------------------------------------------------+
              | SIM                                   | NÃO
              v                                       v
+------------------------------+        +------------------------------+
| Carrega checkpoint           |        | Continua execução nova       |
| - stage                      |        | - resume = false             |
| - done_points                |        |                              |
| - runsCompletedGlobal        |        |                              |
+------------------------------+        +------------------------------+
              \___________________________   __________________________/
                                          \ /
                                           v
+--------------------------------------------------------------+
| Define parâmetros globais                                     |
| - nomes dos parâmetros geométricos                            |
| - alpha / m                                                   |
| - índices de refração                                         |
| - grids COARSE                                                |
| - deltas e steps FINE/SUPER                                   |
| - flags de plot / snapshot / save figs                        |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| Carrega modelo COMSOL (.mph)                                  |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| Planejamento global de runs                                   |
| - calcula runs por ponto (= 4)                                |
| - calcula custo COARSE                                        |
| - estima FINE / SUPER                                         |
| - calcula extras                                              |
| - imprime ETA inicial                                         |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| Executar COARSE ou restaurar?                                 |
+--------------------------------------------------------------+
              | RESTAURAR                               | EXECUTAR
              v                                         v
+------------------------------+        +------------------------------+
| Recupera coarseResultsTable  |        | Laço 5D COARSE              |
| Recupera coarseSeedCandidates|        | Para cada geometria:         |
+------------------------------+        | 1) solve n=1.33             |
                                        | 2) acha pico TMOKE          |
                                        | 3) solve n=1.39             |
                                        | 4) estima S_est             |
                                        | 5) salva linha              |
                                        | 6) atualiza ETA             |
                                        | 7) checkpoint periódico     |
                                        +------------------------------+
                                                      |
                                                      v
                                        +------------------------------+
                                        | Final do COARSE              |
                                        | - cria coarseResultsTable    |
                                        | - rank trade-off             |
                                        | - escolhe TOP-K              |
                                        | - promove para FINE          |
                                        +------------------------------+
                                                      |
                                                      v
+--------------------------------------------------------------+
| Planejamento exato do FINE                                    |
| - constrói janelas refinadas                                  |
| - calcula fineTotalRuns                                       |
| - torna ETA global exato                                      |
| - valida MAX_RUNS                                             |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| Executar FINE ou restaurar?                                   |
+--------------------------------------------------------------+
              | RESTAURAR                               | EXECUTAR
              v                                         v
+------------------------------+        +------------------------------+
| Recupera fineResultsTable    |        | Para cada semente COARSE:    |
| Recupera superSeedCandidates |        | 1) define centro angular     |
+------------------------------+        | 2) cria janelas refinadas    |
                                        | 3) laço 5D local             |
                                        | 4) solve n=1.33              |
                                        | 5) solve n=1.39              |
                                        | 6) calcula S_est             |
                                        | 7) salva linha               |
                                        | 8) ETA + checkpoint          |
                                        +------------------------------+
                                                      |
                                                      v
                                        +------------------------------+
                                        | Final do FINE                |
                                        | - cria fineResultsTable      |
                                        | - rank trade-off             |
                                        | - escolhe TOP-K              |
                                        | - promove para SUPER         |
                                        +------------------------------+
                                                      |
                                                      v
+--------------------------------------------------------------+
| Planejamento exato do SUPER                                   |
| - constrói janelas superfina                                  |
| - calcula superTotalRuns                                      |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| Executar SUPER ou restaurar?                                  |
+--------------------------------------------------------------+
              | RESTAURAR                               | EXECUTAR
              v                                         v
+------------------------------+        +------------------------------+
| Recupera superResultsTable   |        | Para cada semente FINE:      |
| Recupera melhores candidatos |        | 1) define centro angular     |
+------------------------------+        | 2) cria janelas superfina    |
                                        | 3) laço 5D local             |
                                        | 4) solve n=1.33              |
                                        | 5) solve n=1.39              |
                                        | 6) calcula S_est             |
                                        | 7) salva linha               |
                                        | 8) ETA + checkpoint          |
                                        +------------------------------+
                                                      |
                                                      v
                                        +------------------------------+
                                        | Final do SUPER               |
                                        | Seleciona:                   |
                                        | - bestTradeoffCandidate      |
                                        | - bestTmokeCandidate         |
                                        | - bestSensitivityCandidate   |
                                        | Promove para VALID           |
                                        +------------------------------+
                                                      |
                                                      v
+--------------------------------------------------------------+
| VALID (validação densa)                                       |
| - fixa bestTradeoffCandidate                                  |
| - para n = [1.33, 1.36, 1.39]:                               |
|   * sweep alpha denso                                         |
|   * calcula TMOKE(alpha)                                      |
|   * encontra alpha_peak                                       |
| - ajusta alpha_peak vs n                                      |
| - extrai sensitivityDense                                     |
| - monta bestDenseTable                                        |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| Snapshot habilitado?                                          |
+--------------------------------------------------------------+
              | SIM                                   | NÃO
              v                                       v
+------------------------------+        +------------------------------+
| Salva snapshot .mph          |        | Pula snapshot                |
| na geometria ótima           |        |                              |
+------------------------------+        +------------------------------+
              \___________________________   __________________________/
                                          \ /
                                           v
+--------------------------------------------------------------+
| FULL (curva final baseline)                                   |
| - fixa n = baseline                                           |
| - sweep alpha completo                                         |
| - obtém Rplus, Rminus, TMOKE                                  |
| - monta bestFullTable                                         |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| Exporta resultados                                            |
| - CSVs COARSE / FINE / SUPER                                  |
| - CSV denso final                                             |
| - CSV baseline final                                          |
| - atualiza XLSX                                               |
| - apaga checkpoint final                                      |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
| MAKE_PLOTS = true ?                                           |
+--------------------------------------------------------------+
              | SIM                                   | NÃO
              v                                       v
+------------------------------+        +------------------------------+
| Gera figuras finais          |        | Pula geração de plots        |
| 1) TMOKE(alpha) por n        |        |                              |
| 2) alpha_peak vs n           |        |                              |
| 3) |TMOKE|max vs n           |        |                              |
| 4) mapa de trade-off         |        |                              |
+------------------------------+        +------------------------------+
              \___________________________   __________________________/
                                          \ /
                                           v
+--------------------------------------------------------------+
| SAVE_FIGS = true ?                                            |
+--------------------------------------------------------------+
              | SIM                                   | NÃO
              v                                       v
+------------------------------+        +------------------------------+
| Salva todas as figuras       |        | Não salva figuras            |
| em PNG e PDF                 |        |                              |
+------------------------------+        +------------------------------+
              \___________________________   __________________________/
                                          \ /
                                           v
+--------------------------------------------------------------+
| Resumo final                                                  |
| - imprime melhores candidatos                                 |
| - imprime sensitivityDense                                    |
| - imprime runs totais                                         |
| - imprime tempo total                                         |
| - diary off                                                   |
+--------------------------------------------------------------+
                              |
                              v
+--------------------------------------------------------------+
|                             FIM                              |
+--------------------------------------------------------------+
````

---

# Mini-fluxo interno do solver principal (`solveAndGetRplusRminus`)

```text
+--------------------------------------------------+
| solveAndGetRplusRminus(...)                      |
+--------------------------------------------------+
                    |
                    v
+--------------------------------------------------+
| Garante existência das tabelas                   |
| - tblRplus                                       |
| - tblRminus                                      |
+--------------------------------------------------+
                    |
                    v
+--------------------------------------------------+
| Calcula número esperado de pontos angulares      |
+--------------------------------------------------+
                    |
                    v
+--------------------------------------------------+
| Executa sweep com m = +1                         |
| - redireciona numericals para tblRplus           |
| - limpa tabela                                   |
| - configura sweep(alpha, m=+1)                   |
| - roda estudo                                    |
| - atualiza resultados                            |
| - lê alpha e Rplus                               |
+--------------------------------------------------+
                    |
                    v
+--------------------------------------------------+
| Executa sweep com m = -1                         |
| - redireciona numericals para tblRminus          |
| - limpa tabela                                   |
| - configura sweep(alpha, m=-1)                   |
| - roda estudo                                    |
| - atualiza resultados                            |
| - lê alpha e Rminus                              |
+--------------------------------------------------+
                    |
                    v
+--------------------------------------------------+
| Valida consistência                              |
| - mesmo número de pontos                         |
| - mesma grade de alpha                           |
+--------------------------------------------------+
                    |
                    v
+--------------------------------------------------+
| Calcula TMOKE                                    |
| TMOKE = (Rplus - Rminus) / (Rplus + Rminus)      |
| com proteção contra divisão por zero             |
+--------------------------------------------------+
                    |
                    v
+--------------------------------------------------+
| Retorna                                           |
| - alpha_deg                                       |
| - Rplus                                           |
| - Rminus                                          |
| - TMOKE                                           |
+--------------------------------------------------+
```