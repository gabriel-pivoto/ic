# Fluxograma detalhado do script MATLAB + COMSOL  
## 5D Hierarchical Search para otimização de TMOKE e sensibilidade angular
https://drive.google.com/drive/folders/1J6QObtXH-x2FjjVxw428po6ZvkUaICzW?usp=sharing
---

## Estrutura do repositório
# sensitivity.m — Angular sensitivity optimization (MATLAB + COMSOL)

4D hierarchical search that drives a COMSOL model to find the geometry that maximizes angular sensitivity

```
S = d(alpha_peak)/dn   [deg/RIU]
```

where `alpha_peak` is the incidence angle at which `|TMOKE(alpha)|` is maximal, and `n` is the refractive index of the external medium. The script is single-objective: there is no optimization or selection by `|TMOKE|`, only by `|S|`. TMOKE(alpha) is still computed internally only because `alpha_peak` is defined as the peak of that curve — it is the means to locate the tracked resonance, not a competing objective.

## Repository

```
ic/
├── README.md
└── sensitivity.m
```

The rest of the project (other MATLAB scripts, docs, manuscript) exists locally but is out of version control (`.gitignore`: `docs/`, `manuscript/`, `src/`).

## Requirements

- MATLAB with COMSOL's LiveLink for MATLAB configured (`com.comsol.model.*`).
- A `.mph` model with:
  - parameters `h_au`, `L_domain`, `l_dente`, `h_si`, `n`;
  - a parametric sweep on `std1` over `alpha` [deg] and `m` (magnetization, `1`/`-1`);
  - numerical Derived Values redirectable to tables `tblTplus` / `tblTminus`, containing angle and total-transmission columns.
- Set `projectRootDir` (line ~55) and `comsolModelFile` to the local project/`.mph` path.

## How sensitivity is computed

For a set of refractive indices `n` (e.g. `[1.30, 1.33, 1.36]`), with `n = 1.33` as reference:

1. Solve the reference curve (`n = 1.33`) and locate the global peak of `|TMOKE(alpha)|` → this defines `alpha_ref` and the sign of the peak.
2. For each non-reference `n`, solve `alpha` only inside the window `alpha_ref ± trackingHalfWindowDeg` (20°, clipped to the global alpha limits) and pick the peak within that window — preferring the same sign as the reference peak.
3. Fit a line `alpha_peak(n)` by linear regression (`polyfit`); the slope is `S`.

This "tracking" keeps the metric from jumping between different resonances as `n` changes — which would happen if each `n` were solved independently and the global peak of each curve were used instead.

TMOKE itself is computed as:

```
TMOKE(alpha) = 2 * (T+(alpha) - T-(alpha)) / (T+(alpha) + T-(alpha))
```

with `T+`/`T-` being the total transmission for `m = +1` and `m = -1`, read from the COMSOL tables `tblTplus`/`tblTminus`.

## Geometry (4D search)

| Parameter | Meaning | Stage 1 grid |
|---|---|---|
| `L_domain` | domain period [nm] | `800:50:850` (2 pts) |
| `l_dente` | tooth width [nm] | `500:50:600` (3 pts) |
| `h_si` | silicon height [nm] | `[220, 240, 260]` (3 pts) |
| `h_au` | gold height [nm] | `20:10:60` (5 pts) |

Stage 1 total: `2×3×3×5 = 90` geometry points. Each point costs `2 × numel(fastRefractiveIndexSamples)` COMSOL runs (2 magnetizations × 3 refractive indices = 6 runs/point), so `540` runs for Stage 1 alone.

## Refractive indices

- `fastRefractiveIndexSamples = [1.30, 1.33, 1.36]` — used in Stages 1–3.
- `trackingReferenceRefractiveIndex = 1.33` — reference curve used to lock the resonance.
- `trackingHalfWindowDeg = 20` — half-width of the search window around the reference peak.
- `baselineRefractiveIndex = 1.33` — used for the final snapshot.
- `validationRefractiveIndexList = [1.30, 1.33, 1.36]` — used in Stage 4 (same list, much denser alpha).

## Hierarchical pipeline

```
Stage 1 → pick TOP-K (topKCoarse=1) by |S|
  → Stage 2 → pick TOP-K (topKFine=1) by |S|
    → Stage 3 → bestSensitivityCandidate (largest |S|)
      → Stage 4 (dense curves, sensitivityDense)
        → Snapshot .mph (optional)
          → Export CSV/XLSX + figures
```

| Stage | Angular step | Geometry window around the seed | Angular window |
|---|---|---|---|
| Stage 1 | 1.0° (`0`–`89`) | full grid above | `0`–`89` |
| Stage 2 | 0.1° | `±10 nm` (L_domain/l_dente), `±5 nm` (h_si), `±2 nm` (h_au) | `alpha_peak_base ± 6°` |
| Stage 3 | 0.01° | `±4 nm` (L_domain/l_dente), `±2 nm` (h_si), `±1 nm` (h_au) | `alpha_peak_base ± 4°` |
| Stage 4 | 0.01° | fixed geometry (`bestSensitivityCandidate`) | full `0`–`89`, all `n` |

Seed selection at every stage: `selectTopK_single_abs` — sorts by `abs(S_est_deg_per_RIU)` descending and takes the top `K`. There is no combined score and no selection by TMOKE.

## Checkpoint and resume

- `checkpointFilePath`: `<projectRootDir>/checkpoints/sensitivity_4d_checkpoint.mat`
- `progressWorkbookPath`: `<projectRootDir>/checkpoints/sensitivity_4d_progress.xlsx`
- `checkpointSchemaTag = 'tracked_same_peak_sensitivity_v1'` — checkpoints from another version/metric are discarded (with a warning) instead of reused.
- Saved every `checkpointEveryPoints = 10` evaluated points.
- On restart, the script detects the saved `stage` (`Stage 2`/`Stage 3`/`Stage 4`/`FINAL`) and skips completed stages, restoring the matching result tables and seeds.
- The checkpoint file is deleted at the end of a full run.
- `MAX_RUNS = 20000`: the script aborts (`error`) if the planned total run count (exact from Stage 2 planning onward) exceeds this cap.

## Outputs

- **CSV** (in `projectRootDir`): `sensitivity_4d_coarse.csv`, `sensitivity_4d_fine.csv`, `sensitivity_4d_super.csv`, `sensitivity_4d_bestSensitivity_dense_ALLn.csv`.
- **Progress XLSX**: sheets `coarse`, `fine`, `super`, `dense`.
- **Snapshot `.mph`** (if `SAVE_SNAPSHOT = true`): geometry of `bestSensitivityCandidate` at `n = baselineRefractiveIndex`, `alpha = alphaBestDeg`.
- **Figures** (if `MAKE_PLOTS`/`SAVE_FIGS = true`) under `<projectRootDir>/plots/sensitivity_4d_<timestamp>/`:
  - `phase_outputs/`: `|S|` metric per candidate at each stage, dense curves from Stage 4.
  - `iteration_tmoke/`: one PNG per evaluated point, if `PLOT_LIVE` and `SAVE_ITER_PLOTS` are enabled.
  - Final figures: `TMOKE(alpha)` per `n`, `alpha_peak vs n` with linear fit, `|TMOKE|` at the tracked peak vs `n`.

## Key functions

- `evaluateTrackedSensitivityAndCurves(...)` — orchestrates the resonance-tracking logic described above; returns `alphaPeakDegreesByN`, `trackedTmokeAbsByN`, `linearFit`, and `sensitivitySlope` (= `S`).
- `solveAndGetTplusTminus(...)` — runs COMSOL for `m=+1` and `m=-1` over the requested `alpha` sweep, reads `tblTplus`/`tblTminus`, and computes `TMOKE(alpha)`.
- `selectTopK_single_abs(T, K, col)` — the sole selection criterion for seeds/final candidate, by descending `abs(T.(col))`.
- `save_checkpoint` / `maybe_checkpoint` / `write_progress_xlsx` — state and progress persistence.
