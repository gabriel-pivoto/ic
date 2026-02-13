# TMOKE + sensitivity search (MATLAB + COMSOL)

Organized runner for the TMOKE + angular-sensitivity search with packages and portable config.

## Setup
1. Copy `config/tmoke_config_template.m` to `config/tmoke_config_local.m` (or edit the existing stub) and set:
   - `comsolModelFile` (path to your `.mph` model)
   - `outputDir` (where CSV/checkpoints/plots go)
   - Any grids/alpha steps/top-K/window sizes you want to tweak.
2. Ensure COMSOL with LiveLink for MATLAB is available on the MATLAB path.

## Run
In MATLAB, from this folder (`pasta para essa organizaçaõ/`):
```matlab
scripts/run_senseAndTmoke
```
The script loads config, ensures output/checkpoint folders, runs COARSE/FINE/SUPER/VALID/FULL, and logs to `checkpoints/runlog_*.txt`.

## Outputs
- Checkpoints + progress workbook: `config` values `checkpointDir` and `progressWorkbookPath` (default under `output/checkpoints`).
- CSV exports and final baseline/dense curves: under `output/`.
- Plots: `output/plots/tmoke_sens_5d_<timestamp>/` when `SAVE_FIGS=true`.
- Optional COMSOL snapshot of the best trade-off when `SAVE_SNAPSHOT=true`.

## Code layout
- `scripts/run_senseAndTmoke.m` — entrypoint.
- `src/+tmoke/runSenseAndTmoke.m` — main flow (no local functions).
- `src/+tmoke/+comsol` — COMSOL parameter/table helpers.
- `src/+tmoke/+selection` — candidate ranking helpers.
- `src/+tmoke/+checkpoint` — checkpoint/progress utilities.
- `src/+tmoke/+plot` — live plot + figure export helpers.
- `src/+tmoke/+util` — logging/formatting/sanitization helpers.
