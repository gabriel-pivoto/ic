# senseAndTmoke.m

MATLAB + COMSOL workflow that performs a 5D hierarchical search to maximize both the TMOKE response and the angular sensitivity `S = d(alpha_peak)/dn`. The script drives COMSOL sweeps, checkpoints progress, and exports the best candidate with validation curves and plots.

## Workflow (per stage)
- **COARSE**: Sweep the coarse grids (`domainPeriodGridNm`, `toothWidthGridNm`, `siliconHeightGridNm`, `ceyigHeightGridNm`, `goldHeightGridNm`). For each geometry, sweep two refractive indices from `fastRefractiveIndexSamples`, compute `maxAbsTMOKE_base` (at `baselineRefractiveIndex`) and the quick sensitivity estimate `sensitivityEstimateFast`. Keep top `topKCoarse` seeds by combined rank.
- **FINE**: Clamp small windows around each coarse seed using the `fine*Delta/Step` parameters, shrink the alpha window, recompute the same metrics, and keep `topKFine` seeds.
- **SUPER**: Final refinement around each fine seed using the `super*Delta/Step` windows and the tighter alpha step.
- **VALID**: On the best trade-off candidate, run dense alpha sweeps for every `validationRefractiveIndexList` value and fit `sensitivityDense`.
- **SNAPSHOT (optional)**: Save the COMSOL model at the best trade-off geometry/angle (`snapshotFilePath`).
- **FULL**: Baseline alpha sweep at `baselineRefractiveIndex`, export CSV/XLSX, and create plots (optionally saved to `figureOutputDirectory`).

## Candidate selection (who wins)
- **Best TMOKE only**: In SUPER, `selectTopK_single` picks the row with the largest `maxAbsTMOKE_base` and stores it as `bestTmokeCandidate`.
- **Best sensitivity only**: In SUPER, `selectTopK_single_abs` picks the row with the largest `|S_est_deg_per_RIU|` and stores it as `bestSensitivityCandidate`.
- **Best trade-off (TMOKE + sensitivity)**: In SUPER, `selectTopK_tradeoff` ranks candidates by `|TMOKE|` and by `|S_est|` (rank 1 = best in each list) and sums those ranks. The lowest sum wins because it signals a candidate that scores well on both metrics (high TMOKE without sacrificing sensitivity). The winner is stored as `bestTradeoffCandidate`.

## Key variables
- `projectRootDir`: Root folder holding the COMSOL model, checkpoints, and outputs.  
- `comsolModelFile`: Path to the `.mph` model loaded by the script.  
- `MAX_RUNS`: Safety guard; abort if planned runs exceed this budget.  
- `checkpointDirectory`, `checkpointFilePath`, `progressWorkbookPath`, `checkpointEveryPoints`: Checkpoint/resume settings and progress workbook.  
- `fastRefractiveIndexSamples`: Two `n` values used during search to compute `sensitivityEstimateFast`.  
- `baselineRefractiveIndex`: Baseline `n` (first element of `fastRefractiveIndexSamples`) used for TMOKE scoring and final sweep.  
- `validationRefractiveIndexList`: Refractive indices used in the VALID stage to compute `sensitivityDense`.  
- Search grids (nm): `domainPeriodGridNm`, `toothWidthGridNm`, `siliconHeightGridNm`, `ceyigHeightGridNm`, `goldHeightGridNm`.  
- Alpha sweeps: `alphaCoarseRange` (start/step/stop), `alphaFineStep`, `alphaSuperStep`, `alphaDenseStep`, `alphaFullStep`.  
- Promotion strategy: `topKCoarse`, `topKFine`.  
- Refinement windows: `fineGoldHeightDelta/Step`, `fineCeyigHeightDelta/Step`, `fineDomainPeriodDelta/Step`, `fineToothWidthDelta/Step`, `fineSiliconHeightDelta/Step`; similarly `super*` for the SUPER stage.  
- Alpha window sizes: `fineAlphaHalfSpan`, `superAlphaHalfSpan`, and the wider `fineAlphaHalfSpanSensitivity` / `superAlphaHalfSpanSensitivity` used when estimating sensitivity to avoid missing the peak.  
- Plot/export flags: `SAVE_SNAPSHOT`, `MAKE_PLOTS`, `PLOT_LIVE`, `SAVE_FIGS`, `figureOutputDirectory`, `FIG_FORMATS`, `runTimestamp`.  
- COMSOL parameter/tag names: `STUDY_TAG`, `PARAM_HAU`, `PARAM_HCEY`, `PARAM_LDOM`, `PARAM_LDEN`, `PARAM_HSI`, `PARAM_N`, `ALPHA_NAME`, `MSIGN_NAME`, plus the reflectance tables `RPLUS_TABLE_TAG` and `RMINUS_TABLE_TAG`.

## Outputs
- Checkpoint: `checkpointFilePath` stores `checkpointData` (stage, run counters, partial tables); auto-saved every `checkpointEveryPoints` points.
- Progress workbook: `progressWorkbookPath` gathers coarse/fine/super/dense/full tables for quick inspection.
- CSV exports: coarse/fine/super tables plus `tmoke_sens_bestTradeoff_dense_ALLn.csv` and `tmoke_sens_bestTradeoff_full_baseline.csv`.
- Snapshot: optional `.mph` copy of the best trade-off geometry at `alphaBestDeg` and `baselineRefractiveIndex`.
- Plots: live plotting during runs (if `PLOT_LIVE`) and final figures saved under `figureOutputDirectory` when `SAVE_FIGS` is true.

## Functions (what each helper does)
- `setParamNm` / `setParamScalar`: set COMSOL parameters with or without units (nm vs unitless).
- `selectTopK_tradeoff` / `selectTopK_single` / `selectTopK_single_abs`: choose the best rows from a table using combined TMOKE+S rank, raw max, or abs max.
- `solveAndGetRplusRminus`: run COMSOL parametric sweeps for `alpha` and magnetization `m`, read reflectance tables, and compute TMOKE.
- `setTwoParamSweep` / `setAlphaMSweep`: configure the COMSOL study for two-parameter sweeps (alpha, m) with the desired range/list.
- `getParametricTag`: find the parametric sweep feature tag inside the loaded COMSOL study.
- `refreshDerivedValues`: force COMSOL to refresh numerical results after a run.
- `ensureTable` / `clearTable` / `redirectAllNumericalsToTable`: create/clear tables and redirect COMSOL numerical outputs to the chosen table tag.
- `readAlphaAndRFromNamedTable`: parse alpha and reflectance columns from a COMSOL result table and return them as vectors.
- `fmt_time_long` / `frac_eta` / `fmt_eta_flag` / `fmt_pct`: format elapsed/ETA/time strings for logs.
- `logline`: write formatted progress lines and flush the UI.
- `save_checkpoint` / `write_progress_xlsx` / `maybe_checkpoint`: persist intermediate progress to `.mat` and the Excel workbook, with periodic checkpoints.
- `updateLivePlot`: live TMOKE/reflectance plot used during sweeps.
- `isvalidHandle`: small guard to test MATLAB graphics handles.
- `saveAllOpenFigures`: export every open figure to disk in the configured formats.
- `sanitizeFilename`: strip unsafe characters from figure names when exporting.

## Notes
- Resuming: if `checkpointFilePath` exists, the script restores `checkpointData` and resumes the appropriate stage.
- Units: geometry grids are in nanometers; alpha is in degrees; refractive indices are unitless.
- Tuning: adjust the `*Delta/Step` window sizes, alpha steps, or TOP-K values to widen/narrow the search without changing the overall flow.
