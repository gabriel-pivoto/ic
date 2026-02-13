# How to run the organized TMOKE workflow

1) Open MATLAB with COMSOL LiveLink available and `cd` into `organizado/`.
2) Configure:
   - Copy `config/tmoke_config_template.m` to `config/tmoke_config_local.m`.
   - Edit the local file to point to your `.mph` (`comsolModelFile`), set `outputDir`, and adjust grids/steps/top-K as needed.
3) Run the entrypoint:
   ```matlab
   scripts/run_senseAndTmoke
   ```
   The runner loads config, prepares output/checkpoint folders, and executes the COARSE → FINE → SUPER → VALID → FULL stages.
4) Results:
   - Checkpoints/logs/xlsx: under `output/checkpoints/` (or the dirs you set in config).
   - CSVs for coarse/fine/super + dense/full curves: under `output/`.
   - Plots (when `MAKE_PLOTS=true`; saved when `SAVE_FIGS=true`): `output/plots/tmoke_sens_5d_<timestamp>/`.
   - Optional COMSOL snapshot of the best trade-off when `SAVE_SNAPSHOT=true`.
