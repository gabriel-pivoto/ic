%% =====================================================================
% MATLAB + COMSOL - 5D Hierarchical Search (TMOKE + angular sensitivity)
% ---------------------------------------------------------------------
% Goal (one line):
%   Find the 5D geometry that simultaneously maximizes (A) |TMOKE| and
%   (B) angular sensitivity S = d(alpha_peak)/dn [deg/RIU].
%
% How this script merges the two original pipelines:
%   1) TMOKE optimization (max |TMOKE|)      -> cheap (2 runs per point)
%   2) Sensitivity S(n) optimization         -> expensive (>= 2*#n runs)
%
% For every geometry point (COARSE/FINE/SUPER) we compute both metrics:
%   - TMOKE metric:    maxAbsTMOKE_base = max(|TMOKE(alpha)|) at baselineRefractiveIndex
%   - Sensitivity est: sensitivityEstimateFast = (alpha_peak(n2) - alpha_peak(n1)) / (n2 - n1)
%                      using only two refractive indices (fast estimate).
%
% In the VALID stage, on the best candidate, we recompute sensitivity with
% a denser set of refractive indices:
%   - validationRefractiveIndexList = [1.33 1.36 1.39] (configurable)
%   - sensitivityDense = slope of alpha_peak vs n (deg/RIU) from linear fit.
%
% Trade-off selection:
%   In each stage (COARSE/FINE/SUPER) candidates are ranked by:
%     - rank_TM = rank of |TMOKE|max (higher is better)
%     - rank_S  = rank of |sensitivityEstimateFast| (higher is better)
%     - score_comb = rank_TM + rank_S  (smaller score = better trade-off)
%   This yields:
%     - bestTradeoffCandidate: best compromise between TMOKE and S
%     - bestTmokeCandidate: best only in TMOKE
%     - bestSensitivityCandidate: best only in sensitivity
%
% High-level flow:
%   Start -> Load .mph -> (optional resume from checkpoint)
%     -> COARSE (sweep wide grid, compute TMOKE + sensitivityEstimateFast)
%     -> pick TOP-K seeds (trade-off score)
%     -> FINE   (refine geometry + alpha window, recompute metrics)
%     -> pick TOP-K seeds (trade-off score)
%     -> SUPER  (final refinement, recompute metrics)
%     -> choose bestTradeoffCandidate + bestTmokeCandidate + bestSensitivityCandidate
%     -> VALID  (dense curves on bestTradeoffCandidate to compute sensitivityDense)
%     -> FULL   (final baseline curve + CSV/XLSX export)
%     -> SAVE FIGURES
%     -> End
% =====================================================================
clear; clc; close all; format long; tic   % Reset state, close figures, use long format, start overall timer
import com.comsol.model.*;               % Import COMSOL model classes
import com.comsol.model.util.*;          % Import COMSOL utility helpers

%% --------------------------- Paths/Project -------------------------
projectRootDir  = 'C:\Users\gabri\Documents\projetoIC';       % Root folder that holds the COMSOL project and outputs
comsolModelFile = fullfile(projectRootDir,'usandoMatlab.mph');          % COMSOL model file to open
addpath(genpath(projectRootDir));                                       % Expose all helper scripts in the project to MATLAB

%% ------------------------ Run budget (optional) ---------------------
% Keep the spirit of the original TMOKE script: stop if the total planned
% number of COMSOL runs exceeds this budget (fast sensitivity doubles cost).
MAX_RUNS = 20000;

%% ------------------------ Checkpoint config -------------------------
checkpointDirectory   = fullfile(projectRootDir,'checkpoints');         % Folder where checkpoints and logs are written
if ~exist(checkpointDirectory,'dir'), mkdir(checkpointDirectory); end   % Create checkpoint folder if it does not exist
checkpointFilePath    = fullfile(checkpointDirectory,'tmoke_sens_5d_checkpoint.mat'); % Binary checkpoint file
progressWorkbookPath  = fullfile(checkpointDirectory,'tmoke_sens_progress.xlsx');     % Human-friendly Excel progress log

checkpointEveryPoints = 10;                                             % Save progress every N evaluated geometry points
pointsSinceCheckpoint = 0;                                              % Counter since the last checkpoint write

%% ------------------------ Resume / Runlog ---------------------------
resumeFromCheckpoint = false;                                           % Flag: resume from previous checkpoint+/-
checkpointData = struct;                                                % Struct holding restored checkpoint payload
resumeStageTag = "";                                                    % Stage name captured when resuming
if isfile(checkpointFilePath)
    % If a checkpoint exists, restore run counters/tables to continue where the last run stopped.
    try
        checkpointFileContents = load(checkpointFilePath,'checkpointData'); % Load only the checkpointData variable
        checkpointData = checkpointFileContents.checkpointData;             % Extract saved struct from disk
        resumeFromCheckpoint = true;                                        % Flip resume flag because checkpoint exists
        resumeStageTag = string(checkpointData.stage);                      % Stage name recorded when the checkpoint was saved
        fprintf('>>> resumeFromCheckpoint detected | stage=%s | done_points=%d | runsCompletedGlobal=%d\\n', ...
            resumeStageTag, checkpointData.done_points, checkpointData.runsCompletedGlobal);
    catch
        warning('Checkpoint exists but could not be loaded. Starting fresh.');
    end
end

diary(fullfile(checkpointDirectory, sprintf('runlog_%s.txt', datestr(now,'yyyymmdd_HHMMss')))); % Persistent run log

%% ------------------------ Model Tags/Params -------------------------
STUDY_TAG   = 'std1';

PARAM_HAU   = 'h_au';
PARAM_HCEY  = 'h_ceyig';
PARAM_LDOM  = 'L_domain';
PARAM_LDEN  = 'l_dente';
PARAM_HSI   = 'h_si';
PARAM_N     = 'n';               % <<< vem do seu script de sensibilidade

ALPHA_NAME  = 'alpha';
MSIGN_NAME  = 'm';
M_PLUS      = '1';
M_MINUS     = '-1';

RPLUS_TABLE_TAG  = 'tblRplus';
RMINUS_TABLE_TAG = 'tblRminus';

%% --------------------- Refractive indices (n) ----------------------
% During the search stages (COARSE/FINE/SUPER):
% - use 2 refractive indices to estimate S quickly (sensitivityEstimateFast)
% - baselineRefractiveIndex is also used for the TMOKE metric (maxAbsTMOKE_base)
fastRefractiveIndexSamples = [1.33, 1.39];           % Two n values used for the quick sensitivity estimate
baselineRefractiveIndex    = fastRefractiveIndexSamples(1); % Baseline n for TMOKE sweeps (change if needed)

% During the final validation stage (VALID):
% - use 3 (or more) refractive indices to compute sensitivityDense reliably
validationRefractiveIndexList = [1.33, 1.36, 1.39];


%% ------------------------ COARSE grids (exact) ----------------------
domainPeriodGridNm = 800:50:850;         % 3
toothWidthGridNm  = 500:50:600;         % 3
siliconHeightGridNm     = [220, 240, 260];    % 3
ceyigHeightGridNm    = [100, 140];         % 2
goldHeightGridNm     = 20:10:60;           % 5
% COARSE points = 270

%% ---------------------- Alpha ranges/steps --------------------------
alphaCoarseRange = [0, 1.0, 89];  % Coarse sweep: start [deg], step [deg], stop [deg]
alphaFineStep  = 0.1;             % Fine alpha step [deg]
alphaSuperStep = 0.01;            % Superfine alpha step [deg]
alphaDenseStep = 0.01;            % Validation alpha step [deg]
alphaFullStep  = 0.01;            % Final baseline alpha step [deg]

%% ----------------------- TOP-K strategy ----------------------------
topKCoarse = 1;   % Number of coarse candidates promoted to the fine stage
topKFine   = 1;   % Number of fine candidates promoted to the superfine stage

%% --------------------- Refinement windows --------------------------
fineGoldHeightDelta   = 2;    fineGoldHeightStep   = 1;
fineCeyigHeightDelta  = 5;    fineCeyigHeightStep  = 1;
fineDomainPeriodDelta = 10;   fineDomainPeriodStep = 5;
fineToothWidthDelta   = 10;   fineToothWidthStep   = 5;
fineSiliconHeightDelta= 5;    fineSiliconHeightStep= 5;

superGoldHeightDelta  = 1;    superGoldHeightStep  = 0.5;
superCeyigHeightDelta = 2;    superCeyigHeightStep = 1;
superDomainPeriodDelta = 4;   superDomainPeriodStep = 2;
superToothWidthDelta   = 4;   superToothWidthStep   = 2;
superSiliconHeightDelta= 2;   superSiliconHeightStep= 1;

%% ------------------ Alpha windows (centered at last peak) ----------
fineAlphaHalfSpan   = 5;     % Nominal TMOKE window half-span for FINE stage [deg]
superAlphaHalfSpan  = 2;     % Nominal TMOKE window half-span for SUPER stage [deg]

% The alpha_peak can shift slightly with refractive index, so use slightly
% wider windows when estimating sensitivity to avoid missing the peak.
fineAlphaHalfSpanSensitivity  = fineAlphaHalfSpan  + 1;   % +/-6 deg margin for sensitivity sweeps
superAlphaHalfSpanSensitivity = superAlphaHalfSpan + 2;   % +/-4 deg margin for sensitivity sweeps

%% ----------------------- Outputs / Plots ----------------------------
SAVE_SNAPSHOT = true;   % Save a COMSOL snapshot of the final bestTradeoffCandidate
MAKE_PLOTS    = true;   % Render final plots after computations
PLOT_LIVE     = true;   % Update a live TMOKE figure during sweeps (slower)

%% ----------------------- Save figures (FINAL) -----------------------
SAVE_FIGS   = true;                                     % Persist all open figures to disk at the end
FIG_FORMATS = {'png','pdf'};                            % Formats to export for each figure
runTimestamp   = datestr(now,'yyyy-mm-dd_HHMMSS');      % Unique run tag for plot folders
figureOutputDirectory     = fullfile(projectRootDir,'plots', ['tmoke_sens_5d_' runTimestamp]);
if SAVE_FIGS && ~exist(figureOutputDirectory,'dir'), mkdir(figureOutputDirectory); end

%% ------------------------- Load model ------------------------------
ModelUtil.clear;                        % Clear any previous COMSOL models from memory
model = mphload(comsolModelFile);       % Load the COMSOL model to be driven
ModelUtil.showProgress(true);           % Show COMSOL progress in the MATLAB console

%% ==================================================================
%                    GLOBAL PLANNING (runs / ETA)
% ------------------------------------------------------------------
% Compute exact/estimated run counts for each stage to:
%   - enforce the MAX_RUNS budget,
%   - show stage/global ETAs while sweeping,
%   - size refinement windows up front.
% These counters are updated as each point is solved.
% ==================================================================
globalTimerStart = tic;
runsCompletedGlobal = 0;
isTotalRunEstimateExact   = false;

numDomainPeriodPoints = numel(domainPeriodGridNm);      % Number of values in the domain-period grid
numToothWidthPoints   = numel(toothWidthGridNm);        % Number of tooth-width samples
numSiliconHeightPoints= numel(siliconHeightGridNm);     % Number of silicon-height samples
numCeyigHeightPoints  = numel(ceyigHeightGridNm);       % Number of Ce:YIG-height samples
numGoldHeightPoints   = numel(goldHeightGridNm);        % Number of gold-height samples

% Each refractive index requires 2 runs (m = +/-1). Searching uses two n values.
runsPerSearchPoint = 2 * numel(fastRefractiveIndexSamples);     % 4 runs per geometry point

coarseTotalPoints  = numDomainPeriodPoints*numToothWidthPoints*numSiliconHeightPoints*numCeyigHeightPoints*numGoldHeightPoints;
coarseTotalRuns = runsPerSearchPoint * coarseTotalPoints;

% Fine-stage estimates (before exact clamping)
fineGoldHeightCount = 1 + floor(2*fineGoldHeightDelta   / fineGoldHeightStep);
fineCeyigHeightCount = 1 + floor(2*fineCeyigHeightDelta  / fineCeyigHeightStep);
fineDomainPeriodCount = 1 + floor(2*fineDomainPeriodDelta  / fineDomainPeriodStep);
fineToothWidthCount = 1 + floor(2*fineToothWidthDelta  / fineToothWidthStep);
fineSiliconHeightCount = 1 + floor(2*fineSiliconHeightDelta   / fineSiliconHeightStep);
finePointsPerSeedEstimate = fineDomainPeriodCount*fineToothWidthCount*fineSiliconHeightCount*fineCeyigHeightCount*fineGoldHeightCount;
fineTotalPointsEstimate    = topKCoarse * finePointsPerSeedEstimate;
fineTotalRunsEstimate   = runsPerSearchPoint * fineTotalPointsEstimate;

% Superfine counts (exact upfront)
superGoldHeightCount = 1 + floor(2*superGoldHeightDelta   / superGoldHeightStep);
superCeyigHeightCount = 1 + floor(2*superCeyigHeightDelta  / superCeyigHeightStep);
superDomainPeriodCount = 1 + floor(2*superDomainPeriodDelta  / superDomainPeriodStep);
superToothWidthCount = 1 + floor(2*superToothWidthDelta  / superToothWidthStep);
superSiliconHeightCount = 1 + floor(2*superSiliconHeightDelta   / superSiliconHeightStep);
superPointsPerSeed = superDomainPeriodCount*superToothWidthCount*superSiliconHeightCount*superCeyigHeightCount*superGoldHeightCount;
superTotalPoints    = topKFine * superPointsPerSeed;
superTotalRuns   = runsPerSearchPoint * superTotalPoints;

% Extra runs that always happen:
% - VALID: 2 runs per validationRefractiveIndexList element (m = +/-1) => 2*len(validationRefractiveIndexList)
% - snapshot: 2 runs (m = +/-1) at a fixed alpha in baselineRefractiveIndex
% - FULL: 2 runs (m = +/-1) full sweep at baselineRefractiveIndex
extraRunsFixed = 2*numel(validationRefractiveIndexList) + (SAVE_SNAPSHOT*2) + 2;

globalRunTargetEstimate = coarseTotalRuns + fineTotalRunsEstimate + superTotalRuns + extraRunsFixed;
isTotalRunEstimateExact = false;

fprintf('START\n');
fprintf('  COARSE (exact):    %d runs\n', coarseTotalRuns);
fprintf('  FINE   (estimate): %d runs (TOP-%d)\n', fineTotalRunsEstimate, topKCoarse);
fprintf('  SUPER  (exact):    %d runs (TOP-%d)\n', superTotalRuns, topKFine);
fprintf('  EXTRAS (exact):    %d runs\n', extraRunsFixed);
fprintf('  TOTAL  (estimate): %d runs\n\n', globalRunTargetEstimate);

%% ==================================================================
%                               COARSE
% ==================================================================
% Sweep the full coarse geometry grid. For each geometry point:
%   1) Set geometry parameters.
%   2) Sweep alpha at fastRefractiveIndexSamples(1) to find TMOKE peak and peak angle.
%   3) Sweep alpha at fastRefractiveIndexSamples(2) to estimate sensitivity from the peak shift.
%   4) Log metrics, update ETAs, and checkpoint periodically.
% Outputs:
%   - coarseResultsTable: metrics for every coarse point.
%   - coarseSeedCandidates: TOP-K seeds ranked by combined TMOKE/sensitivity score.
coarseResultsTable = [];
coarseSeedCandidates   = [];

skipCoarseStage = resumeFromCheckpoint && any(strcmp(resumeStageTag, ["FINE","SUPER","VALID","FULL"]));
if skipCoarseStage
    % When resuming past COARSE, restore previous coarse results/seeds instead of recomputing.
    if isfield(checkpointData.payload,'coarseResultsTable'), coarseResultsTable = checkpointData.payload.coarseResultsTable; end
    if isfield(checkpointData.payload,'coarseSeedCandidates'),   coarseSeedCandidates   = checkpointData.payload.coarseSeedCandidates;   end
    fprintf('SKIP COARSE -> restored coarseResultsTable (rows=%d), coarseSeedCandidates (rows=%d)\n', ...
        size(coarseResultsTable,1), height(coarseSeedCandidates));
else
    fprintf('STAGE COARSE - EXACT: %d runs\n', coarseTotalRuns);  % Full coarse grid, no estimates here
    stageRunsStart = runsCompletedGlobal;
    stageTotalRuns = coarseTotalRuns;
    stageTimerStart = tic;

    % Colunas:
    % [Ldom, Lden, hsi, hcey, hau,
    %  maxAbsTMOKE_base, alpha_peak_base, TMOKE_at_peak_base,
    %  alpha_peak_n2, sensitivityEstimateFast]
    coarseRows = [];
    coarsePointIndex = 0;

    if resumeFromCheckpoint && resumeStageTag=="COARSE"
        runsCompletedGlobal = checkpointData.runsCompletedGlobal;
        if isfield(checkpointData.payload,'coarseRows'), coarseRows = checkpointData.payload.coarseRows; end
        coarsePointIndex = checkpointData.done_points;
        fprintf('>>> COARSE resume: skipping %d points already computed.\n', coarsePointIndex);
    end

    for domainPeriodIdx = 1:numDomainPeriodPoints
        setParamNm(model, PARAM_LDOM, domainPeriodGridNm(domainPeriodIdx));
        for toothWidthIdx = 1:numToothWidthPoints
            setParamNm(model, PARAM_LDEN, toothWidthGridNm(toothWidthIdx));
            for siliconHeightIdx = 1:numSiliconHeightPoints
                setParamNm(model, PARAM_HSI, siliconHeightGridNm(siliconHeightIdx));
                for ceyigHeightIdx = 1:numCeyigHeightPoints
                    setParamNm(model, PARAM_HCEY, ceyigHeightGridNm(ceyigHeightIdx));
                    for goldHeightIdx = 1:numGoldHeightPoints
                        setParamNm(model, PARAM_HAU, goldHeightGridNm(goldHeightIdx));

                        coarsePointIndex = coarsePointIndex + 1;
                        if resumeFromCheckpoint && resumeStageTag=="COARSE" && coarsePointIndex <= checkpointData.done_points
                            continue;
                        end

                        % =======================
                        % (1) n = fastRefractiveIndexSamples(1) (base)
                        % =======================
                        setParamScalar(model, PARAM_N, fastRefractiveIndexSamples(1)); % n sem unidade
                        [alphaGridFastN1, reflectancePlusFastN1, reflectanceMinusFastN1, tmokeCurveFastN1] = solveAndGetRplusRminus( ...
                            model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                            alphaCoarseRange(1), alphaCoarseRange(2), alphaCoarseRange(3), ...
                            M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                        [tmokePeakMagnitudeFastN1, tmokePeakIndexFastN1] = max(abs(tmokeCurveFastN1));
                        alphaAtPeakFastN1 = alphaGridFastN1(tmokePeakIndexFastN1);
                        tmokeAtPeakFastN1 = tmokeCurveFastN1(tmokePeakIndexFastN1);

                        if tmokePeakIndexFastN1==1 || tmokePeakIndexFastN1==numel(tmokeCurveFastN1)
                            logline('WARN COARSE (n=%.2f) peak at alpha-edge (alpha*=%.3f in [%.3f, %.3f])\n', ...
                                fastRefractiveIndexSamples(1), alphaAtPeakFastN1, alphaGridFastN1(1), alphaGridFastN1(end));
                        end

                        if PLOT_LIVE
                            updateLivePlot('COARSE', ...
                                domainPeriodGridNm(domainPeriodIdx), toothWidthGridNm(toothWidthIdx), siliconHeightGridNm(siliconHeightIdx), ...
                                ceyigHeightGridNm(ceyigHeightIdx), goldHeightGridNm(goldHeightIdx), ...
                                alphaGridFastN1, tmokeCurveFastN1, alphaAtPeakFastN1, tmokeAtPeakFastN1, reflectancePlusFastN1, reflectanceMinusFastN1);
                        end

                        % ======================
                        % (2) n = fastRefractiveIndexSamples(2) (alphaVsNLinearFit/ sensitivityEstimateFast)
                        % ======================
                        setParamScalar(model, PARAM_N, fastRefractiveIndexSamples(2));
                        [alphaGridFastN2, ~, ~, tmokeCurveFastN2] = solveAndGetRplusRminus( ...
                            model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                            alphaCoarseRange(1), alphaCoarseRange(2), alphaCoarseRange(3), ...
                            M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                        [~, tmokePeakIndexFastN2] = max(abs(tmokeCurveFastN2));
                        alphaAtPeakFastN2 = alphaGridFastN2(tmokePeakIndexFastN2);

% Quick sensitivity estimate using 2 points
                        sensitivityEstimateFast = (alphaAtPeakFastN2 - alphaAtPeakFastN1) / (fastRefractiveIndexSamples(2) - fastRefractiveIndexSamples(1));

                        % Store candidate point data
                        coarseRows = [coarseRows; ...
                            domainPeriodGridNm(domainPeriodIdx), toothWidthGridNm(toothWidthIdx), siliconHeightGridNm(siliconHeightIdx), ...
                            ceyigHeightGridNm(ceyigHeightIdx), goldHeightGridNm(goldHeightIdx), ...
                            tmokePeakMagnitudeFastN1, alphaAtPeakFastN1, tmokeAtPeakFastN1, ...
                            alphaAtPeakFastN2, sensitivityEstimateFast]; %#ok<AGROW>

% Point cost: 2*fastRefractiveIndexSamples runs (m=+/-1 for each n)
                        runsCompletedGlobal = runsCompletedGlobal + runsPerSearchPoint;

                        % ETA
                        stageRunsCompleted = runsCompletedGlobal - stageRunsStart;
                        [stageCompletionFraction, stageEtaSeconds]   = frac_eta(stageRunsCompleted, stageTotalRuns, toc(stageTimerStart));
                        [globalCompletionFraction, globalEtaSeconds] = frac_eta(runsCompletedGlobal, globalRunTargetEstimate, toc(globalTimerStart));

                        logline(['COARSE | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                 ' | |TM|=%.5f @ alpha=%.2f deg | sensitivityEstimateFast=%.4f deg/RIU' ...
                                 ' [Stage %5.1f%% | ETA %s | Global %s%% | ETA %s]\n'], ...
                            domainPeriodGridNm(domainPeriodIdx), toothWidthGridNm(toothWidthIdx), siliconHeightGridNm(siliconHeightIdx), ...
                            ceyigHeightGridNm(ceyigHeightIdx), goldHeightGridNm(goldHeightIdx), ...
                            tmokePeakMagnitudeFastN1, alphaAtPeakFastN1, sensitivityEstimateFast, ...
                            100*stageCompletionFraction, fmt_time_long(stageEtaSeconds), fmt_pct(globalCompletionFraction), fmt_eta_flag(globalEtaSeconds, isTotalRunEstimateExact));

                        % Checkpoint
                        pointsSinceCheckpoint = pointsSinceCheckpoint + 1;
                        payload = struct('coarseRows',coarseRows);
                        pointsSinceCheckpoint = maybe_checkpoint( ...
                            'COARSE', checkpointEveryPoints, checkpointFilePath, progressWorkbookPath, ...
                            runsCompletedGlobal, pointsSinceCheckpoint, coarsePointIndex, payload, ...
                            [], [], [], [], []);
                    end
                end
            end
        end
    end

    coarseResultsTable = array2table(coarseRows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm', ...
         'maxAbsTMOKE_base','alpha_peak_base_deg','TMOKE_at_peak_base', ...
         'alpha_peak_n2_deg','S_est_deg_per_RIU'});

% TOP-K selection using combined score (rank TMOKE + rank |S|)
    coarseSeedCandidates = selectTopK_tradeoff(coarseResultsTable, topKCoarse, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');

    save_checkpoint(checkpointFilePath, 'FINE', runsCompletedGlobal, 0, struct('coarseResultsTable', coarseResultsTable, 'coarseSeedCandidates', coarseSeedCandidates));
    write_progress_xlsx(progressWorkbookPath, 'coarse', coarseResultsTable);
end

%% ==================================================================
%                        FINE PLANNING (EXACT NOW)
% ------------------------------------------------------------------
% Build exact FINE counts by clamping search windows around the coarse
% seeds, so the run budget and ETA become exact before starting FINE.
% Also guards against exceeding MAX_RUNS before any fine sweep is run.
% ==================================================================
if isempty(coarseSeedCandidates) && resumeFromCheckpoint && resumeStageTag=="FINE"
    if isfield(checkpointData.payload,'coarseSeedCandidates')
        coarseSeedCandidates = checkpointData.payload.coarseSeedCandidates;
    else
        error('Checkpoint at FINE stage has no "coarseSeedCandidates". Rerun COARSE.');
    end
end

% Exact FINE point count (after clamping to search limits)
fineTotalPoints = 0;
for s = 1:height(coarseSeedCandidates)
    goldHeightList = max(min(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   - fineGoldHeightDelta)  : fineGoldHeightStep  : min(max(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   + fineGoldHeightDelta);
    ceyigHeightList = max(min(ceyigHeightGridNm), coarseSeedCandidates.h_ceyig_nm(s)- fineCeyigHeightDelta) : fineCeyigHeightStep : min(max(ceyigHeightGridNm), coarseSeedCandidates.h_ceyig_nm(s)+ fineCeyigHeightDelta);
    domainPeriodList = max(min(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)- fineDomainPeriodDelta) : fineDomainPeriodStep : min(max(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)+ fineDomainPeriodDelta);
    toothWidthList = max(min(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) - fineToothWidthDelta) : fineToothWidthStep : min(max(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) + fineToothWidthDelta);
    siliconHeightList = max(min(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    - fineSiliconHeightDelta)  : fineSiliconHeightStep  : min(max(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    + fineSiliconHeightDelta);
    fineTotalPoints = fineTotalPoints + numel(domainPeriodList)*numel(toothWidthList)*numel(siliconHeightList)*numel(ceyigHeightList)*numel(goldHeightList);
end
fineTotalRuns = runsPerSearchPoint * fineTotalPoints;

% Global total becomes exact from this point forward
globalRunTargetEstimate = coarseTotalRuns + fineTotalRuns + superTotalRuns + extraRunsFixed;
isTotalRunEstimateExact = true;

fprintf('\nSTAGE FINE - EXACT: %d runs (TOP-%d)\n', fineTotalRuns, height(coarseSeedCandidates));
fprintf('GLOBAL TOTAL (EXACT): %d runs\n\n', globalRunTargetEstimate);
if globalRunTargetEstimate > MAX_RUNS
    error('Planned total exceeds %d runs (%d). Adjust deltas/steps/TOPK.', MAX_RUNS, globalRunTargetEstimate);
end

%% ==================================================================
%                                 FINE
% ------------------------------------------------------------------
% For each coarse seed, sweep refined geometry windows and a narrower
% alpha range. Recompute TMOKE peak and sensitivity, log ETAs, and
% checkpoint. Select TOP-K seeds to feed SUPER.
% ==================================================================
fineResultsTable = [];
superSeedCandidates = [];

if resumeFromCheckpoint && any(strcmp(resumeStageTag, ["SUPER","VALID","FULL"]))
    if isfield(checkpointData.payload,'fineResultsTable'),      fineResultsTable = checkpointData.payload.fineResultsTable; end
    if isfield(checkpointData.payload,'superSeedCandidates'), superSeedCandidates = checkpointData.payload.superSeedCandidates; end
    fprintf('SKIP FINE -> restored fineResultsTable (rows=%d), superSeedCandidates (rows=%d)\n', ...
        size(fineResultsTable,1), height(superSeedCandidates));
else
    fineRows = [];
    stageRunsStart = runsCompletedGlobal;
    stageTotalRuns = fineTotalRuns;
    stageTimerStart = tic;

    finePointIndex = 0;
    if resumeFromCheckpoint && resumeStageTag=="FINE" && checkpointData.done_points > 0
        runsCompletedGlobal = checkpointData.runsCompletedGlobal;
        finePointIndex   = checkpointData.done_points;
        if isfield(checkpointData.payload,'fineRows'), fineRows = checkpointData.payload.fineRows; end
        if isfield(checkpointData.payload,'coarseSeedCandidates'),    coarseSeedCandidates    = checkpointData.payload.coarseSeedCandidates;    end
        fprintf('>>> FINE resume: skipping %d points already computed.\n', finePointIndex);
    end

    for s = 1:height(coarseSeedCandidates)
        % Alpha window center comes from the previous stage peak at baselineRefractiveIndex
        alphaWindowCenterDeg = coarseSeedCandidates.alpha_peak_base_deg(s);

        % Clamp parameters to remain inside allowed ranges
        goldHeightList = max(min(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   - fineGoldHeightDelta)  : fineGoldHeightStep  : min(max(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   + fineGoldHeightDelta);
        ceyigHeightList = max(min(ceyigHeightGridNm), coarseSeedCandidates.h_ceyig_nm(s)- fineCeyigHeightDelta) : fineCeyigHeightStep : min(max(ceyigHeightGridNm), coarseSeedCandidates.h_ceyig_nm(s)+ fineCeyigHeightDelta);
        domainPeriodList = max(min(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)- fineDomainPeriodDelta) : fineDomainPeriodStep : min(max(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)+ fineDomainPeriodDelta);
        toothWidthList = max(min(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) - fineToothWidthDelta) : fineToothWidthStep : min(max(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) + fineToothWidthDelta);
        siliconHeightList = max(min(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    - fineSiliconHeightDelta)  : fineSiliconHeightStep  : min(max(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    + fineSiliconHeightDelta);

        % Alpha sweep window (slightly wider for sensitivity to avoid missing the peak)
        alphaStartDeg = max(0,  alphaWindowCenterDeg - fineAlphaHalfSpanSensitivity);
        alphaStopDeg  = min(89, alphaWindowCenterDeg + fineAlphaHalfSpanSensitivity);

        for domainPeriodIdx = 1:numel(domainPeriodList)
            setParamNm(model, PARAM_LDOM, domainPeriodList(domainPeriodIdx));
            for toothWidthIdx = 1:numel(toothWidthList)
                setParamNm(model, PARAM_LDEN, toothWidthList(toothWidthIdx));
                for siliconHeightIdx = 1:numel(siliconHeightList)
                    setParamNm(model, PARAM_HSI, siliconHeightList(siliconHeightIdx));
                    for ceyigHeightIdx = 1:numel(ceyigHeightList)
                        setParamNm(model, PARAM_HCEY, ceyigHeightList(ceyigHeightIdx));
                        for goldHeightIdx = 1:numel(goldHeightList)
                            setParamNm(model, PARAM_HAU, goldHeightList(goldHeightIdx));

                            finePointIndex = finePointIndex + 1;
                            if resumeFromCheckpoint && resumeStageTag=="FINE" && finePointIndex <= checkpointData.done_points
                                continue;
                            end

                            % n = fastRefractiveIndexSamples(1)
                            setParamScalar(model, PARAM_N, fastRefractiveIndexSamples(1));
                            [alphaGridFastN1, reflectancePlusFastN1, reflectanceMinusFastN1, tmokeCurveFastN1] = solveAndGetRplusRminus( ...
                                model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                                alphaStartDeg, alphaFineStep, alphaStopDeg, ...
                                M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                            [tmokePeakMagnitudeFastN1, tmokePeakIndexFastN1] = max(abs(tmokeCurveFastN1));
                            alphaAtPeakFastN1 = alphaGridFastN1(tmokePeakIndexFastN1);
                            tmokeAtPeakFastN1 = tmokeCurveFastN1(tmokePeakIndexFastN1);

                            if PLOT_LIVE
                                updateLivePlot('FINE', domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), ceyigHeightList(ceyigHeightIdx), goldHeightList(goldHeightIdx), ...
                                    alphaGridFastN1, tmokeCurveFastN1, alphaAtPeakFastN1, tmokeAtPeakFastN1, reflectancePlusFastN1, reflectanceMinusFastN1);
                            end

                            % n = fastRefractiveIndexSamples(2) (for sensitivityEstimateFast calculation)
                            setParamScalar(model, PARAM_N, fastRefractiveIndexSamples(2));
                            [alphaGridFastN2, ~, ~, tmokeCurveFastN2] = solveAndGetRplusRminus( ...
                                model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                                alphaStartDeg, alphaFineStep, alphaStopDeg, ...
                                M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                            [~, tmokePeakIndexFastN2] = max(abs(tmokeCurveFastN2));
                            alphaAtPeakFastN2 = alphaGridFastN2(tmokePeakIndexFastN2);

                            sensitivityEstimateFast = (alphaAtPeakFastN2 - alphaAtPeakFastN1) / (fastRefractiveIndexSamples(2) - fastRefractiveIndexSamples(1));

                            fineRows = [fineRows; ...
                                domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), ceyigHeightList(ceyigHeightIdx), goldHeightList(goldHeightIdx), ...
                                tmokePeakMagnitudeFastN1, alphaAtPeakFastN1, tmokeAtPeakFastN1, alphaAtPeakFastN2, sensitivityEstimateFast]; %#ok<AGROW>

                            runsCompletedGlobal = runsCompletedGlobal + runsPerSearchPoint;

                            stageRunsCompleted = runsCompletedGlobal - stageRunsStart;
                            [stageCompletionFraction, stageEtaSeconds]   = frac_eta(stageRunsCompleted, stageTotalRuns, toc(stageTimerStart));
                            [globalCompletionFraction, globalEtaSeconds] = frac_eta(runsCompletedGlobal, globalRunTargetEstimate, toc(globalTimerStart));
                            logline(['FINE   | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                     ' | |TM|=%.5f @ alpha=%.3f deg | sensitivityEstimateFast=%.4f' ...
                                     ' [Stage %5.1f%% | ETA %s | Global %5.1f%% | ETA %s]\n'], ...
                                domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), ceyigHeightList(ceyigHeightIdx), goldHeightList(goldHeightIdx), ...
                                tmokePeakMagnitudeFastN1, alphaAtPeakFastN1, sensitivityEstimateFast, ...
                                100*stageCompletionFraction, fmt_time_long(stageEtaSeconds), 100*globalCompletionFraction, fmt_time_long(globalEtaSeconds));

                            pointsSinceCheckpoint = pointsSinceCheckpoint + 1;
                            payload = struct('fineRows',fineRows, 'coarseSeedCandidates', coarseSeedCandidates);
                            pointsSinceCheckpoint = maybe_checkpoint( ...
                                'FINE', checkpointEveryPoints, checkpointFilePath, progressWorkbookPath, ...
                                runsCompletedGlobal, pointsSinceCheckpoint, finePointIndex, payload, ...
                                coarseResultsTable, [], [], [], []);
                        end
                    end
                end
            end
        end
    end

    fineResultsTable = array2table(fineRows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm', ...
         'maxAbsTMOKE_base','alpha_peak_base_deg','TMOKE_at_peak_base', ...
         'alpha_peak_n2_deg','S_est_deg_per_RIU'});

    superSeedCandidates = selectTopK_tradeoff(fineResultsTable, topKFine, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');

    save_checkpoint(checkpointFilePath, 'SUPER', runsCompletedGlobal, 0, struct( ...
        'coarseResultsTable', coarseResultsTable, 'fineResultsTable', fineResultsTable, 'superSeedCandidates', superSeedCandidates));
    write_progress_xlsx(progressWorkbookPath, 'fine', fineResultsTable);
end

%% ==================================================================
%                     SUPERFINE PLANNING (EXACT)
% ------------------------------------------------------------------
% Compute exact SUPER run counts around each fine seed before launching
% the final refinement stage, keeping the run budget/ETA precise.
% ==================================================================
if isempty(superSeedCandidates) && resumeFromCheckpoint && resumeStageTag=="SUPER"
    if isfield(checkpointData.payload,'superSeedCandidates')
        superSeedCandidates = checkpointData.payload.superSeedCandidates;
    else
        error('Checkpoint at SUPER stage has no "superSeedCandidates". Rerun FINE.');
    end
end

superTotalPoints = 0;
for s = 1:height(superSeedCandidates)
    goldHeightList = (superSeedCandidates.h_au_nm(s)    - superGoldHeightDelta)  : superGoldHeightStep  : (superSeedCandidates.h_au_nm(s)    + superGoldHeightDelta);
    ceyigHeightList = (superSeedCandidates.h_ceyig_nm(s) - superCeyigHeightDelta) : superCeyigHeightStep : (superSeedCandidates.h_ceyig_nm(s) + superCeyigHeightDelta);
    domainPeriodList = (superSeedCandidates.L_domain_nm(s)- superDomainPeriodDelta) : superDomainPeriodStep : (superSeedCandidates.L_domain_nm(s)+ superDomainPeriodDelta);
    toothWidthList = (superSeedCandidates.l_dente_nm(s) - superToothWidthDelta) : superToothWidthStep : (superSeedCandidates.l_dente_nm(s) + superToothWidthDelta);
    siliconHeightList = (superSeedCandidates.h_si_nm(s)    - superSiliconHeightDelta)  : superSiliconHeightStep  : (superSeedCandidates.h_si_nm(s)    + superSiliconHeightDelta);
    superTotalPoints = superTotalPoints + numel(domainPeriodList)*numel(toothWidthList)*numel(siliconHeightList)*numel(ceyigHeightList)*numel(goldHeightList);
end
superTotalRuns = runsPerSearchPoint * superTotalPoints;

fprintf('\nSTAGE SUPER - EXACT: %d runs (TOP-%d)\n', superTotalRuns, height(superSeedCandidates));
fprintf('GLOBAL TOTAL (EXACT): %d runs\n\n', globalRunTargetEstimate);

%% ==================================================================
%                             SUPERFINE
% ------------------------------------------------------------------
% Final refinement around the fine seeds with tighter alpha steps.
% Recompute TMOKE peak and sensitivity, log ETAs, checkpoint, and pick:
%   - bestTradeoffCandidate (combined score),
%   - bestTmokeCandidate (max |TMOKE|),
%   - bestSensitivityCandidate (max |S_est|).
% ==================================================================
superResultsTable = [];
bestTradeoffCandidate = table(); bestTmokeCandidate = table(); bestSensitivityCandidate = table();

if resumeFromCheckpoint && any(strcmp(resumeStageTag, ["VALID","FULL"]))
    if isfield(checkpointData.payload,'superResultsTable'), superResultsTable = checkpointData.payload.superResultsTable; end
    if isfield(checkpointData.payload,'bestTradeoffCandidate'), bestTradeoffCandidate = checkpointData.payload.bestTradeoffCandidate; end
    if isfield(checkpointData.payload,'bestTmokeCandidate'),    bestTmokeCandidate    = checkpointData.payload.bestTmokeCandidate;    end
    if isfield(checkpointData.payload,'bestSensitivityCandidate'),        bestSensitivityCandidate        = checkpointData.payload.bestSensitivityCandidate;        end
    fprintf('SKIP SUPER -> restored superResultsTable (rows=%d) and best selections.\n', size(superResultsTable,1));
else
    superRows = [];
    stageRunsStart = runsCompletedGlobal;
    stageTotalRuns = superTotalRuns;
    stageTimerStart = tic;

    superPointIndex = 0;
    if resumeFromCheckpoint && resumeStageTag=="SUPER" && checkpointData.done_points > 0
        runsCompletedGlobal = checkpointData.runsCompletedGlobal;
        superPointIndex  = checkpointData.done_points;
        if isfield(checkpointData.payload,'superRows'),  superRows  = checkpointData.payload.superRows;  end
        if isfield(checkpointData.payload,'superSeedCandidates'), superSeedCandidates = checkpointData.payload.superSeedCandidates; end
        fprintf('>>> SUPER resume: skipping %d points already computed.\n', superPointIndex);
    end

    for s = 1:height(superSeedCandidates)
        % Alpha window center comes from the previous stage peak at baselineRefractiveIndex
        alphaWindowCenterDeg = superSeedCandidates.alpha_peak_base_deg(s);

        goldHeightList = (superSeedCandidates.h_au_nm(s)    - superGoldHeightDelta)  : superGoldHeightStep  : (superSeedCandidates.h_au_nm(s)    + superGoldHeightDelta);
        ceyigHeightList = (superSeedCandidates.h_ceyig_nm(s) - superCeyigHeightDelta) : superCeyigHeightStep : (superSeedCandidates.h_ceyig_nm(s) + superCeyigHeightDelta);
        domainPeriodList = (superSeedCandidates.L_domain_nm(s)- superDomainPeriodDelta) : superDomainPeriodStep : (superSeedCandidates.L_domain_nm(s)+ superDomainPeriodDelta);
        toothWidthList = (superSeedCandidates.l_dente_nm(s) - superToothWidthDelta) : superToothWidthStep : (superSeedCandidates.l_dente_nm(s) + superToothWidthDelta);
        siliconHeightList = (superSeedCandidates.h_si_nm(s)    - superSiliconHeightDelta)  : superSiliconHeightStep  : (superSeedCandidates.h_si_nm(s)    + superSiliconHeightDelta);

        % Alpha sweep window (slightly wider for sensitivity to avoid missing the peak)
        alphaStartDeg = max(0,  alphaWindowCenterDeg - superAlphaHalfSpanSensitivity);
        alphaStopDeg  = min(89, alphaWindowCenterDeg + superAlphaHalfSpanSensitivity);

        for domainPeriodIdx = 1:numel(domainPeriodList)
            setParamNm(model, PARAM_LDOM, domainPeriodList(domainPeriodIdx));
            for toothWidthIdx = 1:numel(toothWidthList)
                setParamNm(model, PARAM_LDEN, toothWidthList(toothWidthIdx));
                for siliconHeightIdx = 1:numel(siliconHeightList)
                    setParamNm(model, PARAM_HSI, siliconHeightList(siliconHeightIdx));
                    for ceyigHeightIdx = 1:numel(ceyigHeightList)
                        setParamNm(model, PARAM_HCEY, ceyigHeightList(ceyigHeightIdx));
                        for goldHeightIdx = 1:numel(goldHeightList)
                            setParamNm(model, PARAM_HAU, goldHeightList(goldHeightIdx));

                            superPointIndex = superPointIndex + 1;
                            if resumeFromCheckpoint && resumeStageTag=="SUPER" && superPointIndex <= checkpointData.done_points
                                continue;
                            end

                            setParamScalar(model, PARAM_N, fastRefractiveIndexSamples(1));
                            [alphaGridFastN1, reflectancePlusFastN1, reflectanceMinusFastN1, tmokeCurveFastN1] = solveAndGetRplusRminus( ...
                                model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                                alphaStartDeg, alphaSuperStep, alphaStopDeg, ...
                                M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                            [tmokePeakMagnitudeFastN1, tmokePeakIndexFastN1] = max(abs(tmokeCurveFastN1));
                            alphaAtPeakFastN1 = alphaGridFastN1(tmokePeakIndexFastN1);
                            tmokeAtPeakFastN1 = tmokeCurveFastN1(tmokePeakIndexFastN1);

                            if PLOT_LIVE
                                updateLivePlot('SUPER', domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), ceyigHeightList(ceyigHeightIdx), goldHeightList(goldHeightIdx), ...
                                    alphaGridFastN1, tmokeCurveFastN1, alphaAtPeakFastN1, tmokeAtPeakFastN1, reflectancePlusFastN1, reflectanceMinusFastN1);
                            end

                            setParamScalar(model, PARAM_N, fastRefractiveIndexSamples(2));
                            [alphaGridFastN2, ~, ~, tmokeCurveFastN2] = solveAndGetRplusRminus( ...
                                model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                                alphaStartDeg, alphaSuperStep, alphaStopDeg, ...
                                M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                            [~, tmokePeakIndexFastN2] = max(abs(tmokeCurveFastN2));
                            alphaAtPeakFastN2 = alphaGridFastN2(tmokePeakIndexFastN2);

                            sensitivityEstimateFast = (alphaAtPeakFastN2 - alphaAtPeakFastN1) / (fastRefractiveIndexSamples(2) - fastRefractiveIndexSamples(1));

                            superRows = [superRows; ...
                                domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), ceyigHeightList(ceyigHeightIdx), goldHeightList(goldHeightIdx), ...
                                tmokePeakMagnitudeFastN1, alphaAtPeakFastN1, tmokeAtPeakFastN1, alphaAtPeakFastN2, sensitivityEstimateFast]; %#ok<AGROW>

                            runsCompletedGlobal = runsCompletedGlobal + runsPerSearchPoint;

                            stageRunsCompleted = runsCompletedGlobal - stageRunsStart;
                            [stageCompletionFraction, stageEtaSeconds]   = frac_eta(stageRunsCompleted, stageTotalRuns, toc(stageTimerStart));
                            [globalCompletionFraction, globalEtaSeconds] = frac_eta(runsCompletedGlobal, globalRunTargetEstimate, toc(globalTimerStart));
                            logline(['SUPER  | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                     ' | |TM|=%.5f @ alpha=%.4f deg | sensitivityEstimateFast=%.4f' ...
                                     ' [Stage %5.1f%% | ETA %s | Global %5.1f%% | ETA %s]\n'], ...
                                domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), ceyigHeightList(ceyigHeightIdx), goldHeightList(goldHeightIdx), ...
                                tmokePeakMagnitudeFastN1, alphaAtPeakFastN1, sensitivityEstimateFast, ...
                                100*stageCompletionFraction, fmt_time_long(stageEtaSeconds), 100*globalCompletionFraction, fmt_time_long(globalEtaSeconds));

                            pointsSinceCheckpoint = pointsSinceCheckpoint + 1;
                            payload = struct('superRows',superRows, 'superSeedCandidates', superSeedCandidates);
                            pointsSinceCheckpoint = maybe_checkpoint( ...
                                'SUPER', checkpointEveryPoints, checkpointFilePath, progressWorkbookPath, ...
                                runsCompletedGlobal, pointsSinceCheckpoint, superPointIndex, payload, ...
                                coarseResultsTable, fineResultsTable, [], [], []);
                        end
                    end
                end
            end
        end
    end

    superResultsTable = array2table(superRows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm', ...
         'maxAbsTMOKE_base','alpha_peak_base_deg','TMOKE_at_peak_base', ...
         'alpha_peak_n2_deg','S_est_deg_per_RIU'});

    % Escolhas finais dentro do SUPER:
    bestTradeoffCandidate = selectTopK_tradeoff(superResultsTable, 1, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');
    bestTmokeCandidate    = selectTopK_single(superResultsTable, 1, 'maxAbsTMOKE_base');      % max |TMOKE|
    bestSensitivityCandidate        = selectTopK_single_abs(superResultsTable, 1, 'S_est_deg_per_RIU'); % max |S|

    fprintf('\n===== BEST (SUPERFINE) =====\n');
    fprintf('Trade-off (rankTM+rankS):\n'); disp(bestTradeoffCandidate);
    fprintf('Best TMOKE only:\n');         disp(bestTmokeCandidate);
    fprintf('Best |S| only:\n');          disp(bestSensitivityCandidate);

    save_checkpoint(checkpointFilePath, 'VALID', runsCompletedGlobal, 0, struct( ...
        'coarseResultsTable', coarseResultsTable, 'fineResultsTable', fineResultsTable, 'superResultsTable', superResultsTable, ...
        'bestTradeoffCandidate', bestTradeoffCandidate, 'bestTmokeCandidate', bestTmokeCandidate, 'bestSensitivityCandidate', bestSensitivityCandidate));
    write_progress_xlsx(progressWorkbookPath, 'super', superResultsTable);
end

%% ==================================================================
%                 VALID: dense curves + sensitivityDense (3 n values)
% ------------------------------------------------------------------
% Fix geometry to bestTradeoffCandidate and:
%   - sweep dense alpha for each validationRefractiveIndexList value,
%   - find alpha_peak for each n,
%   - fit alpha_peak vs n to compute sensitivityDense,
%   - store dense tables for later export/plots.
% ==================================================================
% From here, geometry is fixed to bestTradeoffCandidate (best compromise).
if isempty(bestTradeoffCandidate) && resumeFromCheckpoint && any(strcmp(resumeStageTag, ["VALID","FULL"]))
    bestTradeoffCandidate = checkpointData.payload.bestTradeoffCandidate;
end

Ldom_best = bestTradeoffCandidate.L_domain_nm(1);
Lden_best = bestTradeoffCandidate.l_dente_nm(1);
hsi_best  = bestTradeoffCandidate.h_si_nm(1);
hcey_best = bestTradeoffCandidate.h_ceyig_nm(1);
hau_best  = bestTradeoffCandidate.h_au_nm(1);

setParamNm(model, PARAM_LDOM, Ldom_best);
setParamNm(model, PARAM_LDEN, Lden_best);
setParamNm(model, PARAM_HSI,  hsi_best);
setParamNm(model, PARAM_HCEY, hcey_best);
setParamNm(model, PARAM_HAU,  hau_best);

alphaPeakDegreesByN = zeros(size(validationRefractiveIndexList));
tmokeMaxAbsByN   = zeros(size(validationRefractiveIndexList));
tmokeCurvesByN      = cell(size(validationRefractiveIndexList));
alphaGridsByN   = cell(size(validationRefractiveIndexList));
reflectancePlusByN      = cell(size(validationRefractiveIndexList));
reflectanceMinusByN      = cell(size(validationRefractiveIndexList));

for i = 1:numel(validationRefractiveIndexList)
    setParamScalar(model, PARAM_N, validationRefractiveIndexList(i));

    [a, Rp, Rm, TM] = solveAndGetRplusRminus( ...
        model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
        0, alphaDenseStep, 89, ...
        M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

    [tmmax, k] = max(abs(TM));
    alphaPeakDegreesByN(i) = a(k);
    tmokeMaxAbsByN(i)   = tmmax;

    alphaGridsByN{i} = a;
    tmokeCurvesByN{i}    = TM;
    reflectancePlusByN{i}    = Rp;
    reflectanceMinusByN{i}    = Rm;

    runsCompletedGlobal = runsCompletedGlobal + 2;
end

alphaVsNLinearFit = polyfit(validationRefractiveIndexList, alphaPeakDegreesByN, 1);
sensitivityDense = alphaVsNLinearFit(1);

fprintf('\n===== VALID (bestTradeoffCandidate) =====\n');
fprintf('Best geom: Ldom=%g | Lden=%g | hsi=%g | hcey=%g | hau=%g\n', ...
    Ldom_best, Lden_best, hsi_best, hcey_best, hau_best);
fprintf('alpha_peak(n) = ['); fprintf(' %.4f', alphaPeakDegreesByN); fprintf(' ]\n');
fprintf('sensitivityDense ~= %.6f deg/RIU (linear fit)\n', sensitivityDense);

% Define alphaBestDeg / tmokeBestValue as the peak at baselineRefractiveIndex (if inside validationRefractiveIndexList)
[~, idxBase] = min(abs(validationRefractiveIndexList - baselineRefractiveIndex));
alphaBestDeg = alphaPeakDegreesByN(idxBase);

% Recompute tmokeBestValue at baselineRefractiveIndex using cached data
[tmokeBaselineMax, tmokeBaselineIndex] = max(abs(tmokeCurvesByN{idxBase}));
tmokeBestValue = tmokeCurvesByN{idxBase}(tmokeBaselineIndex);

% Consolidate dense results for every alpha point at each n (for export/plots)
bestDenseTable = table();
for i = 1:numel(validationRefractiveIndexList)
    ncol = repmat(validationRefractiveIndexList(i), numel(alphaGridsByN{i}), 1);
    Ttmp = table( ...
        repmat(Ldom_best,numel(alphaGridsByN{i}),1), repmat(Lden_best,numel(alphaGridsByN{i}),1), repmat(hsi_best,numel(alphaGridsByN{i}),1), ...
        repmat(hcey_best,numel(alphaGridsByN{i}),1), repmat(hau_best,numel(alphaGridsByN{i}),1), ...
        ncol, alphaGridsByN{i}(:), reflectancePlusByN{i}(:), reflectanceMinusByN{i}(:), tmokeCurvesByN{i}(:), ...
        'VariableNames', {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm', ...
                          'n','alpha_deg','Rplus','Rminus','TMOKE'});
    bestDenseTable = [bestDenseTable; Ttmp]; %#ok<AGROW>
end

save_checkpoint(checkpointFilePath, 'FULL', runsCompletedGlobal, 0, struct( ...
    'coarseResultsTable', coarseResultsTable, 'fineResultsTable', fineResultsTable, 'superResultsTable', superResultsTable, ...
    'bestTradeoffCandidate', bestTradeoffCandidate, 'bestTmokeCandidate', bestTmokeCandidate, 'bestSensitivityCandidate', bestSensitivityCandidate, ...
    'bestDenseTable', bestDenseTable, 'sensitivityDense', sensitivityDense, ...
    'alphaBestDeg', alphaBestDeg, 'tmokeBestValue', tmokeBestValue));
write_progress_xlsx(progressWorkbookPath, 'dense', bestDenseTable);

%% ==================================================================
%                           Snapshot (optional)
% ==================================================================
% Capture the bestTradeoffCandidate at alphaBestDeg and baselineRefractiveIndex for inspection
if SAVE_SNAPSHOT
    setParamScalar(model, PARAM_N, baselineRefractiveIndex);

    setAlphaMSweep(model, STUDY_TAG, ALPHA_NAME, alphaBestDeg, 0.01, alphaBestDeg, MSIGN_NAME, sprintf('%s %s', M_PLUS, M_MINUS));
    model.study(STUDY_TAG).run;
    refreshDerivedValues(model);

    runsCompletedGlobal = runsCompletedGlobal + 2;

    snapshotTimestamp = datestr(now,'yyyymmdd_HHMMSS');
    snapshotFilePath = fullfile(projectRootDir, sprintf( ...
        'snapshot_bestTradeoff_Ldom%4.0f_Lden%4.0f_hsi%3.0f_hcey%.3f_hau%.3f_n%.3f_alpha%.4f_%s.mph', ...
        Ldom_best, Lden_best, hsi_best, hcey_best, hau_best, baselineRefractiveIndex, alphaBestDeg, snapshotTimestamp));
    mphsave(model, snapshotFilePath);
    logline('Snapshot saved: %s\n', snapshotFilePath);
end

%% ==================================================================
%                     FULL: final baselineRefractiveIndex curve
% ==================================================================
% FULL keeps the standard single-n (baselineRefractiveIndex) curve for a
% clean export and baseline comparison with the VALID dense curves.
setParamScalar(model, PARAM_N, baselineRefractiveIndex);

[alphaFull, reflectancePlusFull, reflectanceMinusFull, tmokeFull] = solveAndGetRplusRminus( ...
    model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
    0, alphaFullStep, 89, ...
    M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

runsCompletedGlobal = runsCompletedGlobal + 2;

if PLOT_LIVE
    updateLivePlot('FULL', Ldom_best, Lden_best, hsi_best, hcey_best, hau_best, ...
        alphaFull, tmokeFull, alphaBestDeg, tmokeBestValue, reflectancePlusFull, reflectanceMinusFull);
end

bestFullTable = table( ...
    repmat(Ldom_best,numel(alphaFull),1), repmat(Lden_best,numel(alphaFull),1), repmat(hsi_best,numel(alphaFull),1), ...
    repmat(hcey_best,numel(alphaFull),1), repmat(hau_best,numel(alphaFull),1), ...
    repmat(baselineRefractiveIndex,numel(alphaFull),1), ...
    alphaFull(:), reflectancePlusFull(:), reflectanceMinusFull(:), tmokeFull(:), ...
    'VariableNames', {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm', ...
                      'n','alpha_deg','Rplus','Rminus','TMOKE'});

%% ==================================================================
%                               CSV / XLSX
% ==================================================================
% Persist coarse/fine/super candidate tables and the bestTradeoffCandidate
% dense/baseline curves for downstream analysis in Python/Excel.
if ~isempty(coarseResultsTable), writetable(coarseResultsTable, fullfile(projectRootDir,'tmoke_sens_5D_coarse.csv')); end
if ~isempty(fineResultsTable),   writetable(fineResultsTable,   fullfile(projectRootDir,'tmoke_sens_5D_fine.csv')); end
if ~isempty(superResultsTable),  writetable(superResultsTable,  fullfile(projectRootDir,'tmoke_sens_5D_super.csv')); end

writetable(bestDenseTable, fullfile(projectRootDir,'tmoke_sens_bestTradeoff_dense_ALLn.csv'));
writetable(bestFullTable,  fullfile(projectRootDir,'tmoke_sens_bestTradeoff_full_baseline.csv'));

write_progress_xlsx(progressWorkbookPath, 'full', bestFullTable);

% Delete checkpoint at the end of a complete execution
if isfile(checkpointFilePath), delete(checkpointFilePath); end

%% ==================================================================
%                               FINAL PLOTS
% ------------------------------------------------------------------
% Plot the validated TMOKE curves, alpha_peak vs n fit, |TMOKE|max vs n,
% and a scatter map of all candidates vs the chosen bestTradeoffCandidate.
% Saved to figureOutputDirectory when SAVE_FIGS is true.
% ==================================================================
if MAKE_PLOTS
    % (1) TMOKE(alpha) for each n (VALID)
    figure('Name','(1) TMOKE(alpha) for each n (VALID)','NumberTitle','off','Color','w');
    hold on; grid on;
    for i = 1:numel(validationRefractiveIndexList)
        plot(alphaGridsByN{i}, tmokeCurvesByN{i}, 'LineWidth', 1.2, 'DisplayName', sprintf('n=%.2f', validationRefractiveIndexList(i)));
    end
    xlabel('\alpha [deg]'); ylabel('TMOKE(\alpha)');
    title(sprintf('Best trade-off | TMOKE(\\alpha) vs \\alpha (S_{dense}=%.4f deg/RIU)', sensitivityDense));
    legend('Location','best');

    % (2) alpha_peak vs n + fit
    figure('Name','(2) alpha_peak vs n (VALID)','NumberTitle','off','Color','w');
    hold on; grid on;
    plot(validationRefractiveIndexList, alphaPeakDegreesByN, 'o', 'LineWidth', 1.5, 'DisplayName','data');
    nfit = linspace(min(validationRefractiveIndexList), max(validationRefractiveIndexList), 100);
    plot(nfit, polyval(alphaVsNLinearFit,nfit), '-', 'LineWidth', 1.2, 'DisplayName', sprintf('fit (S=%.4f)', sensitivityDense));
    xlabel('n'); ylabel('\alpha_{peak} [deg]');
    title('alpha_{peak} as a function of n');
    legend('Location','best');

    % (3) |TMOKE|max vs n
    figure('Name','(3) |TMOKE|max vs n (VALID)','NumberTitle','off','Color','w');
    grid on; hold on;
    plot(validationRefractiveIndexList, tmokeMaxAbsByN, 'o-', 'LineWidth', 1.5);
    xlabel('n'); ylabel('|TMOKE|_{max}');
    title('|TMOKE|_{max} vs n');

    % (4) Scatter (trade-off) in the (|TMOKE|, |sensitivityEstimateFast|) plane
    allCandidateResults = coarseResultsTable;
    if ~isempty(fineResultsTable),  allCandidateResults = [allCandidateResults; fineResultsTable]; end %#ok<AGROW>
    if ~isempty(superResultsTable), allCandidateResults = [allCandidateResults; superResultsTable]; end %#ok<AGROW>

    figure('Name','(4) Trade-off map: |TMOKE| vs |sensitivityEstimateFast|','NumberTitle','off','Color','w');
    grid on; hold on;
    scatter(abs(allCandidateResults.maxAbsTMOKE_base), abs(allCandidateResults.S_est_deg_per_RIU), 25, 'filled');
    scatter(abs(bestTradeoffCandidate.maxAbsTMOKE_base), abs(bestTradeoffCandidate.S_est_deg_per_RIU), 120, 'filled');
    xlabel('|TMOKE|_{max} (baseline)'); ylabel('|S_{est}| [deg/RIU]');
    title('Candidates (COARSE+FINE+SUPER) with bestTradeoffCandidate highlighted');
end

%% ==================================================================
%                           SAVE ALL FIGURES
% ==================================================================
if SAVE_FIGS
    saveAllOpenFigures(figureOutputDirectory, "tmoke_sens_5d", FIG_FORMATS);
    fprintf('Figures saved in: %s\n', figureOutputDirectory);
end

%% ==================================================================
%                               SUMMARY
% ==================================================================
% Print the final selections and validated sensitivity, plus run counts
% and wall-clock time. The detailed log remains in the runlog_* file.
elapsed_total = toc(globalTimerStart);
fprintf('\n===== SUMMARY =====\n');
fprintf('Best TRADE-OFF (SUPER):\n'); disp(bestTradeoffCandidate);
fprintf('Best TMOKE only (SUPER):\n'); disp(bestTmokeCandidate);
fprintf('Best |sensitivityEstimateFast| only (SUPER):\n'); disp(bestSensitivityCandidate);
fprintf('VALID @ bestTradeoffCandidate: sensitivityDense ~= %.6f deg/RIU\n', sensitivityDense);
fprintf('Runs done: %d | Elapsed: %s\n', runsCompletedGlobal, fmt_time_long(elapsed_total));

diary off;

%% =====================================================================
%                           LOCAL FUNCTIONS
% =====================================================================

function setParamNm(mdl, name, val_nm)
    % Set a COMSOL parameter that carries nanometers as its unit
    mdl.param.set(name, sprintf('%.12g[nm]', val_nm));
end

function setParamScalar(mdl, name, val)
    % Set a unitless COMSOL parameter (e.g., refractive index)
    mdl.param.set(name, sprintf('%.12g', val));
end

function Tsel = selectTopK_tradeoff(T, K, colTM, colS)
% Select TOP-K rows using the combined score: rank(|TM|) + rank(|S|)
% - rank 1 = best
    tm = abs(T.(colTM));
    ss = abs(T.(colS));

    rTM = tiedrank(-tm);
    rS  = tiedrank(-ss);

    score = rTM + rS;
    T.score_tradeoff = score; %#ok<AGROW>

    [~, ord] = sort(score, 'ascend');
    Tsel = T(ord(1:min(K,height(T))), :);
end

function Tsel = selectTopK_single(T, K, col)
% Select TOP-K rows by the largest value in a column (no abs)
    [~, ord] = sort(T.(col), 'descend');
    Tsel = T(ord(1:min(K,height(T))), :);
end

function Tsel = selectTopK_single_abs(T, K, col)
% Select TOP-K rows by the largest absolute value in a column
    [~, ord] = sort(abs(T.(col)), 'descend');
    Tsel = T(ord(1:min(K,height(T))), :);
end

function [alpha_deg, Rplus, Rminus, TMOKE] = solveAndGetRplusRminus( ...
    mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, ...
    mPlusStr, mMinusStr, ttagPlus, ttagMinus)

    ensureTable(mdl, ttagPlus);
    ensureTable(mdl, ttagMinus);

    Npts = 1 + floor((aStopDeg - aStartDeg)/aStepDeg + 1e-9);

    % m = +1
    redirectAllNumericalsToTable(mdl, ttagPlus);
    clearTable(mdl, ttagPlus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mPlusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha1_deg, R1] = readAlphaAndRFromNamedTable(mdl, ttagPlus, Npts);

    % m = -1
    redirectAllNumericalsToTable(mdl, ttagMinus);
    clearTable(mdl, ttagMinus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mMinusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha2_deg, R2] = readAlphaAndRFromNamedTable(mdl, ttagMinus, Npts);

    assert(numel(alpha1_deg)==Npts && numel(alpha2_deg)==Npts, 'Unexpected alpha sweep length.');
    assert(max(abs(alpha1_deg - alpha2_deg)) < 1e-9, 'alpha grids differ between m=+1 and m=-1.');

    alpha_deg = alpha1_deg;
    Rplus  = R1;
    Rminus = R2;

    denom = Rplus + Rminus;
    denom(abs(denom) < 1e-9) = 1e-9;
    TMOKE = (Rplus - Rminus) ./ denom;
end

function setTwoParamSweep(mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mStr)
    ptag = getParametricTag(mdl, studyTag);
    mdl.study(studyTag).feature(ptag).set('pname', {alphaName, mName});
    mdl.study(studyTag).feature(ptag).set('punit', {'deg','1'});
    alist = sprintf('range(%.12g[deg], %.12g[deg], %.12g[deg])', aStartDeg, aStepDeg, aStopDeg);
    mdl.study(studyTag).feature(ptag).set('plistarr', {alist, mStr});
end

function ptag = getParametricTag(mdl, studyTag)
    ptag = 'param';
    try
        fts = cell(mdl.study(studyTag).feature().tags());
        for i = 1:numel(fts)
            typ = char(mdl.study(studyTag).feature(fts{i}).featureType());
            if contains(lower(typ), 'param'), ptag = fts{i}; break; end
        end
    catch
    end
end

function refreshDerivedValues(mdl)
    try
        ntags = cell(mdl.result().numerical().tags());
        for k = 1:numel(ntags)
            mdl.result().numerical(ntags{k}).setResult;
        end
    catch
    end
end

function ensureTable(mdl, ttag)
    try
        mdl.result().table(ttag);
    catch
        try mdl.result().table().create(ttag, 'Table'); catch, end
    end
end

function clearTable(mdl, ttag)
    try mdl.result().table(ttag).clearTableData; catch, end
end

function redirectAllNumericalsToTable(mdl, ttag)
    try
        ntags = cell(mdl.result().numerical().tags());
        for k = 1:numel(ntags)
            try mdl.result().numerical(ntags{k}).set('table', ttag); catch, end
        end
    catch
    end
end

function [alpha_deg, Rcol] = readAlphaAndRFromNamedTable(mdl, ttag, Npts)
    S = mphtable(mdl, ttag);

    heads = [];
    if isfield(S,'colhead') && ~isempty(S.colhead), heads = string(S.colhead); end
    if isempty(heads) && isfield(S,'head') && ~isempty(S.head), heads = string(S.head); end
    if isempty(heads) && isfield(S,'header') && ~isempty(S.header), heads = string(S.header); end
    if isempty(heads) && isfield(S,'colnames') && ~isempty(S.colnames), heads = string(S.colnames); end

    if isempty(heads)
        if ~isfield(S,'data') || isempty(S.data) || size(S.data,2) < 2
            error('Table %s empty after run. Check Derived Values.', ttag);
        end
        a = S.data(:,1); R = S.data(:,2);
    else
        hlow = lower(heads);
        aIdx = find(contains(hlow,'alpha'), 1, 'first');
        rIdx = find(contains(hlow,'reflectance') | contains(hlow,'total reflectance') | contains(hlow,' total r'), 1, 'first');
        if isempty(aIdx) || isempty(rIdx)
            if ~isfield(S,'data') || isempty(S.data) || size(S.data,2) < 2
                error('Alpha/Reflectance not found in %s and not enough columns.', ttag);
            end
            aIdx = 1; rIdx = 2;
        end
        a = S.data(:, aIdx); R = S.data(:, rIdx);
    end

    assert(numel(a) >= Npts && numel(R) >= Npts, ...
        'Table %s has %d rows; expected >= %d.', ttag, numel(a), Npts);

    a = a(end-Npts+1:end);
    R = R(end-Npts+1:end);

    if ~isempty(heads)
        hlow = lower(heads);
        if any(contains(hlow, '(rad)')) || any(contains(hlow, '[rad]')), a = a * 180/pi; end
    end

    alpha_deg = a;
    Rcol      = R;
end

function s = fmt_time_long(sec)
    if ~isfinite(sec) || sec < 0, s = '...'; return; end
    days = floor(sec/86400);
    rem  = sec - 86400*days;
    hrs  = floor(rem/3600);
    rem  = rem - 3600*hrs;
    mins = floor(rem/60);
    secs = floor(rem - 60*mins);
    if days > 0
        s = sprintf('%dd %02d:%02d:%02d', days, hrs, mins, secs);
    else
        s = sprintf('%02d:%02d:%02d', hrs, mins, secs);
    end
end

function [frac, eta] = frac_eta(done, total, elapsed)
    if ~isfinite(total) || total <= 0
        frac = NaN; eta = NaN; return;
    end
    frac = max(0,min(1, double(done)/double(total)));
    if done <= 0 || ~isfinite(elapsed) || elapsed <= 0 || frac <= 0
        eta = NaN;
    else
        eta = elapsed/frac - elapsed;
    end
end

function s = fmt_eta_flag(eta, is_exact)
    if ~isfinite(eta), s = 'N/A'; else, s = fmt_time_long(eta); end
    if ~is_exact, s = [s ' (est.)']; end
end

function s = fmt_pct(frac)
    if ~isfinite(frac), s = 'N/A'; else, s = sprintf('%5.1f', 100*frac); end
end

function logline(varargin)
    fprintf(varargin{:});
    drawnow('limitrate');
end

function setAlphaMSweep(mdl, studyTag, alphaName, aStartDeg, aStepDeg, aStopDeg, mName, mListStr)
    ptag = getParametricTag(mdl, studyTag);
    mdl.study(studyTag).feature(ptag).set('pname', {alphaName, mName});
    mdl.study(studyTag).feature(ptag).set('punit', {'deg','1'});
    alist = sprintf('range(%.12g[deg], %.12g[deg], %.12g[deg])', aStartDeg, aStepDeg, aStopDeg);
    mdl.study(studyTag).feature(ptag).set('plistarr', {alist, mListStr});
end

function save_checkpoint(checkpointFilePath, stage, runsCompletedGlobal, done_points, payload)
    checkpointData.stage            = char(stage);
    checkpointData.runsCompletedGlobal = runsCompletedGlobal;
    checkpointData.done_points      = done_points;
    checkpointData.payload          = payload;
    tmp = [checkpointFilePath '.tmp'];
    save(tmp,'checkpointData','-v7.3');
    movefile(tmp, checkpointFilePath, 'f');
end

function write_progress_xlsx(progressWorkbookPath, varargin)
    for i=1:2:numel(varargin)
        sheet = varargin{i}; T = varargin{i+1};
        if ~isempty(T)
            try
                writetable(T, progressWorkbookPath, 'Sheet', sheet);
            catch
                warning('write_progress_xlsx: failed sheet %s', sheet);
            end
        end
    end
end

function points_since_ckpt_new = maybe_checkpoint( ...
        stage, every_points, checkpointFilePath, progressWorkbookPath, ...
        runsCompletedGlobal, pointsSinceCheckpoint, done_points, payload, ...
        coarseResultsTable, fineResultsTable, superResultsTable, bestDenseTable, bestFullTable)
    points_since_ckpt_new = pointsSinceCheckpoint;
    if pointsSinceCheckpoint >= every_points
        save_checkpoint(checkpointFilePath, stage, runsCompletedGlobal, done_points, payload);
        try
            write_progress_xlsx(progressWorkbookPath, ...
                'coarse', coarseResultsTable, 'fine', fineResultsTable, 'super', superResultsTable, ...
                'dense', bestDenseTable, 'full', bestFullTable);
        catch
        end
        points_since_ckpt_new = 0;
    end
end

function updateLivePlot(stage, Ldom, Lden, hsi, hcey, hau, ...
                        alpha_deg, TMOKE_vec, alpha_peak, tmoke_peak, ...
                        Rplus, Rminus)
    persistent fig ax hTM hPeak hRp hRm inited
    if isempty(inited) || ~isvalidHandle(fig) || ~isvalidHandle(ax) || ...
       isempty(hTM)   || ~isvalidHandle(hTM)   || ...
       isempty(hPeak) || ~isvalidHandle(hPeak) || ...
       isempty(hRp)   || ~isvalidHandle(hRp)   || ...
       isempty(hRm)   || ~isvalidHandle(hRm)

        fig = figure('Name','TMOKE vs \alpha (live)', ...
                     'NumberTitle','off','Tag','TMOKE_LIVE_FIG', 'Color','w');
        ax  = axes('Parent',fig,'Tag','TMOKE_LIVE_AX');
        grid(ax,'on'); hold(ax,'on');

        yyaxis(ax,'left');
        hTM   = plot(ax, nan, nan, '-', 'LineWidth', 1.2, 'DisplayName','TMOKE(\alpha)');
        hPeak = plot(ax, nan, nan, 'o', 'MarkerSize', 6, 'DisplayName','peak TMOKE');
        ylabel(ax,'TMOKE'); xlim(ax,[0 89]); xlabel(ax,'\alpha [deg]');

        yyaxis(ax,'right');
        hRp = plot(ax, nan, nan, '--', 'LineWidth', 1.0, 'DisplayName','R^+(\alpha)');
        hRm = plot(ax, nan, nan, ':',  'LineWidth', 1.0, 'DisplayName','R^-(\alpha)');
        ylabel(ax,'Reflectance');

        legend(ax,'Location','best');
        inited = true;
    end

    yyaxis(ax,'left');
    set(hTM,   'XData', alpha_deg, 'YData', TMOKE_vec);
    set(hPeak, 'XData', alpha_peak, 'YData', tmoke_peak);

    yyaxis(ax,'right');
    if nargin >= 12 && ~isempty(Rplus) && ~isempty(Rminus)
        set(hRp, 'XData', alpha_deg, 'YData', Rplus);
        set(hRm, 'XData', alpha_deg, 'YData', Rminus);
    else
        set(hRp, 'XData', nan, 'YData', nan);
        set(hRm, 'XData', nan, 'YData', nan);
    end

    title(ax, sprintf('%s | Ldom=%g | Lden=%g | hsi=%g | hcey=%.3f | hau=%.3f | |TM|=%.5f @ %.3f deg', ...
        stage, Ldom, Lden, hsi, hcey, hau, abs(tmoke_peak), alpha_peak));

    drawnow('limitrate');
end

function tf = isvalidHandle(h)
    tf = ~isempty(h) && isvalid(h);
end

function saveAllOpenFigures(outDir, prefix, formats)
    if ~exist(outDir,'dir'), mkdir(outDir); end
    figs = findall(0,'Type','figure');
    if isempty(figs), fprintf('>> No figures to save.\n'); return; end

    [~, idx] = sort([figs.Number]);
    figs = figs(idx);

    for i = 1:numel(figs)
        fig = figs(i);
        figName = string(get(fig,'Name'));
        figTag  = string(get(fig,'Tag'));
        if strlength(figName)==0
            if strlength(figTag)>0, figName = figTag; else, figName = "Figure_" + string(fig.Number); end
        end

        base = string(prefix) + "__" + string(sanitizeFilename(figName));
        try set(fig,'Color','w'); catch, end

        for k = 1:numel(formats)
            ext = lower(string(formats{k}));
            filePath = fullfile(outDir, base + "." + ext);
            try
                switch ext
                    case "png"
                        exportgraphics(fig, filePath, 'Resolution', 300);
                    case "pdf"
                        exportgraphics(fig, filePath, 'ContentType','vector');
                    otherwise
                        saveas(fig, filePath);
                end
            catch ME
                warning('Failed to save %s: %s', filePath, ME.message);
            end
        end
    end

    fprintf('>> %d figure(s) saved in: %s\n', numel(figs), outDir);
end

function s = sanitizeFilename(str)
    s = regexprep(char(str), '[^\w\d\-]+', '_');
    s = regexprep(s, '_+', '_');
    s = regexprep(s, '^_|_$', '');
end
