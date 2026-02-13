function runSenseAndTmoke(cfg)
% Main TMOKE + sensitivity runner using modular workflow stages.
import com.comsol.model.*;
import com.comsol.model.util.*;
import tmoke.comsol.*;
import tmoke.selection.*;
import tmoke.checkpoint.*;
import tmoke.plot.*;
import tmoke.util.*;
import tmoke.workflow.*;

if nargin < 1 || isempty(cfg)
    error('Config struct required. Use scripts/run_senseAndTmoke.m');
end

% Paths and folders
projectRootDir       = cfg.outputDir;
checkpointDirectory  = cfg.checkpointDir;
progressWorkbookPath = cfg.progressWorkbookPath;
checkpointFilePath   = fullfile(checkpointDirectory,'tmoke_sens_5d_checkpoint.mat');
runTimestamp = datestr(now,'yyyy-mm-dd_HHMMSS');
figureOutputDirectory = fullfile(cfg.plotsRootDir, ['tmoke_sens_5d_' runTimestamp]);

if ~exist(projectRootDir,'dir'), mkdir(projectRootDir); end
if ~exist(checkpointDirectory,'dir'), mkdir(checkpointDirectory); end
if cfg.SAVE_FIGS && ~exist(figureOutputDirectory,'dir'), mkdir(figureOutputDirectory); end

paths = struct( ...
    'projectRootDir', projectRootDir, ...
    'checkpointDirectory', checkpointDirectory, ...
    'checkpointFilePath', checkpointFilePath, ...
    'progressWorkbookPath', progressWorkbookPath, ...
    'figureOutputDirectory', figureOutputDirectory);

% Model path
comsolModelFile = cfg.comsolModelFile;
if ~isfile(comsolModelFile)
    error('COMSOL model not found at %s. Update config/comsolModelFile.', comsolModelFile);
end
if isfield(cfg,'comsolProjectRootDir') && exist(cfg.comsolProjectRootDir,'dir')
    addpath(genpath(cfg.comsolProjectRootDir));
end

% Tags/params
tags = struct( ...
    'STUDY_TAG', 'std1', ...
    'PARAM_HAU', 'h_au', ...
    'PARAM_HCEY', 'h_ceyig', ...
    'PARAM_LDOM', 'L_domain', ...
    'PARAM_LDEN', 'l_dente', ...
    'PARAM_HSI', 'h_si', ...
    'PARAM_N', 'n', ...
    'ALPHA_NAME', 'alpha', ...
    'MSIGN_NAME', 'm', ...
    'M_PLUS', '1', ...
    'M_MINUS', '-1', ...
    'RPLUS_TABLE_TAG', 'tblRplus', ...
    'RMINUS_TABLE_TAG', 'tblRminus');

% Checkpoint resume
pointsSinceCheckpoint = 0;
resumeFromCheckpoint = false;
checkpointData = struct;
resumeStageTag = "";
if isfile(checkpointFilePath)
    try
        checkpointFileContents = load(checkpointFilePath,'checkpointData');
        checkpointData = checkpointFileContents.checkpointData;
        resumeFromCheckpoint = true;
        resumeStageTag = string(checkpointData.stage);
        fprintf('>>> resumeFromCheckpoint detected | stage=%s | done_points=%d | runsCompletedGlobal=%d\n', ...
            resumeStageTag, checkpointData.done_points, checkpointData.runsCompletedGlobal);
    catch
        warning('Checkpoint exists but could not be loaded. Starting fresh.');
    end
end
resumeInfo = struct('resumeFromCheckpoint', resumeFromCheckpoint, 'resumeStageTag', resumeStageTag, 'checkpointData', checkpointData);

diary(fullfile(checkpointDirectory, sprintf('runlog_%s.txt', datestr(now,'yyyymmdd_HHMMss'))));

% Load model
ModelUtil.clear;
model = mphload(comsolModelFile);
ModelUtil.showProgress(true);

% Global planning
plan = computePlan(cfg);
globalTimerStart = tic;
runState = struct( ...
    'runsCompletedGlobal', 0, ...
    'pointsSinceCheckpoint', pointsSinceCheckpoint, ...
    'globalRunTargetEstimate', plan.globalRunTargetEstimate, ...
    'isTotalRunEstimateExact', plan.isTotalRunEstimateExact, ...
    'globalTimerStart', globalTimerStart, ...
    'checkpointEveryPoints', cfg.checkpointEveryPoints, ...
    'MAX_RUNS', cfg.MAX_RUNS);

fprintf('START\n');
fprintf('  COARSE (exact):    %d runs\n', plan.coarseTotalRuns);
fprintf('  FINE   (estimate): %d runs (TOP-%d)\n', plan.fineTotalRunsEstimate, cfg.topKCoarse);
fprintf('  SUPER  (exact):    %d runs (TOP-%d)\n', plan.superTotalRuns, cfg.topKFine);
fprintf('  EXTRAS (exact):    %d runs\n', plan.extraRunsFixed);
fprintf('  TOTAL  (estimate): %d runs\n\n', plan.globalRunTargetEstimate);

% Stage pipeline
[coarseResultsTable, coarseSeedCandidates, runState] = runCoarseStage(model, cfg, tags, paths, plan, runState, resumeInfo);
[fineResultsTable, superSeedCandidates, runState, plan] = runFineStage(model, cfg, tags, paths, plan, runState, resumeInfo, coarseResultsTable, coarseSeedCandidates);
[superResultsTable, bestTradeoffCandidate, bestTmokeCandidate, bestSensitivityCandidate, runState] = runSuperStage(model, cfg, tags, paths, plan, runState, resumeInfo, coarseResultsTable, fineResultsTable, superSeedCandidates);
[validation, runState, bestTradeoffCandidate, bestTmokeCandidate, bestSensitivityCandidate] = runValidationStage(model, cfg, tags, paths, runState, resumeInfo, coarseResultsTable, fineResultsTable, superResultsTable, bestTradeoffCandidate, bestTmokeCandidate, bestSensitivityCandidate);
[bestFullTable, runState] = runFullStage(model, cfg, tags, paths, runState, bestTradeoffCandidate, validation);

persistResults(paths, coarseResultsTable, fineResultsTable, superResultsTable, validation.bestDenseTable, bestFullTable);
makeFinalPlots(cfg, paths, validation, coarseResultsTable, fineResultsTable, superResultsTable, bestTradeoffCandidate);

elapsed_total = toc(runState.globalTimerStart);
printSummary(bestTradeoffCandidate, bestTmokeCandidate, bestSensitivityCandidate, validation.sensitivityDense, runState.runsCompletedGlobal, elapsed_total);

diary off;
end
