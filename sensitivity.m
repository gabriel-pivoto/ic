clear; clc; close all; format long; tic
import com.comsol.model.*;
import com.comsol.model.util.*;
projectRootDir  = 'D:\Gabriel Pivoto\projetoIC';
comsolModelFile = fullfile(projectRootDir,'modelosimplificado.mph');
addpath(genpath(projectRootDir));
MAX_RUNS = 20000;
checkpointDirectory   = fullfile(projectRootDir,'checkpoints');
if ~exist(checkpointDirectory,'dir'), mkdir(checkpointDirectory); end
checkpointFilePath    = fullfile(checkpointDirectory,'sensitivity_4d_checkpoint.mat');
progressWorkbookPath  = fullfile(checkpointDirectory,'sensitivity_4d_progress.xlsx');
checkpointSchemaTag   = 'tracked_same_peak_sensitivity_v1';
checkpointEveryPoints = 10;
pointsSinceCheckpoint = 0;
resumeFromCheckpoint = false;
checkpointData = struct;
resumeStageTag = "";
if isfile(checkpointFilePath)
    try
        checkpointFileContents = load(checkpointFilePath,'checkpointData');
        checkpointData = checkpointFileContents.checkpointData;
        if isfield(checkpointData,'schema') && strcmp(checkpointData.schema, checkpointSchemaTag)
            resumeFromCheckpoint = true;
            resumeStageTag = string(checkpointData.stage);
            fprintf('>>> resumeFromCheckpoint detected | stage=%s | done_points=%d | runsCompletedGlobal=%d\\n', ...
                resumeStageTag, checkpointData.done_points, checkpointData.runsCompletedGlobal);
        else
            warning('Checkpoint exists but uses an older sensitivity metric. Starting fresh.');
            checkpointData = struct;
        end
    catch
        warning('Checkpoint exists but could not be loaded. Starting fresh.');
    end
end
diary(fullfile(checkpointDirectory, sprintf('runlog_%s.txt', datestr(now,'yyyymmdd_HHMMss'))));
STUDY_TAG   = 'std1';
PARAM_HAU   = 'h_au';
PARAM_LDOM  = 'L_domain';
PARAM_LDEN  = 'l_dente';
PARAM_HSI   = 'h_si';
PARAM_N     = 'n';
ALPHA_NAME  = 'alpha';
MSIGN_NAME  = 'm';
M_PLUS      = '1';
M_MINUS     = '-1';
TPLUS_TABLE_TAG  = 'tblTplus';
TMINUS_TABLE_TAG = 'tblTminus';
fastRefractiveIndexSamples = [1.30, 1.33, 1.36];
trackingReferenceRefractiveIndex = 1.33;
trackingHalfWindowDeg = 20;
baselineRefractiveIndex = trackingReferenceRefractiveIndex;
validationRefractiveIndexList = [1.30, 1.33, 1.36];
domainPeriodGridNm = 800:50:850;
toothWidthGridNm  = 500:50:600;
siliconHeightGridNm     = [220, 240, 260];
goldHeightGridNm     = 20:10:60;
alphaCoarseRange = [0, 1.0, 89];
alphaFineStep  = 0.1;
alphaSuperStep = 0.01;
alphaDenseStep = 0.01;
globalAlphaSweepMinDeg = alphaCoarseRange(1);
globalAlphaSweepMaxDeg = alphaCoarseRange(3);
topKCoarse = 1;
topKFine   = 1;
fineGoldHeightDelta   = 2;    fineGoldHeightStep   = 1;
fineDomainPeriodDelta = 10;   fineDomainPeriodStep = 5;
fineToothWidthDelta   = 10;   fineToothWidthStep   = 5;
fineSiliconHeightDelta= 5;    fineSiliconHeightStep= 5;
superGoldHeightDelta  = 1;    superGoldHeightStep  = 0.5;
superDomainPeriodDelta = 4;   superDomainPeriodStep = 2;
superToothWidthDelta   = 4;   superToothWidthStep   = 2;
superSiliconHeightDelta= 2;   superSiliconHeightStep= 1;
fineAlphaHalfSpan   = 5;
superAlphaHalfSpan  = 2;
fineAlphaHalfSpanSensitivity  = fineAlphaHalfSpan  + 1;
superAlphaHalfSpanSensitivity = superAlphaHalfSpan + 2;
SAVE_SNAPSHOT = true;
MAKE_PLOTS    = true;
PLOT_LIVE     = true;
SAVE_PHASE_PLOTS = true;
SAVE_FIGS   = true;
FIG_FORMATS = {'png','pdf'};
runTimestamp   = datestr(now,'yyyy-mm-dd_HHMMSS');
figureOutputDirectory     = fullfile(projectRootDir,'plots', ['sensitivity_4d_' runTimestamp]);
phaseFigureOutputDirectory = fullfile(figureOutputDirectory, 'phase_outputs');
if SAVE_FIGS && ~exist(figureOutputDirectory,'dir'), mkdir(figureOutputDirectory); end
if SAVE_PHASE_PLOTS && ~exist(phaseFigureOutputDirectory,'dir'), mkdir(phaseFigureOutputDirectory); end
SAVE_ITER_PLOTS = true;
if SAVE_ITER_PLOTS
    iterationFigureOutputDirectory = fullfile(figureOutputDirectory, 'iteration_tmoke');
    if ~exist(iterationFigureOutputDirectory,'dir'), mkdir(iterationFigureOutputDirectory); end
else
    iterationFigureOutputDirectory = '';
end
ModelUtil.clear;
model = mphload(comsolModelFile);
ModelUtil.showProgress(true);
globalTimerStart = tic;
runsCompletedGlobal = 0;
isTotalRunEstimateExact   = false;
numDomainPeriodPoints = numel(domainPeriodGridNm);
numToothWidthPoints   = numel(toothWidthGridNm);
numSiliconHeightPoints= numel(siliconHeightGridNm);
numGoldHeightPoints   = numel(goldHeightGridNm);
runsPerSearchPoint = 2 * numel(fastRefractiveIndexSamples);
coarseTotalPoints  = numDomainPeriodPoints*numToothWidthPoints*numSiliconHeightPoints*numGoldHeightPoints;
coarseTotalRuns = runsPerSearchPoint * coarseTotalPoints;
fineGoldHeightCount = 1 + floor(2*fineGoldHeightDelta   / fineGoldHeightStep);
fineDomainPeriodCount = 1 + floor(2*fineDomainPeriodDelta  / fineDomainPeriodStep);
fineToothWidthCount = 1 + floor(2*fineToothWidthDelta  / fineToothWidthStep);
fineSiliconHeightCount = 1 + floor(2*fineSiliconHeightDelta   / fineSiliconHeightStep);
finePointsPerSeedEstimate = fineDomainPeriodCount*fineToothWidthCount*fineSiliconHeightCount*fineGoldHeightCount;
fineTotalPointsEstimate    = topKCoarse * finePointsPerSeedEstimate;
fineTotalRunsEstimate   = runsPerSearchPoint * fineTotalPointsEstimate;
superGoldHeightCount = 1 + floor(2*superGoldHeightDelta   / superGoldHeightStep);
superDomainPeriodCount = 1 + floor(2*superDomainPeriodDelta  / superDomainPeriodStep);
superToothWidthCount = 1 + floor(2*superToothWidthDelta  / superToothWidthStep);
superSiliconHeightCount = 1 + floor(2*superSiliconHeightDelta   / superSiliconHeightStep);
superPointsPerSeed = superDomainPeriodCount*superToothWidthCount*superSiliconHeightCount*superGoldHeightCount;
superTotalPoints    = topKFine * superPointsPerSeed;
superTotalRuns   = runsPerSearchPoint * superTotalPoints;
extraRunsFixed = 2*numel(validationRefractiveIndexList) + (SAVE_SNAPSHOT*2);
globalRunTargetEstimate = coarseTotalRuns + fineTotalRunsEstimate + superTotalRuns + extraRunsFixed;
isTotalRunEstimateExact = false;
fprintf('START\n');
fprintf('  Stage 1 (exact):    %d runs\n', coarseTotalRuns);
fprintf('  Stage 2 (estimate): %d runs (TOP-%d)\n', fineTotalRunsEstimate, topKCoarse);
fprintf('  Stage 3 (exact):    %d runs (TOP-%d)\n', superTotalRuns, topKFine);
fprintf('  EXTRAS (exact):    %d runs\n', extraRunsFixed);
fprintf('  TOTAL  (estimate): %d runs\n\n', globalRunTargetEstimate);
coarseResultsTable = [];
coarseSeedCandidates   = [];
skipCoarseStage = resumeFromCheckpoint && any(strcmp(resumeStageTag, ["Stage 2","Stage 3","Stage 4","FINAL"]));
if skipCoarseStage
    if isfield(checkpointData.payload,'coarseResultsTable'), coarseResultsTable = checkpointData.payload.coarseResultsTable; end
    if isfield(checkpointData.payload,'coarseSeedCandidates'),   coarseSeedCandidates   = checkpointData.payload.coarseSeedCandidates;   end
    if (isempty(coarseSeedCandidates) || height(coarseSeedCandidates)==0) && ~isempty(coarseResultsTable)
        coarseSeedCandidates = selectTopK_single_abs(coarseResultsTable, topKCoarse, 'S_est_deg_per_RIU');
        warning('SKIP Stage 1: coarseSeedCandidates missing in checkpoint -> rebuilt from coarseResultsTable.');
    end
    fprintf('SKIP Stage 1 -> restored coarseResultsTable (rows=%d), coarseSeedCandidates (rows=%d)\n', ...
        size(coarseResultsTable,1), height(coarseSeedCandidates));
    if MAKE_PLOTS && SAVE_PHASE_PLOTS
        saveStageCandidatePlots(coarseResultsTable, 'Stage 1', phaseFigureOutputDirectory, FIG_FORMATS);
    end
else
    fprintf('STAGE 1 - EXACT: %d runs\n', coarseTotalRuns);
    stageRunsStart = runsCompletedGlobal;
    stageTotalRuns = coarseTotalRuns;
    stageTimerStart = tic;
    coarseRows = [];
    coarsePointIndex = 0;
    if resumeFromCheckpoint && resumeStageTag=="Stage 1"
        runsCompletedGlobal = checkpointData.runsCompletedGlobal;
        if isfield(checkpointData.payload,'coarseRows'), coarseRows = checkpointData.payload.coarseRows; end
        coarsePointIndex = checkpointData.done_points;
        fprintf('>>> Stage 1 resume: skipping %d points already computed.\n', coarsePointIndex);
    end
    for domainPeriodIdx = 1:numDomainPeriodPoints
        setParamNm(model, PARAM_LDOM, domainPeriodGridNm(domainPeriodIdx));
        for toothWidthIdx = 1:numToothWidthPoints
            setParamNm(model, PARAM_LDEN, toothWidthGridNm(toothWidthIdx));
            for siliconHeightIdx = 1:numSiliconHeightPoints
                setParamNm(model, PARAM_HSI, siliconHeightGridNm(siliconHeightIdx));
                for goldHeightIdx = 1:numGoldHeightPoints
                    setParamNm(model, PARAM_HAU, goldHeightGridNm(goldHeightIdx));
                        coarsePointIndex = coarsePointIndex + 1;
                        if resumeFromCheckpoint && resumeStageTag=="Stage 1" && coarsePointIndex <= checkpointData.done_points
                            continue;
                        end
                        fastTrackedEval = evaluateTrackedSensitivityAndCurves( ...
                            model, STUDY_TAG, PARAM_N, fastRefractiveIndexSamples, ...
                            trackingReferenceRefractiveIndex, trackingHalfWindowDeg, ...
                            ALPHA_NAME, MSIGN_NAME, ...
                            alphaCoarseRange(1), alphaCoarseRange(2), alphaCoarseRange(3), ...
                            globalAlphaSweepMinDeg, globalAlphaSweepMaxDeg, false, ...
                            M_PLUS, M_MINUS, TPLUS_TABLE_TAG, TMINUS_TABLE_TAG);
                        referenceFastIndex = fastTrackedEval.referenceIndex;
                        highestFastIndex = fastTrackedEval.highestNIndex;
                        alphaGridFastBase = fastTrackedEval.alphaGridsByN{referenceFastIndex};
                        transmissionPlusFastBase = fastTrackedEval.transmissionPlusByN{referenceFastIndex};
                        transmissionMinusFastBase = fastTrackedEval.transmissionMinusByN{referenceFastIndex};
                        tmokeCurveFastBase = fastTrackedEval.tmokeCurvesByN{referenceFastIndex};
                        tmokePeakMagnitudeFastBase = fastTrackedEval.trackedTmokeAbsByN(referenceFastIndex);
                        alphaAtPeakFastBase = fastTrackedEval.alphaPeakDegreesByN(referenceFastIndex);
                        tmokeAtPeakFastBase = fastTrackedEval.tmokeAtTrackedPeakByN(referenceFastIndex);
                        alphaAtPeakFastHigh = fastTrackedEval.alphaPeakDegreesByN(highestFastIndex);
                        sensitivityEstimateFast = fastTrackedEval.sensitivitySlope;
                        tmokePeakIndexFastBase = find(abs(alphaGridFastBase - alphaAtPeakFastBase) < 1e-9, 1, 'first');
                        if ~isempty(tmokePeakIndexFastBase) && ...
                                (tmokePeakIndexFastBase == 1 || tmokePeakIndexFastBase == numel(tmokeCurveFastBase))
                            logline('WARN Stage 1 (n=%.2f) reference peak at alpha-edge (alpha*=%.3f in [%.3f, %.3f])\n', ...
                                baselineRefractiveIndex, alphaAtPeakFastBase, alphaGridFastBase(1), alphaGridFastBase(end));
                        end
                        if PLOT_LIVE
                            updateLivePlot('Stage 1', ...
                                domainPeriodGridNm(domainPeriodIdx), toothWidthGridNm(toothWidthIdx), siliconHeightGridNm(siliconHeightIdx), ...
                                goldHeightGridNm(goldHeightIdx), ...
                                alphaGridFastBase, tmokeCurveFastBase, alphaAtPeakFastBase, tmokeAtPeakFastBase, transmissionPlusFastBase, transmissionMinusFastBase, ...
                                iterationFigureOutputDirectory);
                        end
                        coarseRows = [coarseRows; ...
                            domainPeriodGridNm(domainPeriodIdx), toothWidthGridNm(toothWidthIdx), siliconHeightGridNm(siliconHeightIdx), ...
                            goldHeightGridNm(goldHeightIdx), ...
                            tmokePeakMagnitudeFastBase, alphaAtPeakFastBase, tmokeAtPeakFastBase, ...
                            alphaAtPeakFastHigh, sensitivityEstimateFast];
                        runsCompletedGlobal = runsCompletedGlobal + runsPerSearchPoint;
                        stageRunsCompleted = runsCompletedGlobal - stageRunsStart;
                        [stageCompletionFraction, stageEtaSeconds]   = frac_eta(stageRunsCompleted, stageTotalRuns, toc(stageTimerStart));
                        [globalCompletionFraction, globalEtaSeconds] = frac_eta(runsCompletedGlobal, globalRunTargetEstimate, toc(globalTimerStart));
                        logline(['Stage 1 | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_au=%5.1f' ...
                                 ' | |TM|=%.5f @ alpha=%.2f deg | sensitivityEstimateFast=%.4f deg/RIU' ...
                                 ' [Stage %5.1f%% | ETA %s | Global %s%% | ETA %s]\n'], ...
                            domainPeriodGridNm(domainPeriodIdx), toothWidthGridNm(toothWidthIdx), siliconHeightGridNm(siliconHeightIdx), ...
                            goldHeightGridNm(goldHeightIdx), ...
                            tmokePeakMagnitudeFastBase, alphaAtPeakFastBase, sensitivityEstimateFast, ...
                            100*stageCompletionFraction, fmt_time_long(stageEtaSeconds), fmt_pct(globalCompletionFraction), fmt_eta_flag(globalEtaSeconds, isTotalRunEstimateExact));
                        pointsSinceCheckpoint = pointsSinceCheckpoint + 1;
                        payload = struct('coarseRows',coarseRows);
                        pointsSinceCheckpoint = maybe_checkpoint( ...
                            'Stage 1', checkpointEveryPoints, checkpointFilePath, progressWorkbookPath, ...
                            runsCompletedGlobal, pointsSinceCheckpoint, coarsePointIndex, payload, ...
                            [], [], [], []);
                end
            end
        end
    end
    coarseResultsTable = array2table(coarseRows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_au_nm', ...
         'maxAbsTMOKE_base','alpha_peak_base_deg','TMOKE_at_peak_base', ...
         'alpha_peak_high_deg','S_est_deg_per_RIU'});
    coarseSeedCandidates = selectTopK_single_abs(coarseResultsTable, topKCoarse, 'S_est_deg_per_RIU');
    save_checkpoint(checkpointFilePath, 'Stage 2', runsCompletedGlobal, 0, struct('coarseResultsTable', coarseResultsTable, 'coarseSeedCandidates', coarseSeedCandidates));
    write_progress_xlsx(progressWorkbookPath, 'coarse', coarseResultsTable);
    if MAKE_PLOTS && SAVE_PHASE_PLOTS
        saveStageCandidatePlots(coarseResultsTable, 'Stage 1', phaseFigureOutputDirectory, FIG_FORMATS);
    end
end
if isempty(coarseSeedCandidates) && resumeFromCheckpoint && resumeStageTag=="Stage 2"
    if isfield(checkpointData.payload,'coarseSeedCandidates')
        coarseSeedCandidates = checkpointData.payload.coarseSeedCandidates;
    else
        error('Checkpoint at Stage 2 has no "coarseSeedCandidates". Rerun Stage 1.');
    end
end
fineTotalPoints = 0;
for s = 1:height(coarseSeedCandidates)
    goldHeightList = max(min(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   - fineGoldHeightDelta)  : fineGoldHeightStep  : min(max(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   + fineGoldHeightDelta);
    domainPeriodList = max(min(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)- fineDomainPeriodDelta) : fineDomainPeriodStep : min(max(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)+ fineDomainPeriodDelta);
    toothWidthList = max(min(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) - fineToothWidthDelta) : fineToothWidthStep : min(max(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) + fineToothWidthDelta);
    siliconHeightList = max(min(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    - fineSiliconHeightDelta)  : fineSiliconHeightStep  : min(max(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    + fineSiliconHeightDelta);
    fineTotalPoints = fineTotalPoints + numel(domainPeriodList)*numel(toothWidthList)*numel(siliconHeightList)*numel(goldHeightList);
end
fineTotalRuns = runsPerSearchPoint * fineTotalPoints;
globalRunTargetEstimate = coarseTotalRuns + fineTotalRuns + superTotalRuns + extraRunsFixed;
isTotalRunEstimateExact = true;
fprintf('\nSTAGE 2 - EXACT: %d runs (TOP-%d)\n', fineTotalRuns, height(coarseSeedCandidates));
fprintf('GLOBAL TOTAL (EXACT): %d runs\n\n', globalRunTargetEstimate);
if globalRunTargetEstimate > MAX_RUNS
    error('Planned total exceeds %d runs (%d). Adjust deltas/steps/TOPK.', MAX_RUNS, globalRunTargetEstimate);
end
fineResultsTable = [];
superSeedCandidates = [];
if resumeFromCheckpoint && any(strcmp(resumeStageTag, ["Stage 3","Stage 4","FINAL"]))
    if isfield(checkpointData.payload,'fineResultsTable'),      fineResultsTable = checkpointData.payload.fineResultsTable; end
    if isfield(checkpointData.payload,'superSeedCandidates'), superSeedCandidates = checkpointData.payload.superSeedCandidates; end
    if (isempty(superSeedCandidates) || height(superSeedCandidates)==0) && ~isempty(fineResultsTable)
        superSeedCandidates = selectTopK_single_abs(fineResultsTable, topKFine, 'S_est_deg_per_RIU');
        warning('SKIP Stage 2: superSeedCandidates missing in checkpoint -> rebuilt from fineResultsTable.');
    end
    fprintf('SKIP Stage 2 -> restored fineResultsTable (rows=%d), superSeedCandidates (rows=%d)\n', ...
        size(fineResultsTable,1), height(superSeedCandidates));
    if MAKE_PLOTS && SAVE_PHASE_PLOTS
        saveStageCandidatePlots(fineResultsTable, 'Stage 2', phaseFigureOutputDirectory, FIG_FORMATS);
    end
else
    fineRows = [];
    stageRunsStart = runsCompletedGlobal;
    stageTotalRuns = fineTotalRuns;
    stageTimerStart = tic;
    finePointIndex = 0;
    if resumeFromCheckpoint && resumeStageTag=="Stage 2" && checkpointData.done_points > 0
        runsCompletedGlobal = checkpointData.runsCompletedGlobal;
        finePointIndex   = checkpointData.done_points;
        if isfield(checkpointData.payload,'fineRows'), fineRows = checkpointData.payload.fineRows; end
        if isfield(checkpointData.payload,'coarseSeedCandidates'),    coarseSeedCandidates    = checkpointData.payload.coarseSeedCandidates;    end
        fprintf('>>> Stage 2 resume: skipping %d points already computed.\n', finePointIndex);
    end
    for s = 1:height(coarseSeedCandidates)
        alphaWindowCenterDeg = coarseSeedCandidates.alpha_peak_base_deg(s);
        goldHeightList = max(min(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   - fineGoldHeightDelta)  : fineGoldHeightStep  : min(max(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   + fineGoldHeightDelta);
        domainPeriodList = max(min(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)- fineDomainPeriodDelta) : fineDomainPeriodStep : min(max(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)+ fineDomainPeriodDelta);
        toothWidthList = max(min(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) - fineToothWidthDelta) : fineToothWidthStep : min(max(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) + fineToothWidthDelta);
        siliconHeightList = max(min(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    - fineSiliconHeightDelta)  : fineSiliconHeightStep  : min(max(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    + fineSiliconHeightDelta);
        alphaStartDeg = max(0,  alphaWindowCenterDeg - fineAlphaHalfSpanSensitivity);
        alphaStopDeg  = min(89, alphaWindowCenterDeg + fineAlphaHalfSpanSensitivity);
        for domainPeriodIdx = 1:numel(domainPeriodList)
            setParamNm(model, PARAM_LDOM, domainPeriodList(domainPeriodIdx));
            for toothWidthIdx = 1:numel(toothWidthList)
                setParamNm(model, PARAM_LDEN, toothWidthList(toothWidthIdx));
                for siliconHeightIdx = 1:numel(siliconHeightList)
                    setParamNm(model, PARAM_HSI, siliconHeightList(siliconHeightIdx));
                    for goldHeightIdx = 1:numel(goldHeightList)
                        setParamNm(model, PARAM_HAU, goldHeightList(goldHeightIdx));
                            finePointIndex = finePointIndex + 1;
                            if resumeFromCheckpoint && resumeStageTag=="Stage 2" && finePointIndex <= checkpointData.done_points
                                continue;
                            end
                            fastTrackedEval = evaluateTrackedSensitivityAndCurves( ...
                                model, STUDY_TAG, PARAM_N, fastRefractiveIndexSamples, ...
                                trackingReferenceRefractiveIndex, trackingHalfWindowDeg, ...
                                ALPHA_NAME, MSIGN_NAME, ...
                                alphaStartDeg, alphaFineStep, alphaStopDeg, ...
                                globalAlphaSweepMinDeg, globalAlphaSweepMaxDeg, false, ...
                                M_PLUS, M_MINUS, TPLUS_TABLE_TAG, TMINUS_TABLE_TAG);
                            referenceFastIndex = fastTrackedEval.referenceIndex;
                            highestFastIndex = fastTrackedEval.highestNIndex;
                            alphaGridFastBase = fastTrackedEval.alphaGridsByN{referenceFastIndex};
                            transmissionPlusFastBase = fastTrackedEval.transmissionPlusByN{referenceFastIndex};
                            transmissionMinusFastBase = fastTrackedEval.transmissionMinusByN{referenceFastIndex};
                            tmokeCurveFastBase = fastTrackedEval.tmokeCurvesByN{referenceFastIndex};
                            tmokePeakMagnitudeFastBase = fastTrackedEval.trackedTmokeAbsByN(referenceFastIndex);
                            alphaAtPeakFastBase = fastTrackedEval.alphaPeakDegreesByN(referenceFastIndex);
                            tmokeAtPeakFastBase = fastTrackedEval.tmokeAtTrackedPeakByN(referenceFastIndex);
                            alphaAtPeakFastHigh = fastTrackedEval.alphaPeakDegreesByN(highestFastIndex);
                            sensitivityEstimateFast = fastTrackedEval.sensitivitySlope;
                            if PLOT_LIVE
                                updateLivePlot('Stage 2', domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), goldHeightList(goldHeightIdx), ...
                                    alphaGridFastBase, tmokeCurveFastBase, alphaAtPeakFastBase, tmokeAtPeakFastBase, transmissionPlusFastBase, transmissionMinusFastBase, ...
                                    iterationFigureOutputDirectory);
                            end
                            fineRows = [fineRows; ...
                                domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), goldHeightList(goldHeightIdx), ...
                                tmokePeakMagnitudeFastBase, alphaAtPeakFastBase, tmokeAtPeakFastBase, alphaAtPeakFastHigh, sensitivityEstimateFast];
                            runsCompletedGlobal = runsCompletedGlobal + runsPerSearchPoint;
                            stageRunsCompleted = runsCompletedGlobal - stageRunsStart;
                            [stageCompletionFraction, stageEtaSeconds]   = frac_eta(stageRunsCompleted, stageTotalRuns, toc(stageTimerStart));
                            [globalCompletionFraction, globalEtaSeconds] = frac_eta(runsCompletedGlobal, globalRunTargetEstimate, toc(globalTimerStart));
                            logline(['Stage 2 | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_au=%5.1f' ...
                                     ' | |TM|=%.5f @ alpha=%.3f deg | sensitivityEstimateFast=%.4f' ...
                                     ' [Stage %5.1f%% | ETA %s | Global %5.1f%% | ETA %s]\n'], ...
                                domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), goldHeightList(goldHeightIdx), ...
                                tmokePeakMagnitudeFastBase, alphaAtPeakFastBase, sensitivityEstimateFast, ...
                                100*stageCompletionFraction, fmt_time_long(stageEtaSeconds), 100*globalCompletionFraction, fmt_time_long(globalEtaSeconds));
                            pointsSinceCheckpoint = pointsSinceCheckpoint + 1;
                            payload = struct('fineRows',fineRows, ...
                                'coarseResultsTable', coarseResultsTable, 'coarseSeedCandidates', coarseSeedCandidates);
                            pointsSinceCheckpoint = maybe_checkpoint( ...
                                'Stage 2', checkpointEveryPoints, checkpointFilePath, progressWorkbookPath, ...
                                runsCompletedGlobal, pointsSinceCheckpoint, finePointIndex, payload, ...
                                coarseResultsTable, [], [], []);
                    end
                end
            end
        end
    end
    fineResultsTable = array2table(fineRows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_au_nm', ...
         'maxAbsTMOKE_base','alpha_peak_base_deg','TMOKE_at_peak_base', ...
         'alpha_peak_high_deg','S_est_deg_per_RIU'});
    superSeedCandidates = selectTopK_single_abs(fineResultsTable, topKFine, 'S_est_deg_per_RIU');
    save_checkpoint(checkpointFilePath, 'Stage 3', runsCompletedGlobal, 0, struct( ...
        'coarseResultsTable', coarseResultsTable, 'coarseSeedCandidates', coarseSeedCandidates, ...
        'fineResultsTable', fineResultsTable, 'superSeedCandidates', superSeedCandidates));
    write_progress_xlsx(progressWorkbookPath, 'fine', fineResultsTable);
    if MAKE_PLOTS && SAVE_PHASE_PLOTS
        saveStageCandidatePlots(fineResultsTable, 'Stage 2', phaseFigureOutputDirectory, FIG_FORMATS);
    end
end
if isempty(superSeedCandidates) && resumeFromCheckpoint && resumeStageTag=="Stage 3"
    if isfield(checkpointData.payload,'superSeedCandidates')
        superSeedCandidates = checkpointData.payload.superSeedCandidates;
    else
        error('Checkpoint at Stage 3 has no "superSeedCandidates". Rerun Stage 2.');
    end
end
superTotalPoints = 0;
for s = 1:height(superSeedCandidates)
    goldHeightList = (superSeedCandidates.h_au_nm(s)    - superGoldHeightDelta)  : superGoldHeightStep  : (superSeedCandidates.h_au_nm(s)    + superGoldHeightDelta);
    domainPeriodList = (superSeedCandidates.L_domain_nm(s)- superDomainPeriodDelta) : superDomainPeriodStep : (superSeedCandidates.L_domain_nm(s)+ superDomainPeriodDelta);
    toothWidthList = (superSeedCandidates.l_dente_nm(s) - superToothWidthDelta) : superToothWidthStep : (superSeedCandidates.l_dente_nm(s) + superToothWidthDelta);
    siliconHeightList = (superSeedCandidates.h_si_nm(s)    - superSiliconHeightDelta)  : superSiliconHeightStep  : (superSeedCandidates.h_si_nm(s)    + superSiliconHeightDelta);
    superTotalPoints = superTotalPoints + numel(domainPeriodList)*numel(toothWidthList)*numel(siliconHeightList)*numel(goldHeightList);
end
superTotalRuns = runsPerSearchPoint * superTotalPoints;
fprintf('\nSTAGE 3 - EXACT: %d runs (TOP-%d)\n', superTotalRuns, height(superSeedCandidates));
fprintf('GLOBAL TOTAL (EXACT): %d runs\n\n', globalRunTargetEstimate);
superResultsTable = [];
bestSensitivityCandidate = table();
if resumeFromCheckpoint && any(strcmp(resumeStageTag, ["Stage 4","FINAL"]))
    if isfield(checkpointData.payload,'superResultsTable'), superResultsTable = checkpointData.payload.superResultsTable; end
    if isfield(checkpointData.payload,'bestSensitivityCandidate'), bestSensitivityCandidate = checkpointData.payload.bestSensitivityCandidate; end
    fprintf('SKIP Stage 3 -> restored superResultsTable (rows=%d) and best selection.\n', size(superResultsTable,1));
    if MAKE_PLOTS && SAVE_PHASE_PLOTS
        saveStageCandidatePlots(superResultsTable, 'Stage 3', phaseFigureOutputDirectory, FIG_FORMATS);
    end
else
    superRows = [];
    stageRunsStart = runsCompletedGlobal;
    stageTotalRuns = superTotalRuns;
    stageTimerStart = tic;
    superPointIndex = 0;
    if resumeFromCheckpoint && resumeStageTag=="Stage 3" && checkpointData.done_points > 0
        runsCompletedGlobal = checkpointData.runsCompletedGlobal;
        superPointIndex  = checkpointData.done_points;
        if isfield(checkpointData.payload,'superRows'),  superRows  = checkpointData.payload.superRows;  end
        if isfield(checkpointData.payload,'superSeedCandidates'), superSeedCandidates = checkpointData.payload.superSeedCandidates; end
        fprintf('>>> Stage 3 resume: skipping %d points already computed.\n', superPointIndex);
    end
    for s = 1:height(superSeedCandidates)
        alphaWindowCenterDeg = superSeedCandidates.alpha_peak_base_deg(s);
        goldHeightList = (superSeedCandidates.h_au_nm(s)    - superGoldHeightDelta)  : superGoldHeightStep  : (superSeedCandidates.h_au_nm(s)    + superGoldHeightDelta);
        domainPeriodList = (superSeedCandidates.L_domain_nm(s)- superDomainPeriodDelta) : superDomainPeriodStep : (superSeedCandidates.L_domain_nm(s)+ superDomainPeriodDelta);
        toothWidthList = (superSeedCandidates.l_dente_nm(s) - superToothWidthDelta) : superToothWidthStep : (superSeedCandidates.l_dente_nm(s) + superToothWidthDelta);
        siliconHeightList = (superSeedCandidates.h_si_nm(s)    - superSiliconHeightDelta)  : superSiliconHeightStep  : (superSeedCandidates.h_si_nm(s)    + superSiliconHeightDelta);
        alphaStartDeg = max(0,  alphaWindowCenterDeg - superAlphaHalfSpanSensitivity);
        alphaStopDeg  = min(89, alphaWindowCenterDeg + superAlphaHalfSpanSensitivity);
        for domainPeriodIdx = 1:numel(domainPeriodList)
            setParamNm(model, PARAM_LDOM, domainPeriodList(domainPeriodIdx));
            for toothWidthIdx = 1:numel(toothWidthList)
                setParamNm(model, PARAM_LDEN, toothWidthList(toothWidthIdx));
                for siliconHeightIdx = 1:numel(siliconHeightList)
                    setParamNm(model, PARAM_HSI, siliconHeightList(siliconHeightIdx));
                    for goldHeightIdx = 1:numel(goldHeightList)
                        setParamNm(model, PARAM_HAU, goldHeightList(goldHeightIdx));
                            superPointIndex = superPointIndex + 1;
                            if resumeFromCheckpoint && resumeStageTag=="Stage 3" && superPointIndex <= checkpointData.done_points
                                continue;
                            end
                            fastTrackedEval = evaluateTrackedSensitivityAndCurves( ...
                                model, STUDY_TAG, PARAM_N, fastRefractiveIndexSamples, ...
                                trackingReferenceRefractiveIndex, trackingHalfWindowDeg, ...
                                ALPHA_NAME, MSIGN_NAME, ...
                                alphaStartDeg, alphaSuperStep, alphaStopDeg, ...
                                globalAlphaSweepMinDeg, globalAlphaSweepMaxDeg, false, ...
                                M_PLUS, M_MINUS, TPLUS_TABLE_TAG, TMINUS_TABLE_TAG);
                            referenceFastIndex = fastTrackedEval.referenceIndex;
                            highestFastIndex = fastTrackedEval.highestNIndex;
                            alphaGridFastBase = fastTrackedEval.alphaGridsByN{referenceFastIndex};
                            transmissionPlusFastBase = fastTrackedEval.transmissionPlusByN{referenceFastIndex};
                            transmissionMinusFastBase = fastTrackedEval.transmissionMinusByN{referenceFastIndex};
                            tmokeCurveFastBase = fastTrackedEval.tmokeCurvesByN{referenceFastIndex};
                            tmokePeakMagnitudeFastBase = fastTrackedEval.trackedTmokeAbsByN(referenceFastIndex);
                            alphaAtPeakFastBase = fastTrackedEval.alphaPeakDegreesByN(referenceFastIndex);
                            tmokeAtPeakFastBase = fastTrackedEval.tmokeAtTrackedPeakByN(referenceFastIndex);
                            alphaAtPeakFastHigh = fastTrackedEval.alphaPeakDegreesByN(highestFastIndex);
                            sensitivityEstimateFast = fastTrackedEval.sensitivitySlope;
                            if PLOT_LIVE
                                updateLivePlot('Stage 3', domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), goldHeightList(goldHeightIdx), ...
                                    alphaGridFastBase, tmokeCurveFastBase, alphaAtPeakFastBase, tmokeAtPeakFastBase, transmissionPlusFastBase, transmissionMinusFastBase, ...
                                    iterationFigureOutputDirectory);
                            end
                            superRows = [superRows; ...
                                domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), goldHeightList(goldHeightIdx), ...
                                tmokePeakMagnitudeFastBase, alphaAtPeakFastBase, tmokeAtPeakFastBase, alphaAtPeakFastHigh, sensitivityEstimateFast];
                            runsCompletedGlobal = runsCompletedGlobal + runsPerSearchPoint;
                            stageRunsCompleted = runsCompletedGlobal - stageRunsStart;
                            [stageCompletionFraction, stageEtaSeconds]   = frac_eta(stageRunsCompleted, stageTotalRuns, toc(stageTimerStart));
                            [globalCompletionFraction, globalEtaSeconds] = frac_eta(runsCompletedGlobal, globalRunTargetEstimate, toc(globalTimerStart));
                            logline(['Stage 3 | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_au=%5.1f' ...
                                     ' | |TM|=%.5f @ alpha=%.4f deg | sensitivityEstimateFast=%.4f' ...
                                     ' [Stage %5.1f%% | ETA %s | Global %5.1f%% | ETA %s]\n'], ...
                                domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), goldHeightList(goldHeightIdx), ...
                                tmokePeakMagnitudeFastBase, alphaAtPeakFastBase, sensitivityEstimateFast, ...
                                100*stageCompletionFraction, fmt_time_long(stageEtaSeconds), 100*globalCompletionFraction, fmt_time_long(globalEtaSeconds));
                            pointsSinceCheckpoint = pointsSinceCheckpoint + 1;
                            payload = struct('superRows',superRows, ...
                                'coarseResultsTable', coarseResultsTable, 'coarseSeedCandidates', coarseSeedCandidates, ...
                                'fineResultsTable', fineResultsTable, 'superSeedCandidates', superSeedCandidates);
                            pointsSinceCheckpoint = maybe_checkpoint( ...
                                'Stage 3', checkpointEveryPoints, checkpointFilePath, progressWorkbookPath, ...
                                runsCompletedGlobal, pointsSinceCheckpoint, superPointIndex, payload, ...
                                coarseResultsTable, fineResultsTable, [], []);
                    end
                end
            end
        end
    end
    superResultsTable = array2table(superRows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_au_nm', ...
         'maxAbsTMOKE_base','alpha_peak_base_deg','TMOKE_at_peak_base', ...
         'alpha_peak_high_deg','S_est_deg_per_RIU'});
    bestSensitivityCandidate = selectTopK_single_abs(superResultsTable, 1, 'S_est_deg_per_RIU');
    fprintf('\n===== BEST (Stage 3) =====\n');
    fprintf('Best |S| only:\n'); disp(bestSensitivityCandidate);
    save_checkpoint(checkpointFilePath, 'Stage 4', runsCompletedGlobal, 0, struct( ...
        'coarseResultsTable', coarseResultsTable, 'coarseSeedCandidates', coarseSeedCandidates, ...
        'fineResultsTable', fineResultsTable, 'superSeedCandidates', superSeedCandidates, 'superResultsTable', superResultsTable, ...
        'bestSensitivityCandidate', bestSensitivityCandidate));
    write_progress_xlsx(progressWorkbookPath, 'super', superResultsTable);
    if MAKE_PLOTS && SAVE_PHASE_PLOTS
        saveStageCandidatePlots(superResultsTable, 'Stage 3', phaseFigureOutputDirectory, FIG_FORMATS);
    end
end
if isempty(bestSensitivityCandidate) && resumeFromCheckpoint && any(strcmp(resumeStageTag, ["Stage 4","FINAL"]))
    bestSensitivityCandidate = checkpointData.payload.bestSensitivityCandidate;
end
Ldom_best = bestSensitivityCandidate.L_domain_nm(1);
Lden_best = bestSensitivityCandidate.l_dente_nm(1);
hsi_best  = bestSensitivityCandidate.h_si_nm(1);
hau_best  = bestSensitivityCandidate.h_au_nm(1);
setParamNm(model, PARAM_LDOM, Ldom_best);
setParamNm(model, PARAM_LDEN, Lden_best);
setParamNm(model, PARAM_HSI,  hsi_best);
setParamNm(model, PARAM_HAU,  hau_best);
validationTrackedEval = evaluateTrackedSensitivityAndCurves( ...
    model, STUDY_TAG, PARAM_N, validationRefractiveIndexList, ...
    trackingReferenceRefractiveIndex, trackingHalfWindowDeg, ...
    ALPHA_NAME, MSIGN_NAME, ...
    0, alphaDenseStep, 89, ...
    globalAlphaSweepMinDeg, globalAlphaSweepMaxDeg, true, ...
    M_PLUS, M_MINUS, TPLUS_TABLE_TAG, TMINUS_TABLE_TAG);
alphaPeakDegreesByN = validationTrackedEval.alphaPeakDegreesByN;
trackedTmokeAbsByN = validationTrackedEval.trackedTmokeAbsByN;
tmokeCurvesByN = validationTrackedEval.tmokeCurvesByN;
alphaGridsByN = validationTrackedEval.alphaGridsByN;
transmissionPlusByN = validationTrackedEval.transmissionPlusByN;
transmissionMinusByN = validationTrackedEval.transmissionMinusByN;
alphaVsNLinearFit = validationTrackedEval.linearFit;
sensitivityDense = validationTrackedEval.sensitivitySlope;
runsCompletedGlobal = runsCompletedGlobal + 2*numel(validationRefractiveIndexList);
fprintf('\n===== Stage 4 (bestSensitivityCandidate) =====\n');
fprintf('Best geom: Ldom=%g | Lden=%g | hsi=%g | hau=%g\n', ...
    Ldom_best, Lden_best, hsi_best, hau_best);
fprintf('alpha_peak(n) = ['); fprintf(' %.4f', alphaPeakDegreesByN); fprintf(' ]\n');
fprintf('sensitivityDense ~= %.6f deg/RIU (linear fit)\n', sensitivityDense);
[~, idxBase] = min(abs(validationRefractiveIndexList - baselineRefractiveIndex));
alphaBestDeg = alphaPeakDegreesByN(idxBase);
bestDenseTable = table();
for i = 1:numel(validationRefractiveIndexList)
    ncol = repmat(validationRefractiveIndexList(i), numel(alphaGridsByN{i}), 1);
    Ttmp = table( ...
        repmat(Ldom_best,numel(alphaGridsByN{i}),1), repmat(Lden_best,numel(alphaGridsByN{i}),1), repmat(hsi_best,numel(alphaGridsByN{i}),1), ...
        repmat(hau_best,numel(alphaGridsByN{i}),1), ...
        ncol, alphaGridsByN{i}(:), transmissionPlusByN{i}(:), transmissionMinusByN{i}(:), tmokeCurvesByN{i}(:), ...
        'VariableNames', {'L_domain_nm','l_dente_nm','h_si_nm','h_au_nm', ...
                          'n','alpha_deg','Tplus','Tminus','TMOKE'});
    bestDenseTable = [bestDenseTable; Ttmp];
end
save_checkpoint(checkpointFilePath, 'FINAL', runsCompletedGlobal, 0, struct( ...
    'coarseResultsTable', coarseResultsTable, 'coarseSeedCandidates', coarseSeedCandidates, ...
    'fineResultsTable', fineResultsTable, 'superSeedCandidates', superSeedCandidates, 'superResultsTable', superResultsTable, ...
    'bestSensitivityCandidate', bestSensitivityCandidate, ...
    'bestDenseTable', bestDenseTable, 'sensitivityDense', sensitivityDense, 'alphaBestDeg', alphaBestDeg));
write_progress_xlsx(progressWorkbookPath, 'dense', bestDenseTable);
if MAKE_PLOTS && SAVE_PHASE_PLOTS
    saveValidationPhasePlots(validationRefractiveIndexList, alphaGridsByN, tmokeCurvesByN, ...
        alphaPeakDegreesByN, trackedTmokeAbsByN, alphaVsNLinearFit, sensitivityDense, ...
        phaseFigureOutputDirectory, FIG_FORMATS);
end
if SAVE_SNAPSHOT
    setParamScalar(model, PARAM_N, baselineRefractiveIndex);
    setAlphaMSweep(model, STUDY_TAG, ALPHA_NAME, alphaBestDeg, 0.01, alphaBestDeg, MSIGN_NAME, sprintf('%s %s', M_PLUS, M_MINUS));
    model.study(STUDY_TAG).run;
    refreshDerivedValues(model);
    runsCompletedGlobal = runsCompletedGlobal + 2;
    snapshotTimestamp = datestr(now,'yyyymmdd_HHMMSS');
    snapshotFilePath = fullfile(projectRootDir, sprintf( ...
        'snapshot_4d_sensitivity_bestSensitivity_Ldom%4.0f_Lden%4.0f_hsi%3.0f_hau%.3f_n%.3f_alpha%.4f_%s.mph', ...
        Ldom_best, Lden_best, hsi_best, hau_best, baselineRefractiveIndex, alphaBestDeg, snapshotTimestamp));
    mphsave(model, snapshotFilePath);
    logline('Snapshot saved: %s\n', snapshotFilePath);
end
if ~isempty(coarseResultsTable), writetable(coarseResultsTable, fullfile(projectRootDir,'sensitivity_4d_coarse.csv')); end
if ~isempty(fineResultsTable),   writetable(fineResultsTable,   fullfile(projectRootDir,'sensitivity_4d_fine.csv')); end
if ~isempty(superResultsTable),  writetable(superResultsTable,  fullfile(projectRootDir,'sensitivity_4d_super.csv')); end
writetable(bestDenseTable, fullfile(projectRootDir,'sensitivity_4d_bestSensitivity_dense_ALLn.csv'));
if isfile(checkpointFilePath), delete(checkpointFilePath); end
if MAKE_PLOTS
    figure('Name','(1) TMOKE(alpha) for each n (Stage 4)','NumberTitle','off','Color','w');
    hold on; grid on;
    for i = 1:numel(validationRefractiveIndexList)
        plot(alphaGridsByN{i}, tmokeCurvesByN{i}, 'LineWidth', 1.2, 'DisplayName', sprintf('n=%.2f', validationRefractiveIndexList(i)));
    end
    xlabel('\alpha [deg]'); ylabel('TMOKE(\alpha)');
    title(sprintf('Best sensitivity candidate | TMOKE(\\alpha) vs \\alpha (S_{dense}=%.4f deg/RIU)', sensitivityDense));
    legend('Location','best');
    figure('Name','(2) alpha_peak vs n (Stage 4)','NumberTitle','off','Color','w');
    hold on; grid on;
    plot(validationRefractiveIndexList, alphaPeakDegreesByN, 'o', 'LineWidth', 1.5, 'DisplayName','data');
    nfit = linspace(min(validationRefractiveIndexList), max(validationRefractiveIndexList), 100);
    plot(nfit, polyval(alphaVsNLinearFit,nfit), '-', 'LineWidth', 1.2, 'DisplayName', sprintf('fit (S=%.4f)', sensitivityDense));
    xlabel('n'); ylabel('\alpha_{peak} [deg]');
    title('alpha_{peak} as a function of n');
    legend('Location','best');
    figure('Name','(3) |TMOKE| at tracked peak vs n (Stage 4)','NumberTitle','off','Color','w');
    grid on; hold on;
    plot(validationRefractiveIndexList, trackedTmokeAbsByN, 'o-', 'LineWidth', 1.5);
    xlabel('n'); ylabel('|TMOKE| at tracked peak');
    title('|TMOKE| at the tracked peak vs n');
end
if SAVE_FIGS
    saveAllOpenFigures(figureOutputDirectory, "sensitivity_4d", FIG_FORMATS);
    fprintf('Figures saved in: %s\n', figureOutputDirectory);
end
elapsed_total = toc(globalTimerStart);
fprintf('\n===== SUMMARY =====\n');
fprintf('Best |sensitivityEstimateFast| only (Stage 3):\n'); disp(bestSensitivityCandidate);
fprintf('Stage 4 @ bestSensitivityCandidate: sensitivityDense ~= %.6f deg/RIU\n', sensitivityDense);
fprintf('Runs done: %d | Elapsed: %s\n', runsCompletedGlobal, fmt_time_long(elapsed_total));
diary off;
function setParamNm(mdl, name, val_nm)
    mdl.param.set(name, sprintf('%.12g[nm]', val_nm));
end
function setParamScalar(mdl, name, val)
    mdl.param.set(name, sprintf('%.12g', val));
end
function Tsel = selectTopK_single_abs(T, K, col)
    [~, ord] = sort(abs(T.(col)), 'descend');
    Tsel = T(ord(1:min(K,height(T))), :);
end
function trackedEval = evaluateTrackedSensitivityAndCurves( ...
        mdl, studyTag, PARAM_N, refractiveIndexList, ...
        trackingReferenceRefractiveIndex, trackingHalfWindowDeg, ...
        alphaName, mName, referenceSweepStartDeg, referenceSweepStepDeg, referenceSweepStopDeg, ...
        globalAlphaMinDeg, globalAlphaMaxDeg, sweepAllNonReferenceCurves, ...
        mPlusStr, mMinusStr, ttagPlus, ttagMinus)
    refractiveIndexList = refractiveIndexList(:).';
    tol = 1e-9;
    referenceIndex = find(abs(refractiveIndexList - trackingReferenceRefractiveIndex) < tol, 1, 'first');
    if isempty(referenceIndex)
        error('Reference refractive index %.4f is not present in the configured list.', ...
            trackingReferenceRefractiveIndex);
    end
    [~, highestNIndex] = max(refractiveIndexList);
    trackedEval = struct();
    trackedEval.referenceIndex = referenceIndex;
    trackedEval.highestNIndex = highestNIndex;
    trackedEval.alphaPeakDegreesByN = zeros(size(refractiveIndexList));
    trackedEval.tmokeAtTrackedPeakByN = zeros(size(refractiveIndexList));
    trackedEval.trackedTmokeAbsByN = zeros(size(refractiveIndexList));
    trackedEval.alphaGridsByN = cell(size(refractiveIndexList));
    trackedEval.transmissionPlusByN = cell(size(refractiveIndexList));
    trackedEval.transmissionMinusByN = cell(size(refractiveIndexList));
    trackedEval.tmokeCurvesByN = cell(size(refractiveIndexList));
    evaluationOrder = buildTrackingEvaluationOrder(refractiveIndexList, trackingReferenceRefractiveIndex);
    referencePeakAlphaDeg = NaN;
    referencePeakSign = 0;
    for orderPos = 1:numel(evaluationOrder)
        idx = evaluationOrder(orderPos);
        nVal = refractiveIndexList(idx);
        if idx == referenceIndex || sweepAllNonReferenceCurves
            sweepStartDeg = referenceSweepStartDeg;
            sweepStopDeg = referenceSweepStopDeg;
        else
            sweepStartDeg = max(globalAlphaMinDeg, referencePeakAlphaDeg - trackingHalfWindowDeg);
            sweepStopDeg  = min(globalAlphaMaxDeg, referencePeakAlphaDeg + trackingHalfWindowDeg);
        end
        setParamScalar(mdl, PARAM_N, nVal);
        [alphaDeg, transmissionPlus, transmissionMinus, tmokeCurve] = solveAndGetTplusTminus( ...
            mdl, studyTag, alphaName, mName, ...
            sweepStartDeg, referenceSweepStepDeg, sweepStopDeg, ...
            mPlusStr, mMinusStr, ttagPlus, ttagMinus);
        trackedEval.alphaGridsByN{idx} = alphaDeg;
        trackedEval.transmissionPlusByN{idx} = transmissionPlus;
        trackedEval.transmissionMinusByN{idx} = transmissionMinus;
        trackedEval.tmokeCurvesByN{idx} = tmokeCurve;
        if idx == referenceIndex
            [~, trackedPeakIndex] = max(abs(tmokeCurve));
            referencePeakAlphaDeg = alphaDeg(trackedPeakIndex);
            referencePeakSign = sign(tmokeCurve(trackedPeakIndex));
        else
            trackedPeakIndex = selectTrackedPeakIndex( ...
                alphaDeg, tmokeCurve, referencePeakAlphaDeg, referencePeakSign, ...
                trackingHalfWindowDeg, globalAlphaMinDeg, globalAlphaMaxDeg);
        end
        trackedEval.alphaPeakDegreesByN(idx) = alphaDeg(trackedPeakIndex);
        trackedEval.tmokeAtTrackedPeakByN(idx) = tmokeCurve(trackedPeakIndex);
        trackedEval.trackedTmokeAbsByN(idx) = abs(tmokeCurve(trackedPeakIndex));
    end
    trackedEval.referencePeakAlphaDeg = referencePeakAlphaDeg;
    trackedEval.referencePeakSign = referencePeakSign;
    trackedEval.linearFit = polyfit(refractiveIndexList, trackedEval.alphaPeakDegreesByN, 1);
    trackedEval.sensitivitySlope = trackedEval.linearFit(1);
end
function evaluationOrder = buildTrackingEvaluationOrder(refractiveIndexList, trackingReferenceRefractiveIndex)
    tol = 1e-9;
    idxReference = find(abs(refractiveIndexList - trackingReferenceRefractiveIndex) < tol, 1, 'first');
    idxLower = find(refractiveIndexList < trackingReferenceRefractiveIndex - tol);
    idxUpper = find(refractiveIndexList > trackingReferenceRefractiveIndex + tol);
    [~, lowerOrder] = sort(refractiveIndexList(idxLower), 'descend');
    [~, upperOrder] = sort(refractiveIndexList(idxUpper), 'ascend');
    evaluationOrder = [idxReference, idxLower(lowerOrder), idxUpper(upperOrder)];
end
function trackedPeakIndex = selectTrackedPeakIndex( ...
        alphaDeg, tmokeCurve, referencePeakAlphaDeg, referencePeakSign, ...
        trackingHalfWindowDeg, globalAlphaMinDeg, globalAlphaMaxDeg)
    trackedWindowStartDeg = max(globalAlphaMinDeg, referencePeakAlphaDeg - trackingHalfWindowDeg);
    trackedWindowStopDeg = min(globalAlphaMaxDeg, referencePeakAlphaDeg + trackingHalfWindowDeg);
    inTrackedWindow = alphaDeg >= trackedWindowStartDeg - 1e-9 & alphaDeg <= trackedWindowStopDeg + 1e-9;
    if ~any(inTrackedWindow)
        error('The tracked TMOKE window [%.3f, %.3f] contains no alpha points.', ...
            trackedWindowStartDeg, trackedWindowStopDeg);
    end
    candidateIndices = find(inTrackedWindow);
    if referencePeakSign ~= 0
        sameSignMask = sign(tmokeCurve(candidateIndices)) == referencePeakSign;
        if any(sameSignMask)
            candidateIndices = candidateIndices(sameSignMask);
        end
    end
    [~, localIndex] = max(abs(tmokeCurve(candidateIndices)));
    trackedPeakIndex = candidateIndices(localIndex);
end
function [alpha_deg, Tplus, Tminus, TMOKE] = solveAndGetTplusTminus( ...
    mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, ...
    mPlusStr, mMinusStr, ttagPlus, ttagMinus)
    ensureTable(mdl, ttagPlus);
    ensureTable(mdl, ttagMinus);
    Npts = 1 + floor((aStopDeg - aStartDeg)/aStepDeg + 1e-9);
    redirectAllNumericalsToTable(mdl, ttagPlus);
    clearTable(mdl, ttagPlus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mPlusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha1_deg, T1] = readAlphaAndTFromNamedTable(mdl, ttagPlus, Npts);
    redirectAllNumericalsToTable(mdl, ttagMinus);
    clearTable(mdl, ttagMinus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mMinusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha2_deg, T2] = readAlphaAndTFromNamedTable(mdl, ttagMinus, Npts);
    assert(numel(alpha1_deg)==Npts && numel(alpha2_deg)==Npts, 'Unexpected alpha sweep length.');
    assert(max(abs(alpha1_deg - alpha2_deg)) < 1e-9, 'alpha grids differ between m=+1 and m=-1.');
    alpha_deg = alpha1_deg;
    Tplus  = T1;
    Tminus = T2;
    denom = Tplus + Tminus;
    denom(abs(denom) < 1e-9) = 1e-9;
    TMOKE = 2 * (Tplus - Tminus) ./ denom;
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
function [alpha_deg, Tcol] = readAlphaAndTFromNamedTable(mdl, ttag, Npts)
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
        a = S.data(:,1); T = S.data(:,2);
    else
        hlow = lower(heads);
        aIdx = find(contains(hlow,'alpha'), 1, 'first');
        tIdx = find(contains(hlow,'transmission') | contains(hlow,'transmittance') | contains(hlow,'total transmission') | contains(hlow,'total transmittance') | contains(hlow,' total t'), 1, 'first');
        if isempty(aIdx) || isempty(tIdx)
            if ~isfield(S,'data') || isempty(S.data) || size(S.data,2) < 2
                error('Alpha/Transmission not found in %s and not enough columns.', ttag);
            end
            aIdx = 1; tIdx = 2;
        end
        a = S.data(:, aIdx); T = S.data(:, tIdx);
    end
    assert(numel(a) >= Npts && numel(T) >= Npts, ...
        'Table %s has %d rows; expected >= %d.', ttag, numel(a), Npts);
    a = a(end-Npts+1:end);
    T = T(end-Npts+1:end);
    if ~isempty(heads)
        hlow = lower(heads);
        if any(contains(hlow, '(rad)')) || any(contains(hlow, '[rad]')), a = a * 180/pi; end
    end
    alpha_deg = a;
    Tcol      = T;
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
    checkpointData.schema           = 'tracked_same_peak_sensitivity_v1';
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
        coarseResultsTable, fineResultsTable, superResultsTable, bestDenseTable)
    points_since_ckpt_new = pointsSinceCheckpoint;
    if pointsSinceCheckpoint >= every_points
        save_checkpoint(checkpointFilePath, stage, runsCompletedGlobal, done_points, payload);
        try
            write_progress_xlsx(progressWorkbookPath, ...
                'coarse', coarseResultsTable, 'fine', fineResultsTable, 'super', superResultsTable, ...
                'dense', bestDenseTable);
        catch
        end
        points_since_ckpt_new = 0;
    end
end
function updateLivePlot(stage, Ldom, Lden, hsi, hau, ...
                        alpha_deg, TMOKE_vec, alpha_peak, tmoke_peak, ...
                        Tplus, Tminus, saveDir)
    persistent fig ax hTM hPeak hRp hRm inited frameCounters
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
        hRp = plot(ax, nan, nan, '--', 'LineWidth', 1.0, 'DisplayName','T^+(\alpha)');
        hRm = plot(ax, nan, nan, ':',  'LineWidth', 1.0, 'DisplayName','T^-(\alpha)');
        ylabel(ax,'Transmission');
        legend(ax,'Location','best');
        inited = true;
    end
    yyaxis(ax,'left');
    set(hTM,   'XData', alpha_deg, 'YData', TMOKE_vec);
    set(hPeak, 'XData', alpha_peak, 'YData', tmoke_peak);
    yyaxis(ax,'right');
    if nargin >= 11 && ~isempty(Tplus) && ~isempty(Tminus)
        set(hRp, 'XData', alpha_deg, 'YData', Tplus);
        set(hRm, 'XData', alpha_deg, 'YData', Tminus);
    else
        set(hRp, 'XData', nan, 'YData', nan);
        set(hRm, 'XData', nan, 'YData', nan);
    end
    title(ax, sprintf('%s | Ldom=%g | Lden=%g | hsi=%g | hau=%.3f | |TM|=%.5f @ %.3f deg', ...
        stage, Ldom, Lden, hsi, hau, abs(tmoke_peak), alpha_peak));
    drawnow('limitrate');
    if nargin >= 12 && ~isempty(saveDir)
        if isempty(frameCounters)
            frameCounters = containers.Map('KeyType','char','ValueType','double');
        end
        if ~isKey(frameCounters, stage), frameCounters(stage) = 0; end
        frameCounters(stage) = frameCounters(stage) + 1;
        iterIndex = frameCounters(stage);
        stageDir = fullfile(saveDir, stage);
        if ~exist(stageDir,'dir'), mkdir(stageDir); end
        baseName = sprintf('%s_%05d_Ldom%g_Lden%g_hsi%g_hau%.3f', ...
            stage, iterIndex, Ldom, Lden, hsi, hau);
        try
            exportgraphics(fig, fullfile(stageDir, [baseName '.png']), 'Resolution', 150);
        catch ME
            warning('Failed to save iteration plot %s: %s', baseName, ME.message);
        end
    end
end
function saveStageCandidatePlots(T, stageName, outDir, formats)
    if isempty(T), return; end
    if ~exist(outDir,'dir'), mkdir(outDir); end
    fig = figure('Name', sprintf('%s candidate sensitivity', stageName), ...
                 'NumberTitle','off','Color','w','Visible','off');
    candidateIndex = 1:height(T);
    plot(candidateIndex, abs(T.S_est_deg_per_RIU), '-o', 'LineWidth', 1.1);
    xlabel('Candidate index');
    ylabel('|S_{est}| [deg/RIU]');
    title(sprintf('%s sensitivity by candidate', stageName));
    grid on;
    saveFigureAndClose(fig, outDir, sprintf('%s_candidate_sensitivity', lower(stageName)), formats);
end
function saveValidationPhasePlots(nList, alphaGridsByN, tmokeCurvesByN, ...
        alphaPeakDegreesByN, trackedTmokeAbsByN, alphaVsNLinearFit, sensitivityDense, ...
        outDir, formats)
    if ~exist(outDir,'dir'), mkdir(outDir); end
    fig = figure('Name','Stage 4 dense tracked curves', ...
                 'NumberTitle','off','Color','w','Visible','off');
    tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');
    nexttile;
    hold on; grid on;
    for i = 1:numel(nList)
        plot(alphaGridsByN{i}, tmokeCurvesByN{i}, 'LineWidth', 1.1, ...
            'DisplayName', sprintf('n=%.2f', nList(i)));
    end
    xlabel('\alpha [deg]');
    ylabel('TMOKE');
    title('Stage 4 TMOKE curves');
    legend('Location','best');
    nexttile;
    hold on; grid on;
    plot(nList, alphaPeakDegreesByN, 'o', 'LineWidth', 1.4, 'DisplayName','data');
    nfit = linspace(min(nList), max(nList), 100);
    plot(nfit, polyval(alphaVsNLinearFit,nfit), '-', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('S=%.4f', sensitivityDense));
    xlabel('n');
    ylabel('\alpha_{peak} [deg]');
    title('Tracked peak shift');
    legend('Location','best');
    nexttile;
    plot(nList, trackedTmokeAbsByN, 'o-', 'LineWidth', 1.4);
    xlabel('n');
    ylabel('|TMOKE|');
    title('|TMOKE| at tracked peak');
    grid on;
    saveFigureAndClose(fig, outDir, 'valid_dense_tracked_curves', formats);
end
function saveFigureAndClose(fig, outDir, baseName, formats)
    for k = 1:numel(formats)
        ext = lower(string(formats{k}));
        filePath = fullfile(outDir, string(baseName) + "." + ext);
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
    close(fig);
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
