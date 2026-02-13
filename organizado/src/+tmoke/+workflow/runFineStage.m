function [fineResultsTable, superSeedCandidates, runState, plan] = runFineStage(model, cfg, tags, paths, plan, runState, resumeInfo, coarseResultsTable, coarseSeedCandidates)
% Refine around coarse seeds, making run budget exact and checkpointing.
import tmoke.comsol.*;
import tmoke.selection.*;
import tmoke.checkpoint.*;
import tmoke.plot.*;
import tmoke.util.*;

resumeFromCheckpoint = resumeInfo.resumeFromCheckpoint;
resumeStageTag       = resumeInfo.resumeStageTag;
checkpointData       = resumeInfo.checkpointData;

runsPerSearchPoint = plan.runsPerSearchPoint;
coarseTotalRuns    = plan.coarseTotalRuns;
superTotalRuns     = plan.superTotalRuns;
extraRunsFixed     = plan.extraRunsFixed;

runsCompletedGlobal     = runState.runsCompletedGlobal;
pointsSinceCheckpoint   = runState.pointsSinceCheckpoint;
globalRunTargetEstimate = runState.globalRunTargetEstimate;
isTotalRunEstimateExact = runState.isTotalRunEstimateExact;
globalTimerStart        = runState.globalTimerStart;
checkpointEveryPoints   = runState.checkpointEveryPoints;
MAX_RUNS                = runState.MAX_RUNS;

domainPeriodGridNm  = cfg.domainPeriodGridNm;
toothWidthGridNm    = cfg.toothWidthGridNm;
siliconHeightGridNm = cfg.siliconHeightGridNm;
ceyigHeightGridNm   = cfg.ceyigHeightGridNm;
goldHeightGridNm    = cfg.goldHeightGridNm;

fineGoldHeightDelta    = cfg.fineGoldHeightDelta;    fineGoldHeightStep    = cfg.fineGoldHeightStep;
fineCeyigHeightDelta   = cfg.fineCeyigHeightDelta;   fineCeyigHeightStep   = cfg.fineCeyigHeightStep;
fineDomainPeriodDelta  = cfg.fineDomainPeriodDelta;  fineDomainPeriodStep  = cfg.fineDomainPeriodStep;
fineToothWidthDelta    = cfg.fineToothWidthDelta;    fineToothWidthStep    = cfg.fineToothWidthStep;
fineSiliconHeightDelta = cfg.fineSiliconHeightDelta; fineSiliconHeightStep = cfg.fineSiliconHeightStep;

fineAlphaHalfSpanSensitivity = cfg.fineAlphaHalfSpanSensitivity;
alphaFineStep  = cfg.alphaFineStep;

fastRefractiveIndexSamples = cfg.fastRefractiveIndexSamples;
PLOT_LIVE = cfg.PLOT_LIVE;

% --- Exact planning for the fine stage (clamped windows) ---
if isempty(coarseSeedCandidates) && resumeFromCheckpoint && resumeStageTag=="FINE"
    if isfield(checkpointData.payload,'coarseSeedCandidates')
        coarseSeedCandidates = checkpointData.payload.coarseSeedCandidates;
    else
        error('Checkpoint at FINE stage has no "coarseSeedCandidates". Rerun COARSE.');
    end
end

fineTotalPoints = 0;
for s = 1:height(coarseSeedCandidates)
    goldHeightList    = max(min(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   - fineGoldHeightDelta)  : fineGoldHeightStep  : min(max(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   + fineGoldHeightDelta);
    ceyigHeightList   = max(min(ceyigHeightGridNm), coarseSeedCandidates.h_ceyig_nm(s)- fineCeyigHeightDelta) : fineCeyigHeightStep : min(max(ceyigHeightGridNm), coarseSeedCandidates.h_ceyig_nm(s)+ fineCeyigHeightDelta);
    domainPeriodList  = max(min(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)- fineDomainPeriodDelta) : fineDomainPeriodStep : min(max(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)+ fineDomainPeriodDelta);
    toothWidthList    = max(min(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) - fineToothWidthDelta) : fineToothWidthStep : min(max(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) + fineToothWidthDelta);
    siliconHeightList = max(min(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    - fineSiliconHeightDelta)  : fineSiliconHeightStep  : min(max(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    + fineSiliconHeightDelta);
    fineTotalPoints = fineTotalPoints + numel(domainPeriodList)*numel(toothWidthList)*numel(siliconHeightList)*numel(ceyigHeightList)*numel(goldHeightList);
end
fineTotalRuns = runsPerSearchPoint * fineTotalPoints;

globalRunTargetEstimate = coarseTotalRuns + fineTotalRuns + superTotalRuns + extraRunsFixed;
isTotalRunEstimateExact = true;

fprintf('\nSTAGE FINE - EXACT: %d runs (TOP-%d)\n', fineTotalRuns, height(coarseSeedCandidates));
fprintf('GLOBAL TOTAL (EXACT): %d runs\n\n', globalRunTargetEstimate);
if globalRunTargetEstimate > MAX_RUNS
    error('Planned total exceeds %d runs (%d). Adjust deltas/steps/TOPK.', MAX_RUNS, globalRunTargetEstimate);
end

% --- Sweep fine windows ---
fineResultsTable = [];
superSeedCandidates = [];

if resumeFromCheckpoint && any(strcmp(resumeStageTag, ["SUPER","VALID","FULL"]))
    if isfield(checkpointData.payload,'fineResultsTable'),      fineResultsTable = checkpointData.payload.fineResultsTable; end
    if isfield(checkpointData.payload,'superSeedCandidates'), superSeedCandidates = checkpointData.payload.superSeedCandidates; end
    fprintf('SKIP FINE -> restored fineResultsTable (rows=%d), superSeedCandidates (rows=%d)\n', ...
        size(fineResultsTable,1), height(superSeedCandidates));
else
    fineRows = [];
    stageRunsStart  = runsCompletedGlobal;
    stageTotalRuns  = fineTotalRuns;
    stageTimerStart = tic;

    finePointIndex = 0;
    if resumeFromCheckpoint && resumeStageTag=="FINE" && checkpointData.done_points > 0
        runsCompletedGlobal = checkpointData.runsCompletedGlobal;
        finePointIndex   = checkpointData.done_points;
        if isfield(checkpointData.payload,'fineRows'), fineRows = checkpointData.payload.fineRows; end
        if isfield(checkpointData.payload,'coarseSeedCandidates'), coarseSeedCandidates = checkpointData.payload.coarseSeedCandidates; end
        fprintf('>>> FINE resume: skipping %d points already computed.\n', finePointIndex);
    end

    for s = 1:height(coarseSeedCandidates)
        alphaWindowCenterDeg = coarseSeedCandidates.alpha_peak_base_deg(s);

        goldHeightList    = max(min(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   - fineGoldHeightDelta)  : fineGoldHeightStep  : min(max(goldHeightGridNm),  coarseSeedCandidates.h_au_nm(s)   + fineGoldHeightDelta);
        ceyigHeightList   = max(min(ceyigHeightGridNm), coarseSeedCandidates.h_ceyig_nm(s)- fineCeyigHeightDelta) : fineCeyigHeightStep : min(max(ceyigHeightGridNm), coarseSeedCandidates.h_ceyig_nm(s)+ fineCeyigHeightDelta);
        domainPeriodList  = max(min(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)- fineDomainPeriodDelta) : fineDomainPeriodStep : min(max(domainPeriodGridNm), coarseSeedCandidates.L_domain_nm(s)+ fineDomainPeriodDelta);
        toothWidthList    = max(min(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) - fineToothWidthDelta) : fineToothWidthStep : min(max(toothWidthGridNm),  coarseSeedCandidates.l_dente_nm(s) + fineToothWidthDelta);
        siliconHeightList = max(min(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    - fineSiliconHeightDelta)  : fineSiliconHeightStep  : min(max(siliconHeightGridNm),     coarseSeedCandidates.h_si_nm(s)    + fineSiliconHeightDelta);

        alphaStartDeg = max(0,  alphaWindowCenterDeg - fineAlphaHalfSpanSensitivity);
        alphaStopDeg  = min(89, alphaWindowCenterDeg + fineAlphaHalfSpanSensitivity);

        for domainPeriodIdx = 1:numel(domainPeriodList)
            setParamNm(model, tags.PARAM_LDOM, domainPeriodList(domainPeriodIdx));
            for toothWidthIdx = 1:numel(toothWidthList)
                setParamNm(model, tags.PARAM_LDEN, toothWidthList(toothWidthIdx));
                for siliconHeightIdx = 1:numel(siliconHeightList)
                    setParamNm(model, tags.PARAM_HSI, siliconHeightList(siliconHeightIdx));
                    for ceyigHeightIdx = 1:numel(ceyigHeightList)
                        setParamNm(model, tags.PARAM_HCEY, ceyigHeightList(ceyigHeightIdx));
                        for goldHeightIdx = 1:numel(goldHeightList)
                            setParamNm(model, tags.PARAM_HAU, goldHeightList(goldHeightIdx));

                            finePointIndex = finePointIndex + 1;
                            if resumeFromCheckpoint && resumeStageTag=="FINE" && finePointIndex <= checkpointData.done_points
                                continue;
                            end

                            setParamScalar(model, tags.PARAM_N, fastRefractiveIndexSamples(1));
                            [alphaGridFastN1, reflectancePlusFastN1, reflectanceMinusFastN1, tmokeCurveFastN1] = solveAndGetRplusRminus( ...
                                model, tags.STUDY_TAG, tags.ALPHA_NAME, tags.MSIGN_NAME, ...
                                alphaStartDeg, alphaFineStep, alphaStopDeg, ...
                                tags.M_PLUS, tags.M_MINUS, tags.RPLUS_TABLE_TAG, tags.RMINUS_TABLE_TAG);

                            [tmokePeakMagnitudeFastN1, tmokePeakIndexFastN1] = max(abs(tmokeCurveFastN1));
                            alphaAtPeakFastN1 = alphaGridFastN1(tmokePeakIndexFastN1);
                            tmokeAtPeakFastN1 = tmokeCurveFastN1(tmokePeakIndexFastN1);

                            setParamScalar(model, tags.PARAM_N, fastRefractiveIndexSamples(2));
                            [alphaGridFastN2, ~, ~, tmokeCurveFastN2] = solveAndGetRplusRminus( ...
                                model, tags.STUDY_TAG, tags.ALPHA_NAME, tags.MSIGN_NAME, ...
                                alphaStartDeg, alphaFineStep, alphaStopDeg, ...
                                tags.M_PLUS, tags.M_MINUS, tags.RPLUS_TABLE_TAG, tags.RMINUS_TABLE_TAG);

                            [~, tmokePeakIndexFastN2] = max(abs(tmokeCurveFastN2));
                            alphaAtPeakFastN2 = alphaGridFastN2(tmokePeakIndexFastN2);

                            sensitivityEstimateFast = (alphaAtPeakFastN2 - alphaAtPeakFastN1) / (fastRefractiveIndexSamples(2) - fastRefractiveIndexSamples(1));

                            if PLOT_LIVE
                                updateLivePlot('FINE', domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), ceyigHeightList(ceyigHeightIdx), goldHeightList(goldHeightIdx), ...
                                    alphaGridFastN1, tmokeCurveFastN1, alphaAtPeakFastN1, tmokeAtPeakFastN1, reflectancePlusFastN1, reflectanceMinusFastN1);
                            end

                            fineRows = [fineRows; ...
                                domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), ...
                                ceyigHeightList(ceyigHeightIdx), goldHeightList(goldHeightIdx), ...
                                tmokePeakMagnitudeFastN1, alphaAtPeakFastN1, tmokeAtPeakFastN1, ...
                                alphaAtPeakFastN2, sensitivityEstimateFast]; %#ok<AGROW>

                            runsCompletedGlobal = runsCompletedGlobal + runsPerSearchPoint;
                            stageRunsCompleted = runsCompletedGlobal - stageRunsStart;
                            [stageCompletionFraction, stageEtaSeconds]   = frac_eta(stageRunsCompleted, stageTotalRuns, toc(stageTimerStart));
                            [globalCompletionFraction, globalEtaSeconds] = frac_eta(runsCompletedGlobal, globalRunTargetEstimate, toc(globalTimerStart));

                            logline(['FINE   | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                     ' | |TM|=%.5f @ alpha=%.2f deg | sensitivityEstimateFast=%.4f deg/RIU' ...
                                     ' [Stage %5.1f%% | ETA %s | Global %s%% | ETA %s]\n'], ...
                                domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), ceyigHeightList(ceyigHeightIdx), goldHeightList(goldHeightIdx), ...
                                tmokePeakMagnitudeFastN1, alphaAtPeakFastN1, sensitivityEstimateFast, ...
                                100*stageCompletionFraction, fmt_time_long(stageEtaSeconds), fmt_pct(globalCompletionFraction), fmt_eta_flag(globalEtaSeconds, isTotalRunEstimateExact));

                            pointsSinceCheckpoint = pointsSinceCheckpoint + 1;
                            payload = struct('fineRows',fineRows, 'coarseSeedCandidates', coarseSeedCandidates);
                            pointsSinceCheckpoint = maybe_checkpoint( ...
                                'FINE', checkpointEveryPoints, paths.checkpointFilePath, paths.progressWorkbookPath, ...
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

    superSeedCandidates = selectTopK_tradeoff(fineResultsTable, cfg.topKFine, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');

    save_checkpoint(paths.checkpointFilePath, 'SUPER', runsCompletedGlobal, 0, struct( ...
        'coarseResultsTable', coarseResultsTable, 'fineResultsTable', fineResultsTable, 'superSeedCandidates', superSeedCandidates));
    write_progress_xlsx(paths.progressWorkbookPath, 'fine', fineResultsTable);
end

plan.fineTotalPoints = fineTotalPoints;
plan.fineTotalRuns   = fineTotalRuns;
plan.globalRunTargetEstimate = globalRunTargetEstimate;
plan.isTotalRunEstimateExact = isTotalRunEstimateExact;

runState.runsCompletedGlobal   = runsCompletedGlobal;
runState.pointsSinceCheckpoint = pointsSinceCheckpoint;
runState.globalRunTargetEstimate = globalRunTargetEstimate;
runState.isTotalRunEstimateExact = isTotalRunEstimateExact;
runState.globalTimerStart = globalTimerStart;
runState.checkpointEveryPoints = checkpointEveryPoints;

end
