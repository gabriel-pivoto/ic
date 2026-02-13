function [superResultsTable, bestTradeoffCandidate, bestTmokeCandidate, bestSensitivityCandidate, runState] = runSuperStage(model, cfg, tags, paths, plan, runState, resumeInfo, coarseResultsTable, fineResultsTable, superSeedCandidates)
% Final refinement around fine seeds with tighter alpha steps.
import tmoke.comsol.*;
import tmoke.selection.*;
import tmoke.checkpoint.*;
import tmoke.plot.*;
import tmoke.util.*;

resumeFromCheckpoint = resumeInfo.resumeFromCheckpoint;
resumeStageTag       = resumeInfo.resumeStageTag;
checkpointData       = resumeInfo.checkpointData;

runsPerSearchPoint = plan.runsPerSearchPoint;
superTotalRuns     = plan.superTotalRuns;

runsCompletedGlobal     = runState.runsCompletedGlobal;
pointsSinceCheckpoint   = runState.pointsSinceCheckpoint;
globalRunTargetEstimate = runState.globalRunTargetEstimate;
isTotalRunEstimateExact = runState.isTotalRunEstimateExact;
globalTimerStart        = runState.globalTimerStart;
checkpointEveryPoints   = runState.checkpointEveryPoints;

superGoldHeightDelta    = cfg.superGoldHeightDelta;    superGoldHeightStep    = cfg.superGoldHeightStep;
superCeyigHeightDelta   = cfg.superCeyigHeightDelta;   superCeyigHeightStep   = cfg.superCeyigHeightStep;
superDomainPeriodDelta  = cfg.superDomainPeriodDelta;  superDomainPeriodStep  = cfg.superDomainPeriodStep;
superToothWidthDelta    = cfg.superToothWidthDelta;    superToothWidthStep    = cfg.superToothWidthStep;
superSiliconHeightDelta = cfg.superSiliconHeightDelta; superSiliconHeightStep = cfg.superSiliconHeightStep;

superAlphaHalfSpanSensitivity = cfg.superAlphaHalfSpanSensitivity;
alphaSuperStep = cfg.alphaSuperStep;
fastRefractiveIndexSamples = cfg.fastRefractiveIndexSamples;
PLOT_LIVE = cfg.PLOT_LIVE;

% --- Planning (exact) ---
if isempty(superSeedCandidates) && resumeFromCheckpoint && resumeStageTag=="SUPER"
    if isfield(checkpointData.payload,'superSeedCandidates')
        superSeedCandidates = checkpointData.payload.superSeedCandidates;
    else
        error('Checkpoint at SUPER stage has no "superSeedCandidates". Rerun FINE.');
    end
end

superTotalPoints = 0;
for s = 1:height(superSeedCandidates)
    goldHeightList    = (superSeedCandidates.h_au_nm(s)    - superGoldHeightDelta)  : superGoldHeightStep  : (superSeedCandidates.h_au_nm(s)    + superGoldHeightDelta);
    ceyigHeightList   = (superSeedCandidates.h_ceyig_nm(s) - superCeyigHeightDelta) : superCeyigHeightStep : (superSeedCandidates.h_ceyig_nm(s) + superCeyigHeightDelta);
    domainPeriodList  = (superSeedCandidates.L_domain_nm(s)- superDomainPeriodDelta) : superDomainPeriodStep : (superSeedCandidates.L_domain_nm(s)+ superDomainPeriodDelta);
    toothWidthList    = (superSeedCandidates.l_dente_nm(s) - superToothWidthDelta) : superToothWidthStep : (superSeedCandidates.l_dente_nm(s) + superToothWidthDelta);
    siliconHeightList = (superSeedCandidates.h_si_nm(s)    - superSiliconHeightDelta)  : superSiliconHeightStep  : (superSeedCandidates.h_si_nm(s)    + superSiliconHeightDelta);
    superTotalPoints = superTotalPoints + numel(domainPeriodList)*numel(toothWidthList)*numel(siliconHeightList)*numel(ceyigHeightList)*numel(goldHeightList);
end
superTotalRuns = runsPerSearchPoint * superTotalPoints;

fprintf('\nSTAGE SUPER - EXACT: %d runs (TOP-%d)\n', superTotalRuns, height(superSeedCandidates));

% --- Sweep ---
superResultsTable = [];
bestTradeoffCandidate   = [];
bestTmokeCandidate      = [];
bestSensitivityCandidate = [];

if resumeFromCheckpoint && any(strcmp(resumeStageTag, ["VALID","FULL"]))
    if isfield(checkpointData.payload,'superResultsTable'), superResultsTable = checkpointData.payload.superResultsTable; end
    if isfield(checkpointData.payload,'bestTradeoffCandidate'),   bestTradeoffCandidate   = checkpointData.payload.bestTradeoffCandidate;   end
    if isfield(checkpointData.payload,'bestTmokeCandidate'),      bestTmokeCandidate      = checkpointData.payload.bestTmokeCandidate;      end
    if isfield(checkpointData.payload,'bestSensitivityCandidate'), bestSensitivityCandidate = checkpointData.payload.bestSensitivityCandidate; end
    fprintf('SKIP SUPER -> restored superResultsTable (rows=%d) and best selections.\n', size(superResultsTable,1));
else
    superRows = [];
    stageRunsStart  = runsCompletedGlobal;
    stageTotalRuns  = superTotalRuns;
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
        alphaWindowCenterDeg = superSeedCandidates.alpha_peak_base_deg(s);

        goldHeightList    = (superSeedCandidates.h_au_nm(s)    - superGoldHeightDelta)  : superGoldHeightStep  : (superSeedCandidates.h_au_nm(s)    + superGoldHeightDelta);
        ceyigHeightList   = (superSeedCandidates.h_ceyig_nm(s) - superCeyigHeightDelta) : superCeyigHeightStep : (superSeedCandidates.h_ceyig_nm(s) + superCeyigHeightDelta);
        domainPeriodList  = (superSeedCandidates.L_domain_nm(s)- superDomainPeriodDelta) : superDomainPeriodStep : (superSeedCandidates.L_domain_nm(s)+ superDomainPeriodDelta);
        toothWidthList    = (superSeedCandidates.l_dente_nm(s) - superToothWidthDelta) : superToothWidthStep : (superSeedCandidates.l_dente_nm(s) + superToothWidthDelta);
        siliconHeightList = (superSeedCandidates.h_si_nm(s)    - superSiliconHeightDelta)  : superSiliconHeightStep  : (superSeedCandidates.h_si_nm(s)    + superSiliconHeightDelta);

        alphaStartDeg = max(0,  alphaWindowCenterDeg - superAlphaHalfSpanSensitivity);
        alphaStopDeg  = min(89, alphaWindowCenterDeg + superAlphaHalfSpanSensitivity);

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

                            superPointIndex = superPointIndex + 1;
                            if resumeFromCheckpoint && resumeStageTag=="SUPER" && superPointIndex <= checkpointData.done_points
                                continue;
                            end

                            setParamScalar(model, tags.PARAM_N, fastRefractiveIndexSamples(1));
                            [alphaGridFastN1, reflectancePlusFastN1, reflectanceMinusFastN1, tmokeCurveFastN1] = solveAndGetRplusRminus( ...
                                model, tags.STUDY_TAG, tags.ALPHA_NAME, tags.MSIGN_NAME, ...
                                alphaStartDeg, alphaSuperStep, alphaStopDeg, ...
                                tags.M_PLUS, tags.M_MINUS, tags.RPLUS_TABLE_TAG, tags.RMINUS_TABLE_TAG);

                            [tmokePeakMagnitudeFastN1, tmokePeakIndexFastN1] = max(abs(tmokeCurveFastN1));
                            alphaAtPeakFastN1 = alphaGridFastN1(tmokePeakIndexFastN1);
                            tmokeAtPeakFastN1 = tmokeCurveFastN1(tmokePeakIndexFastN1);

                            setParamScalar(model, tags.PARAM_N, fastRefractiveIndexSamples(2));
                            [alphaGridFastN2, ~, ~, tmokeCurveFastN2] = solveAndGetRplusRminus( ...
                                model, tags.STUDY_TAG, tags.ALPHA_NAME, tags.MSIGN_NAME, ...
                                alphaStartDeg, alphaSuperStep, alphaStopDeg, ...
                                tags.M_PLUS, tags.M_MINUS, tags.RPLUS_TABLE_TAG, tags.RMINUS_TABLE_TAG);

                            [~, tmokePeakIndexFastN2] = max(abs(tmokeCurveFastN2));
                            alphaAtPeakFastN2 = alphaGridFastN2(tmokePeakIndexFastN2);

                            sensitivityEstimateFast = (alphaAtPeakFastN2 - alphaAtPeakFastN1) / (fastRefractiveIndexSamples(2) - fastRefractiveIndexSamples(1));

                            if PLOT_LIVE
                                updateLivePlot('SUPER', domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), ceyigHeightList(ceyigHeightIdx), goldHeightList(goldHeightIdx), ...
                                    alphaGridFastN1, tmokeCurveFastN1, alphaAtPeakFastN1, tmokeAtPeakFastN1, reflectancePlusFastN1, reflectanceMinusFastN1);
                            end

                            superRows = [superRows; ...
                                domainPeriodList(domainPeriodIdx), toothWidthList(toothWidthIdx), siliconHeightList(siliconHeightIdx), ...
                                ceyigHeightList(ceyigHeightIdx), goldHeightList(goldHeightIdx), ...
                                tmokePeakMagnitudeFastN1, alphaAtPeakFastN1, tmokeAtPeakFastN1, ...
                                alphaAtPeakFastN2, sensitivityEstimateFast]; %#ok<AGROW>

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
                                'SUPER', checkpointEveryPoints, paths.checkpointFilePath, paths.progressWorkbookPath, ...
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

    bestTradeoffCandidate    = selectTopK_tradeoff(superResultsTable, 1, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');
    bestTmokeCandidate       = selectTopK_single(superResultsTable, 1, 'maxAbsTMOKE_base');
    bestSensitivityCandidate = selectTopK_single_abs(superResultsTable, 1, 'S_est_deg_per_RIU');

    fprintf('\n===== BEST (SUPERFINE) =====\n');
    fprintf('Trade-off (rankTM+rankS):\n'); disp(bestTradeoffCandidate);
    fprintf('Best TMOKE only:\n');         disp(bestTmokeCandidate);
    fprintf('Best |S| only:\n');           disp(bestSensitivityCandidate);

    save_checkpoint(paths.checkpointFilePath, 'VALID', runsCompletedGlobal, 0, struct( ...
        'coarseResultsTable', coarseResultsTable, 'fineResultsTable', fineResultsTable, 'superResultsTable', superResultsTable, ...
        'bestTradeoffCandidate', bestTradeoffCandidate, 'bestTmokeCandidate', bestTmokeCandidate, 'bestSensitivityCandidate', bestSensitivityCandidate));
    write_progress_xlsx(paths.progressWorkbookPath, 'super', superResultsTable);
end

runState.runsCompletedGlobal   = runsCompletedGlobal;
runState.pointsSinceCheckpoint = pointsSinceCheckpoint;
runState.globalRunTargetEstimate = globalRunTargetEstimate;
runState.isTotalRunEstimateExact = isTotalRunEstimateExact;
runState.globalTimerStart = globalTimerStart;
runState.checkpointEveryPoints = checkpointEveryPoints;

end
