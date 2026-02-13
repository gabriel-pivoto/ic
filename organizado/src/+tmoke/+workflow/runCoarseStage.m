function [coarseResultsTable, coarseSeedCandidates, runState] = runCoarseStage(model, cfg, tags, paths, plan, runState, resumeInfo)
% Sweep the coarse grid, compute TMOKE + fast sensitivity, and checkpoint.
import tmoke.comsol.*;
import tmoke.selection.*;
import tmoke.checkpoint.*;
import tmoke.plot.*;
import tmoke.util.*;

resumeFromCheckpoint = resumeInfo.resumeFromCheckpoint;
resumeStageTag       = resumeInfo.resumeStageTag;
checkpointData       = resumeInfo.checkpointData;

numDomainPeriodPoints  = plan.numDomainPeriodPoints;
numToothWidthPoints    = plan.numToothWidthPoints;
numSiliconHeightPoints = plan.numSiliconHeightPoints;
numCeyigHeightPoints   = plan.numCeyigHeightPoints;
numGoldHeightPoints    = plan.numGoldHeightPoints;
runsPerSearchPoint     = plan.runsPerSearchPoint;
coarseTotalRuns        = plan.coarseTotalRuns;

globalRunTargetEstimate = runState.globalRunTargetEstimate;
isTotalRunEstimateExact = runState.isTotalRunEstimateExact;
globalTimerStart        = runState.globalTimerStart;
checkpointEveryPoints   = runState.checkpointEveryPoints;
pointsSinceCheckpoint   = runState.pointsSinceCheckpoint;
runsCompletedGlobal     = runState.runsCompletedGlobal;

domainPeriodGridNm   = cfg.domainPeriodGridNm;
toothWidthGridNm     = cfg.toothWidthGridNm;
siliconHeightGridNm  = cfg.siliconHeightGridNm;
ceyigHeightGridNm    = cfg.ceyigHeightGridNm;
goldHeightGridNm     = cfg.goldHeightGridNm;
alphaCoarseRange     = cfg.alphaCoarseRange;
fastRefractiveIndexSamples = cfg.fastRefractiveIndexSamples;
PLOT_LIVE = cfg.PLOT_LIVE;

coarseResultsTable = [];
coarseSeedCandidates = [];

skipCoarseStage = resumeFromCheckpoint && any(strcmp(resumeStageTag, ["FINE","SUPER","VALID","FULL"]));
if skipCoarseStage
    if isfield(checkpointData.payload,'coarseResultsTable'), coarseResultsTable = checkpointData.payload.coarseResultsTable; end
    if isfield(checkpointData.payload,'coarseSeedCandidates'), coarseSeedCandidates = checkpointData.payload.coarseSeedCandidates; end
    fprintf('SKIP COARSE -> restored coarseResultsTable (rows=%d), coarseSeedCandidates (rows=%d)\n', ...
        size(coarseResultsTable,1), height(coarseSeedCandidates));
else
    fprintf('STAGE COARSE - EXACT: %d runs\n', coarseTotalRuns);
    stageRunsStart  = runsCompletedGlobal;
    stageTotalRuns  = coarseTotalRuns;
    stageTimerStart = tic;

    coarseRows = [];
    coarsePointIndex = 0;

    if resumeFromCheckpoint && resumeStageTag=="COARSE"
        runsCompletedGlobal = checkpointData.runsCompletedGlobal;
        if isfield(checkpointData.payload,'coarseRows'), coarseRows = checkpointData.payload.coarseRows; end
        coarsePointIndex = checkpointData.done_points;
        fprintf('>>> COARSE resume: skipping %d points already computed.\n', coarsePointIndex);
    end

    for domainPeriodIdx = 1:numDomainPeriodPoints
        setParamNm(model, tags.PARAM_LDOM, domainPeriodGridNm(domainPeriodIdx));
        for toothWidthIdx = 1:numToothWidthPoints
            setParamNm(model, tags.PARAM_LDEN, toothWidthGridNm(toothWidthIdx));
            for siliconHeightIdx = 1:numSiliconHeightPoints
                setParamNm(model, tags.PARAM_HSI, siliconHeightGridNm(siliconHeightIdx));
                for ceyigHeightIdx = 1:numCeyigHeightPoints
                    setParamNm(model, tags.PARAM_HCEY, ceyigHeightGridNm(ceyigHeightIdx));
                    for goldHeightIdx = 1:numGoldHeightPoints
                        setParamNm(model, tags.PARAM_HAU, goldHeightGridNm(goldHeightIdx));

                        coarsePointIndex = coarsePointIndex + 1;
                        if resumeFromCheckpoint && resumeStageTag=="COARSE" && coarsePointIndex <= checkpointData.done_points
                            continue;
                        end

                        % Base refractive index sweep (m = +/-1)
                        setParamScalar(model, tags.PARAM_N, fastRefractiveIndexSamples(1));
                        [alphaGridFastN1, reflectancePlusFastN1, reflectanceMinusFastN1, tmokeCurveFastN1] = solveAndGetRplusRminus( ...
                            model, tags.STUDY_TAG, tags.ALPHA_NAME, tags.MSIGN_NAME, ...
                            alphaCoarseRange(1), alphaCoarseRange(2), alphaCoarseRange(3), ...
                            tags.M_PLUS, tags.M_MINUS, tags.RPLUS_TABLE_TAG, tags.RMINUS_TABLE_TAG);

                        [tmokePeakMagnitudeFastN1, tmokePeakIndexFastN1] = max(abs(tmokeCurveFastN1));
                        alphaAtPeakFastN1 = alphaGridFastN1(tmokePeakIndexFastN1);
                        tmokeAtPeakFastN1  = tmokeCurveFastN1(tmokePeakIndexFastN1);

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

                        % Fast sensitivity estimate at second refractive index
                        setParamScalar(model, tags.PARAM_N, fastRefractiveIndexSamples(2));
                        [alphaGridFastN2, ~, ~, tmokeCurveFastN2] = solveAndGetRplusRminus( ...
                            model, tags.STUDY_TAG, tags.ALPHA_NAME, tags.MSIGN_NAME, ...
                            alphaCoarseRange(1), alphaCoarseRange(2), alphaCoarseRange(3), ...
                            tags.M_PLUS, tags.M_MINUS, tags.RPLUS_TABLE_TAG, tags.RMINUS_TABLE_TAG);

                        [~, tmokePeakIndexFastN2] = max(abs(tmokeCurveFastN2));
                        alphaAtPeakFastN2 = alphaGridFastN2(tmokePeakIndexFastN2);

                        sensitivityEstimateFast = (alphaAtPeakFastN2 - alphaAtPeakFastN1) / (fastRefractiveIndexSamples(2) - fastRefractiveIndexSamples(1));

                        coarseRows = [coarseRows; ...
                            domainPeriodGridNm(domainPeriodIdx), toothWidthGridNm(toothWidthIdx), siliconHeightGridNm(siliconHeightIdx), ...
                            ceyigHeightGridNm(ceyigHeightIdx), goldHeightGridNm(goldHeightIdx), ...
                            tmokePeakMagnitudeFastN1, alphaAtPeakFastN1, tmokeAtPeakFastN1, ...
                            alphaAtPeakFastN2, sensitivityEstimateFast]; %#ok<AGROW>

                        runsCompletedGlobal = runsCompletedGlobal + runsPerSearchPoint;

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

                        pointsSinceCheckpoint = pointsSinceCheckpoint + 1;
                        payload = struct('coarseRows',coarseRows);
                        pointsSinceCheckpoint = maybe_checkpoint( ...
                            'COARSE', checkpointEveryPoints, paths.checkpointFilePath, paths.progressWorkbookPath, ...
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

    coarseSeedCandidates = selectTopK_tradeoff(coarseResultsTable, cfg.topKCoarse, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');

    save_checkpoint(paths.checkpointFilePath, 'FINE', runsCompletedGlobal, 0, struct('coarseResultsTable', coarseResultsTable, 'coarseSeedCandidates', coarseSeedCandidates));
    write_progress_xlsx(paths.progressWorkbookPath, 'coarse', coarseResultsTable);
end

runState.runsCompletedGlobal   = runsCompletedGlobal;
runState.pointsSinceCheckpoint = pointsSinceCheckpoint;
runState.globalRunTargetEstimate = globalRunTargetEstimate;
runState.isTotalRunEstimateExact = isTotalRunEstimateExact;
runState.globalTimerStart = globalTimerStart;
runState.checkpointEveryPoints = checkpointEveryPoints;

end
