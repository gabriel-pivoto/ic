function [validation, runState, bestTradeoffCandidate, bestTmokeCandidate, bestSensitivityCandidate] = runValidationStage(model, cfg, tags, paths, runState, resumeInfo, coarseResultsTable, fineResultsTable, superResultsTable, bestTradeoffCandidate, bestTmokeCandidate, bestSensitivityCandidate)
% Dense alpha sweeps at multiple n values for the chosen trade-off candidate.
import tmoke.comsol.*;
import tmoke.checkpoint.*;
import tmoke.util.*;

resumeFromCheckpoint = resumeInfo.resumeFromCheckpoint;
resumeStageTag       = resumeInfo.resumeStageTag;
checkpointData       = resumeInfo.checkpointData;

runsCompletedGlobal = runState.runsCompletedGlobal;

baselineRefractiveIndex     = cfg.baselineRefractiveIndex;
validationRefractiveIndexList = cfg.validationRefractiveIndexList;
alphaDenseStep              = cfg.alphaDenseStep;

% Restore bestTradeoffCandidate if resuming from VALID/FULL
if isempty(bestTradeoffCandidate) && resumeFromCheckpoint && any(strcmp(resumeStageTag, ["VALID","FULL"]))
    bestTradeoffCandidate = checkpointData.payload.bestTradeoffCandidate;
end

Ldom_best = bestTradeoffCandidate.L_domain_nm(1);
Lden_best = bestTradeoffCandidate.l_dente_nm(1);
hsi_best  = bestTradeoffCandidate.h_si_nm(1);
hcey_best = bestTradeoffCandidate.h_ceyig_nm(1);
hau_best  = bestTradeoffCandidate.h_au_nm(1);

setParamNm(model, tags.PARAM_LDOM, Ldom_best);
setParamNm(model, tags.PARAM_LDEN, Lden_best);
setParamNm(model, tags.PARAM_HSI,  hsi_best);
setParamNm(model, tags.PARAM_HCEY, hcey_best);
setParamNm(model, tags.PARAM_HAU,  hau_best);

alphaPeakDegreesByN = zeros(size(validationRefractiveIndexList));
tmokeMaxAbsByN      = zeros(size(validationRefractiveIndexList));
tmokeCurvesByN      = cell(size(validationRefractiveIndexList));
alphaGridsByN       = cell(size(validationRefractiveIndexList));
reflectancePlusByN  = cell(size(validationRefractiveIndexList));
reflectanceMinusByN = cell(size(validationRefractiveIndexList));

for i = 1:numel(validationRefractiveIndexList)
    setParamScalar(model, tags.PARAM_N, validationRefractiveIndexList(i));

    [a, Rp, Rm, TM] = solveAndGetRplusRminus( ...
        model, tags.STUDY_TAG, tags.ALPHA_NAME, tags.MSIGN_NAME, ...
        0, alphaDenseStep, 89, ...
        tags.M_PLUS, tags.M_MINUS, tags.RPLUS_TABLE_TAG, tags.RMINUS_TABLE_TAG);

    [tmmax, k] = max(abs(TM));
    alphaPeakDegreesByN(i) = a(k);
    tmokeMaxAbsByN(i)      = tmmax;

    alphaGridsByN{i}       = a;
    tmokeCurvesByN{i}      = TM;
    reflectancePlusByN{i}  = Rp;
    reflectanceMinusByN{i} = Rm;

    runsCompletedGlobal = runsCompletedGlobal + 2;
end

alphaVsNLinearFit = polyfit(validationRefractiveIndexList, alphaPeakDegreesByN, 1);
sensitivityDense = alphaVsNLinearFit(1);

fprintf('\n===== VALID (bestTradeoffCandidate) =====\n');
fprintf('Best geom: Ldom=%g | Lden=%g | hsi=%g | hcey=%g | hau=%g\n', ...
    Ldom_best, Lden_best, hsi_best, hcey_best, hau_best);
fprintf('alpha_peak(n) = ['); fprintf(' %.4f', alphaPeakDegreesByN); fprintf(' ]\n');
fprintf('sensitivityDense ~= %.6f deg/RIU (linear fit)\n', sensitivityDense);

[~, idxBase] = min(abs(validationRefractiveIndexList - baselineRefractiveIndex));
alphaBestDeg = alphaPeakDegreesByN(idxBase);
[tmokeBaselineMax, tmokeBaselineIndex] = max(abs(tmokeCurvesByN{idxBase}));
tmokeBestValue = tmokeCurvesByN{idxBase}(tmokeBaselineIndex);

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

save_checkpoint(paths.checkpointFilePath, 'FULL', runsCompletedGlobal, 0, struct( ...
    'coarseResultsTable', coarseResultsTable, 'fineResultsTable', fineResultsTable, 'superResultsTable', superResultsTable, ...
    'bestTradeoffCandidate', bestTradeoffCandidate, 'bestTmokeCandidate', bestTmokeCandidate, 'bestSensitivityCandidate', bestSensitivityCandidate, ...
    'bestDenseTable', bestDenseTable, 'sensitivityDense', sensitivityDense, ...
    'alphaBestDeg', alphaBestDeg, 'tmokeBestValue', tmokeBestValue));
write_progress_xlsx(paths.progressWorkbookPath, 'dense', bestDenseTable);

validation = struct( ...
    'alphaPeakDegreesByN', alphaPeakDegreesByN, ...
    'tmokeMaxAbsByN', tmokeMaxAbsByN, ...
    'tmokeCurvesByN', {tmokeCurvesByN}, ...
    'alphaGridsByN', {alphaGridsByN}, ...
    'reflectancePlusByN', {reflectancePlusByN}, ...
    'reflectanceMinusByN', {reflectanceMinusByN}, ...
    'alphaVsNLinearFit', alphaVsNLinearFit, ...
    'sensitivityDense', sensitivityDense, ...
    'alphaBestDeg', alphaBestDeg, ...
    'tmokeBestValue', tmokeBestValue, ...
    'bestDenseTable', bestDenseTable);

runState.runsCompletedGlobal = runsCompletedGlobal;

end
