function [bestFullTable, runState] = runFullStage(model, cfg, tags, paths, runState, bestTradeoffCandidate, validation)
% Optional snapshot + final baseline sweep for the chosen geometry.
import tmoke.comsol.*;
import tmoke.plot.*;
import tmoke.util.*;

runsCompletedGlobal = runState.runsCompletedGlobal;

baselineRefractiveIndex = cfg.baselineRefractiveIndex;
alphaFullStep = cfg.alphaFullStep;
PLOT_LIVE    = cfg.PLOT_LIVE;
SAVE_SNAPSHOT = cfg.SAVE_SNAPSHOT;

Ldom_best = bestTradeoffCandidate.L_domain_nm(1);
Lden_best = bestTradeoffCandidate.l_dente_nm(1);
hsi_best  = bestTradeoffCandidate.h_si_nm(1);
hcey_best = bestTradeoffCandidate.h_ceyig_nm(1);
hau_best  = bestTradeoffCandidate.h_au_nm(1);

alphaBestDeg  = validation.alphaBestDeg;
tmokeBestValue = validation.tmokeBestValue;

setParamNm(model, tags.PARAM_LDOM, Ldom_best);
setParamNm(model, tags.PARAM_LDEN, Lden_best);
setParamNm(model, tags.PARAM_HSI,  hsi_best);
setParamNm(model, tags.PARAM_HCEY, hcey_best);
setParamNm(model, tags.PARAM_HAU,  hau_best);

if SAVE_SNAPSHOT
    setParamScalar(model, tags.PARAM_N, baselineRefractiveIndex);

    setAlphaMSweep(model, tags.STUDY_TAG, tags.ALPHA_NAME, alphaBestDeg, 0.01, alphaBestDeg, tags.MSIGN_NAME, sprintf('%s %s', tags.M_PLUS, tags.M_MINUS));
    model.study(tags.STUDY_TAG).run;
    refreshDerivedValues(model);

    runsCompletedGlobal = runsCompletedGlobal + 2;

    snapshotTimestamp = datestr(now,'yyyymmdd_HHMMSS');
    snapshotFilePath = fullfile(paths.projectRootDir, sprintf( ...
        'snapshot_bestTradeoff_Ldom%4.0f_Lden%4.0f_hsi%3.0f_hcey%.3f_hau%.3f_n%.3f_alpha%.4f_%s.mph', ...
        Ldom_best, Lden_best, hsi_best, hcey_best, hau_best, baselineRefractiveIndex, alphaBestDeg, snapshotTimestamp));
    mphsave(model, snapshotFilePath);
    logline('Snapshot saved: %s\n', snapshotFilePath);
end

setParamScalar(model, tags.PARAM_N, baselineRefractiveIndex);
[alphaFull, reflectancePlusFull, reflectanceMinusFull, tmokeFull] = solveAndGetRplusRminus( ...
    model, tags.STUDY_TAG, tags.ALPHA_NAME, tags.MSIGN_NAME, ...
    0, alphaFullStep, 89, ...
    tags.M_PLUS, tags.M_MINUS, tags.RPLUS_TABLE_TAG, tags.RMINUS_TABLE_TAG);

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

runState.runsCompletedGlobal = runsCompletedGlobal;

end
