%% =====================================================================
% Angular sensitivity vs n: INITIAL geometry vs BEST geometry
% ---------------------------------------------------------------------
% What it does (no search, just the sensitivity metric on two geometries):
%   1) Loads the COMSOL model exactly as saved (.mph).
%   2) For the model's INITIAL geometry, tracks the same TMOKE resonance
%      across a list of refractive indices n, records alpha_peak(n), and
%      fits a line -> sensitivity S = d(alpha_peak)/dn [deg/RIU].
%   3) Sets the BEST geometry (from melhor_dimensao.mat, produced by
%      extrair_melhor_dimensao.m) and repeats the sensitivity measurement.
%   4) Plots alpha_peak(n) + linear fit for each geometry, plus an overlay
%      comparison, and saves figures + CSVs to the "graficos" folder.
%
% Sensitivity method (identical to senseAndTmoke_semCeyig.m):
%   - solve the reference curve at n = 1.33 and lock the TMOKE peak;
%   - solve the other n values inside alpha_ref +/- trackingHalfWindow;
%   - fit alpha_peak vs n; the slope is the sensitivity.
% =====================================================================
clear; clc; close all; format long;
import com.comsol.model.*;
import com.comsol.model.util.*;

%% --------------------------- Config --------------------------------
projectRootDir  = 'D:\Gabriel Pivoto\projetoIC';
comsolModelFile = fullfile(projectRootDir,'modelosimplificado.mph');

% Best geometry saved by extrair_melhor_dimensao.m (path from your run):
bestMatPath = 'C:\Users\NanoPhotonicsGroup\.comsol\v64\llmatlab\melhor_dimensao.mat';
if ~isfile(bestMatPath)
    % Fallbacks: current folder, then project root.
    if isfile(fullfile(pwd,'melhor_dimensao.mat'))
        bestMatPath = fullfile(pwd,'melhor_dimensao.mat');
    elseif isfile(fullfile(projectRootDir,'melhor_dimensao.mat'))
        bestMatPath = fullfile(projectRootDir,'melhor_dimensao.mat');
    end
end

% Refractive indices sampled to measure the sensitivity slope.
% Must include the reference index (1.33). More points -> nicer fit.
refractiveIndexList              = 1.30:0.01:1.36;
trackingReferenceRefractiveIndex = 1.33;   % resonance is locked here
trackingHalfWindowDeg            = 20;     % follow the peak within +/-20 deg

% Alpha sweep used to locate the peaks
alphaStartDeg = 0;
alphaStepDeg  = 0.01;
alphaStopDeg  = 89;
globalAlphaMinDeg = alphaStartDeg;
globalAlphaMaxDeg = alphaStopDeg;

% Sweep every n on the full alpha range (true) or only the tracked window
% around the reference peak (false, faster). Peaks/metric are identical.
sweepAllNonReferenceCurves = false;

% Output folder for the figures AND CSVs (folder named "graficos")
outDir = fullfile(projectRootDir,'graficos');
if ~exist(outDir,'dir'), mkdir(outDir); end
figFormats = {'png','pdf'};

%% ------------------------ Model Tags/Params ------------------------
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

%% ------------------------- Load model ------------------------------
ModelUtil.clear;
model = mphload(comsolModelFile);
ModelUtil.showProgress(true);

%% ============================ INITIAL ==============================
initialGeom = struct();
initialGeom.L_domain_nm = model.param.evaluate(PARAM_LDOM) * 1e9;
initialGeom.l_dente_nm  = model.param.evaluate(PARAM_LDEN)  * 1e9;
initialGeom.h_si_nm     = model.param.evaluate(PARAM_HSI)   * 1e9;
initialGeom.h_au_nm     = model.param.evaluate(PARAM_HAU)   * 1e9;

fprintf('INITIAL geometry (from .mph):\n');
fprintf('  L_domain = %.4g nm | l_dente = %.4g nm | h_si = %.4g nm | h_au = %.4g nm\n\n', ...
    initialGeom.L_domain_nm, initialGeom.l_dente_nm, initialGeom.h_si_nm, initialGeom.h_au_nm);

initEval = evaluateTrackedSensitivityAndCurves( ...
    model, STUDY_TAG, PARAM_N, refractiveIndexList, ...
    trackingReferenceRefractiveIndex, trackingHalfWindowDeg, ...
    ALPHA_NAME, MSIGN_NAME, alphaStartDeg, alphaStepDeg, alphaStopDeg, ...
    globalAlphaMinDeg, globalAlphaMaxDeg, sweepAllNonReferenceCurves, ...
    M_PLUS, M_MINUS, TPLUS_TABLE_TAG, TMINUS_TABLE_TAG);
S_init = initEval.sensitivitySlope;

%% ============================== BEST ===============================
if ~isfile(bestMatPath)
    error('Best geometry file not found. Adjust bestMatPath. Tried: %s', bestMatPath);
end
loaded = load(bestMatPath,'best');
best = loaded.best;

fprintf('BEST geometry (from %s):\n', bestMatPath);
fprintf('  L_domain = %.4g nm | l_dente = %.4g nm | h_si = %.4g nm | h_au = %.4g nm\n\n', ...
    best.L_domain_nm, best.l_dente_nm, best.h_si_nm, best.h_au_nm);

setParamNm(model, PARAM_LDOM, best.L_domain_nm);
setParamNm(model, PARAM_LDEN, best.l_dente_nm);
setParamNm(model, PARAM_HSI,  best.h_si_nm);
setParamNm(model, PARAM_HAU,  best.h_au_nm);

bestEval = evaluateTrackedSensitivityAndCurves( ...
    model, STUDY_TAG, PARAM_N, refractiveIndexList, ...
    trackingReferenceRefractiveIndex, trackingHalfWindowDeg, ...
    ALPHA_NAME, MSIGN_NAME, alphaStartDeg, alphaStepDeg, alphaStopDeg, ...
    globalAlphaMinDeg, globalAlphaMaxDeg, sweepAllNonReferenceCurves, ...
    M_PLUS, M_MINUS, TPLUS_TABLE_TAG, TMINUS_TABLE_TAG);
S_best = bestEval.sensitivitySlope;

%% ============================= PLOTS ===============================
nList = refractiveIndexList(:);

% 1) Initial geometry
figInit = makeSensitivityFigure(nList, initEval.alphaPeakDegreesByN(:), initEval.linearFit, ...
    S_init, 'Angular sensitivity - initial geometry', initialGeom);
saveFigure(figInit, outDir, 'sensitivity_vs_n_initial', figFormats);

% 2) Best geometry
figBest = makeSensitivityFigure(nList, bestEval.alphaPeakDegreesByN(:), bestEval.linearFit, ...
    S_best, 'Angular sensitivity - best parameters', best);
saveFigure(figBest, outDir, 'sensitivity_vs_n_best', figFormats);

% 3) Overlay comparison
figCmp = makeOverlayFigure(nList, ...
    initEval.alphaPeakDegreesByN(:), initEval.linearFit, S_init, ...
    bestEval.alphaPeakDegreesByN(:), bestEval.linearFit, S_best);
saveFigure(figCmp, outDir, 'sensitivity_vs_n_overlay', figFormats);

%% ============================= DATA ================================
writetable(table(nList, initEval.alphaPeakDegreesByN(:), initEval.trackedTmokeAbsByN(:), ...
    'VariableNames', {'n','alpha_peak_deg','absTMOKE_at_peak'}), ...
    fullfile(outDir,'sensitivity_initial.csv'));
writetable(table(nList, bestEval.alphaPeakDegreesByN(:), bestEval.trackedTmokeAbsByN(:), ...
    'VariableNames', {'n','alpha_peak_deg','absTMOKE_at_peak'}), ...
    fullfile(outDir,'sensitivity_best.csv'));

fprintf('\nDONE.\n');
fprintf('  INITIAL: S = %+.5f deg/RIU\n', S_init);
fprintf('  BEST   : S = %+.5f deg/RIU\n', S_best);
fprintf('  Figures/CSVs saved in:\n  %s\n', outDir);

%% ======================== local functions =========================
function fig = makeSensitivityFigure(nList, alphaPeakByN, linFit, S, ttl, geom)
    fig = figure('Color','w','Name',ttl,'NumberTitle','off');
    ax = axes(fig); hold(ax,'on'); grid(ax,'on');
    nFine = linspace(min(nList), max(nList), 200).';
    plot(ax, nFine, polyval(linFit, nFine), '-', 'LineWidth',1.4, ...
        'DisplayName', sprintf('linear fit: S = %+.4f deg/RIU', S));
    plot(ax, nList, alphaPeakByN, 'o', 'MarkerSize',8, 'LineWidth',1.2, ...
        'DisplayName','\alpha_{peak}(n)');
    xlabel(ax,'n (refractive index)'); ylabel(ax,'\alpha_{peak} [deg]');
    title(ax, ttl);
    legend(ax,'Location','best');

    % Annotation box: only the relevant geometry data + the slope.
    annText = {
        sprintf('L\\_domain = %.4g nm', geom.L_domain_nm)
        sprintf('l\\_dente = %.4g nm',  geom.l_dente_nm)
        sprintf('h\\_si = %.4g nm',     geom.h_si_nm)
        sprintf('h\\_au = %.4g nm',     geom.h_au_nm)
        sprintf('S = %+.4f deg/RIU',    S)
    };
    annotation(fig,'textbox',[0.15 0.66 0.30 0.24], 'String',annText, ...
        'FitBoxToText','on','BackgroundColor','w','EdgeColor',[0.6 0.6 0.6], ...
        'Interpreter','tex','FontSize',9);
end

function fig = makeOverlayFigure(nList, aInit, fitInit, S_init, aBest, fitBest, S_best)
    fig = figure('Color','w','Name','Angular sensitivity - initial vs best','NumberTitle','off');
    ax = axes(fig); hold(ax,'on'); grid(ax,'on');
    nFine = linspace(min(nList), max(nList), 200).';
    cInit = [0.00 0.45 0.74];   cBest = [0.85 0.33 0.10];

    plot(ax, nFine, polyval(fitInit,nFine), '-',  'Color',cInit, 'LineWidth',1.4, ...
        'DisplayName', sprintf('initial fit: S = %+.4f deg/RIU', S_init));
    plot(ax, nList, aInit, 'o', 'Color',cInit, 'MarkerFaceColor',cInit, 'MarkerSize',7, ...
        'HandleVisibility','off');
    plot(ax, nFine, polyval(fitBest,nFine), '-',  'Color',cBest, 'LineWidth',1.4, ...
        'DisplayName', sprintf('best fit: S = %+.4f deg/RIU', S_best));
    plot(ax, nList, aBest, 's', 'Color',cBest, 'MarkerFaceColor',cBest, 'MarkerSize',7, ...
        'HandleVisibility','off');

    xlabel(ax,'n (refractive index)'); ylabel(ax,'\alpha_{peak} [deg]');
    title(ax,'Angular sensitivity: initial vs best geometry');
    legend(ax,'Location','best');
end

function saveFigure(fig, outDir, baseName, formats)
    for k = 1:numel(formats)
        ext = lower(string(formats{k}));
        fp = fullfile(outDir, string(baseName) + "." + ext);
        try
            switch ext
                case "png", exportgraphics(fig, fp, 'Resolution',300);
                case "pdf", exportgraphics(fig, fp, 'ContentType','vector');
                otherwise,  saveas(fig, fp);
            end
        catch ME
            warning('Failed to save %s: %s', fp, ME.message);
        end
    end
end

% ------- helpers copied verbatim from senseAndTmoke_semCeyig.m -------
function setParamNm(mdl, name, val_nm)
    mdl.param.set(name, sprintf('%.12g[nm]', val_nm));
end

function setParamScalar(mdl, name, val)
    mdl.param.set(name, sprintf('%.12g', val));
end

function trackedEval = evaluateTrackedSensitivityAndCurves( ...
        mdl, studyTag, PARAM_N, refractiveIndexList, ...
        trackingReferenceRefractiveIndex, trackingHalfWindowDeg, ...
        alphaName, mName, referenceSweepStartDeg, referenceSweepStepDeg, referenceSweepStopDeg, ...
        globalAlphaMinDeg, globalAlphaMaxDeg, sweepAllNonReferenceCurves, ...
        mPlusStr, mMinusStr, ttagPlus, ttagMinus)
% Follow the same TMOKE resonance across n: lock it at the reference index,
% then only look inside alpha_ref +/- trackingHalfWindowDeg for the others.
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

    % m = +1
    redirectAllNumericalsToTable(mdl, ttagPlus);
    clearTable(mdl, ttagPlus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mPlusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha1_deg, T1] = readAlphaAndTFromNamedTable(mdl, ttagPlus, Npts);

    % m = -1
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
