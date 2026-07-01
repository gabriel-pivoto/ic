%% =====================================================================
% TMOKE x alpha: INITIAL geometry vs BEST geometry
% ---------------------------------------------------------------------
% What it does (no search, just two curves):
%   1) Loads the COMSOL model exactly as saved (.mph).
%   2) Computes TMOKE(alpha) for the model's INITIAL geometry.
%   3) Sets the BEST geometry (from melhor_dimensao.mat, produced by
%      extrair_melhor_dimensao.m) and computes TMOKE(alpha) for it.
%   4) Plots each curve in its own figure with only the relevant data
%      (geometry parameters + peak), and also an overlay comparison.
%
% Both sweeps use alpha = 0 : 0.01 : 89 deg at n = baselineRefractiveIndex.
%
% TMOKE definition (same as senseAndTmoke_semCeyig.m):
%   TMOKE = 2*(T+ - T-)/(T+ + T-)
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

% Alpha sweep (dense)
alphaStartDeg = 0;
alphaStepDeg  = 0.01;
alphaStopDeg  = 89;

% Refractive index used for BOTH curves (same baseline used by the search).
baselineRefractiveIndex = 1.33;

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
% Read the geometry stored in the .mph (do NOT change it), then solve.
initialGeom = struct();
initialGeom.L_domain_nm = model.param.evaluate(PARAM_LDOM) * 1e9;
initialGeom.l_dente_nm  = model.param.evaluate(PARAM_LDEN)  * 1e9;
initialGeom.h_si_nm     = model.param.evaluate(PARAM_HSI)   * 1e9;
initialGeom.h_au_nm     = model.param.evaluate(PARAM_HAU)   * 1e9;

fprintf('INITIAL geometry (from .mph):\n');
fprintf('  L_domain = %.4g nm | l_dente = %.4g nm | h_si = %.4g nm | h_au = %.4g nm\n\n', ...
    initialGeom.L_domain_nm, initialGeom.l_dente_nm, initialGeom.h_si_nm, initialGeom.h_au_nm);

setParamScalar(model, PARAM_N, baselineRefractiveIndex);
[alphaInit, TplusInit, TminusInit, tmokeInit] = solveAndGetTplusTminus( ...
    model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
    alphaStartDeg, alphaStepDeg, alphaStopDeg, ...
    M_PLUS, M_MINUS, TPLUS_TABLE_TAG, TMINUS_TABLE_TAG);
[tmokePeakInit, alphaPeakInit] = peakInfo(alphaInit, tmokeInit);

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
setParamScalar(model, PARAM_N, baselineRefractiveIndex);

[alphaBest, TplusBest, TminusBest, tmokeBest] = solveAndGetTplusTminus( ...
    model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
    alphaStartDeg, alphaStepDeg, alphaStopDeg, ...
    M_PLUS, M_MINUS, TPLUS_TABLE_TAG, TMINUS_TABLE_TAG);
[tmokePeakBest, alphaPeakBest] = peakInfo(alphaBest, tmokeBest);

%% ============================= PLOTS ===============================
% 1) Initial geometry
figInit = makeTmokeFigure(alphaInit, tmokeInit, alphaPeakInit, tmokePeakInit, ...
    'TMOKE x \alpha - geometria inicial', initialGeom, baselineRefractiveIndex);
saveFigure(figInit, outDir, 'tmoke_alpha_inicial', figFormats);

% 2) Best geometry
figBest = makeTmokeFigure(alphaBest, tmokeBest, alphaPeakBest, tmokePeakBest, ...
    'TMOKE x \alpha - melhores parametros', best, baselineRefractiveIndex);
saveFigure(figBest, outDir, 'tmoke_alpha_melhor', figFormats);

%% ============================= DATA ================================
% Save the raw curves so a later script can re-plot without COMSOL.
writetable(table(alphaInit(:), TplusInit(:), TminusInit(:), tmokeInit(:), ...
    'VariableNames', {'alpha_deg','Tplus','Tminus','TMOKE'}), ...
    fullfile(outDir,'curva_inicial.csv'));
writetable(table(alphaBest(:), TplusBest(:), TminusBest(:), tmokeBest(:), ...
    'VariableNames', {'alpha_deg','Tplus','Tminus','TMOKE'}), ...
    fullfile(outDir,'curva_melhor.csv'));

fprintf('\nDONE.\n');
fprintf('  INITIAL: |TMOKE|max = %.5f @ alpha = %.3f deg\n', abs(tmokePeakInit), alphaPeakInit);
fprintf('  BEST   : |TMOKE|max = %.5f @ alpha = %.3f deg\n', abs(tmokePeakBest), alphaPeakBest);
fprintf('  Figures/CSVs saved in:\n  %s\n', outDir);

%% ======================== local functions =========================
function [tmokePeak, alphaPeak] = peakInfo(alphaDeg, tmokeCurve)
    [~, k] = max(abs(tmokeCurve));
    tmokePeak = tmokeCurve(k);
    alphaPeak = alphaDeg(k);
end

function fig = makeTmokeFigure(alphaDeg, tmokeCurve, alphaPeak, tmokePeak, ttl, geom, nVal)
    fig = figure('Color','w','Name',ttl,'NumberTitle','off');
    ax = axes(fig); hold(ax,'on'); grid(ax,'on');
    plot(ax, alphaDeg, tmokeCurve, '-', 'LineWidth',1.4, 'DisplayName','TMOKE(\alpha)');
    plot(ax, alphaPeak, tmokePeak, 'o', 'MarkerSize',8, 'LineWidth',1.2, ...
        'DisplayName',sprintf('pico: %.4f @ %.2f\\circ', tmokePeak, alphaPeak));
    xlabel(ax,'\alpha [deg]'); ylabel(ax,'TMOKE');
    xlim(ax,[alphaDeg(1) alphaDeg(end)]);
    title(ax, ttl);
    legend(ax,'Location','best');

    % Annotation box: only the relevant geometry data.
    annText = {
        sprintf('L\\_domain = %.4g nm', geom.L_domain_nm)
        sprintf('l\\_dente = %.4g nm',  geom.l_dente_nm)
        sprintf('h\\_si = %.4g nm',     geom.h_si_nm)
        sprintf('h\\_au = %.4g nm',     geom.h_au_nm)
        sprintf('n = %.3f',             nVal)
        sprintf('|TMOKE|_{max} = %.4f', abs(tmokePeak))
        sprintf('\\alpha_{pico} = %.2f\\circ', alphaPeak)
    };
    annotation(fig,'textbox',[0.15 0.62 0.28 0.28], 'String',annText, ...
        'FitBoxToText','on','BackgroundColor','w','EdgeColor',[0.6 0.6 0.6], ...
        'Interpreter','tex','FontSize',9);
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
