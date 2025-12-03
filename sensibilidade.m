% =====================================================================
% TMOKE 5D — Otimização de Sensibilidade vs n (deg/RIU)
%
% - Varre as 5 camadas: (h_au, h_ceyig, L_domain, l_dente, h_si)
% - Para cada combinação:
%       * para n em [1.33 1.36 1.39]:
%             - seta n no modelo
%             - roda solveAndGetRplusRminus(alpha sweep)
%             - acha alpha_peak(n) onde |TMOKE| é máximo
%       * faz regressão linear alpha_peak(n) => S = d(alpha)/dn [deg/RIU]
% - Objetivo: maximizar |S|
%
% Estágios:
%   1) COARSE: grade grossa 5D
%   2) FINE: janela em torno do melhor ponto (ou TOP-K)
%
% Checkpoint:
%   - Arquivo: /checkpoints/tmoke_sens_5d_checkpoint.mat
%   - Campos: stage, done_points, runs_done_global, Tcoarse, Tfine, best
% =====================================================================

clear; clc; close all; format long;
import com.comsol.model.*
import com.comsol.model.util.*

% --------------------------- Paths/Model ------------------------------
homeDir  = 'C:\Users\gabri\Documents\projetoIC';   % ajuste se precisar
mph_file = fullfile(homeDir,'usandoMatlab.mph');
addpath(genpath(homeDir));

CKPT_DIR  = fullfile(homeDir,'checkpoints');
if ~exist(CKPT_DIR,'dir'), mkdir(CKPT_DIR); end
CKPT_FILE = fullfile(CKPT_DIR,'tmoke_sens_5d_checkpoint.mat');

% ------------------------- Tags/params COMSOL ------------------------
STUDY_TAG   = 'std1';      % Study com sweep em alpha e m
PARAM_HAU   = 'h_au';
PARAM_HCEY  = 'h_ceyig';
PARAM_LDOM  = 'L_domain';
PARAM_LDEN  = 'l_dente';
PARAM_HSI   = 'h_si';
PARAM_N     = 'n';         % índice do meio externo
ALPHA_NAME  = 'alpha';     % ângulo [deg]
MSIGN_NAME  = 'm';         % magnetização
M_PLUS      = '1';
M_MINUS     = '-1';

RPLUS_TABLE_TAG  = 'tblRplus';
RMINUS_TABLE_TAG = 'tblRminus';

% ------------------------ Sensibilidade vs n --------------------------
na = [1.33 1.36 1.39];         % índices para sensibilidade
alpha_start = 0;               % varredura em alpha (COARSE/FINE)
alpha_step_coarse = 1.0;
alpha_step_fine   = 0.1;
alpha_stop  = 89;

% ------------------------ Grade COARSE das camadas -------------------
% Loop order: L_domain -> l_dente -> h_si -> h_ceyig -> h_au
Ldomain_grid_nm = 800:50:850;         % 3
ldente_grid_nm  = 500:50:600;         % 3
hsi_grid_nm     = [220 240 260];      % 3
hcey_grid_nm    = [100 140];          % 2
hau_grid_nm     = 20:10:60;           % 5
% Total pontos COARSE = 3*3*3*2*5 = 270

% ------------------------ Janela FINE em torno do melhor -------------
fine_hau_delta   = 2;    fine_hau_step   = 1;
fine_hcey_delta  = 5;    fine_hcey_step  = 1;
fine_Ldom_delta  = 10;   fine_Ldom_step  = 5;
fine_Lden_delta  = 10;   fine_Lden_step  = 5;
fine_hsi_delta   = 5;    fine_hsi_step   = 5;

TOPK_COARSE = 1;   % quantos seeds do COARSE entram no FINE

% ---------------------------- Checkpoint -----------------------------
RESUME      = false;
RES_STAGE   = "NONE";
Tcoarse     = [];
Tfine       = [];
bestStruct  = struct();

if isfile(CKPT_FILE)
    try
        Sckpt = load(CKPT_FILE,'CKPT');
        CKPT = Sckpt.CKPT;
        RESUME    = true;
        RES_STAGE = string(CKPT.stage);
        if isfield(CKPT,'Tcoarse'),    Tcoarse    = CKPT.Tcoarse;    end
        if isfield(CKPT,'Tfine'),      Tfine      = CKPT.Tfine;      end
        if isfield(CKPT,'bestStruct'), bestStruct = CKPT.bestStruct; end
        fprintf('>>> RESUME detectado: stage=%s | done_points=%d | runs=%d\n', ...
            RES_STAGE, CKPT.done_points, CKPT.runs_done_global);
    catch
        warning('Checkpoint existe mas não pôde ser carregado. Iniciando do zero.');
        RESUME = false;
        RES_STAGE = "NONE";
    end
end

runs_done_global   = 0;
done_points_global = 0;
CKPT_EVERY_POINTS  = 10;  % salva a cada 10 pontos de grade

% --------------------------- Carrega modelo --------------------------
ModelUtil.clear;
model = mphload(mph_file);
ModelUtil.showProgress(true);

t0 = tic;

% =====================================================================
%                             COARSE
% =====================================================================
coarse_done_points = 0;
if RESUME && any(RES_STAGE == ["FINE","VALID"])
    fprintf('Pulando COARSE (já feito).\n');
else
    fprintf('\n===== STAGE COARSE =====\n');
    Rows = [];   % Ldom, Lden, hsi, hcey, hau, S, S_sign

    if RESUME && RES_STAGE == "COARSE"
        runs_done_global   = CKPT.runs_done_global;
        coarse_done_points = CKPT.done_points;
        if isfield(CKPT,'Rows'), Rows = CKPT.Rows; end
        fprintf('>>> RESUME COARSE: já tinha %d pontos.\n', coarse_done_points);
    end

    nL = numel(Ldomain_grid_nm);
    nD = numel(ldente_grid_nm);
    nS = numel(hsi_grid_nm);
    nG = numel(hcey_grid_nm);
    nH = numel(hau_grid_nm);

    flat_idx = 0;

    for iL = 1:nL
        setParamNm(model, PARAM_LDOM, Ldomain_grid_nm(iL));
        for iD = 1:nD
            setParamNm(model, PARAM_LDEN, ldente_grid_nm(iD));
            for iS = 1:nS
                setParamNm(model, PARAM_HSI, hsi_grid_nm(iS));
                for ig = 1:nG
                    setParamNm(model, PARAM_HCEY, hcey_grid_nm(ig));
                    for ih = 1:nH
                        setParamNm(model, PARAM_HAU, hau_grid_nm(ih));

                        flat_idx = flat_idx + 1;
                        if flat_idx <= coarse_done_points
                            continue;   % já foi calculado antes do checkpoint
                        end

                        % ---- avalia sensibilidade S para este ponto 5D ----
                        [Sval, alpha_peaks] = evaluate_sensitivity( ...
                            model, STUDY_TAG, PARAM_N, na, ...
                            ALPHA_NAME, MSIGN_NAME, ...
                            alpha_start, alpha_step_coarse, alpha_stop, ...
                            M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                        Rows = [Rows; ...
                            Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), ...
                            hcey_grid_nm(ig), hau_grid_nm(ih), Sval, sign(Sval)]; %#ok<AGROW>

                        runs_done_global   = runs_done_global + numel(na)*2;  % cada n -> 2 runs (m=±1)
                        coarse_done_points = coarse_done_points + 1;
                        done_points_global = done_points_global + 1;

                        fprintf(['COARSE | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                 ' | S=%.6f deg/RIU (sign=%+d) | ponto=%d\n'], ...
                            Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), ...
                            hcey_grid_nm(ig), hau_grid_nm(ih), Sval, sign(Sval), coarse_done_points);

                        % ---- checkpoint simples ----
                        if mod(coarse_done_points, CKPT_EVERY_POINTS) == 0
                            CKPT.stage            = 'COARSE';
                            CKPT.done_points      = coarse_done_points;
                            CKPT.runs_done_global = runs_done_global;
                            CKPT.Rows             = Rows;
                            CKPT.Tcoarse          = [];
                            CKPT.Tfine            = Tfine;
                            CKPT.bestStruct       = bestStruct;
                            save(CKPT_FILE,'CKPT','-v7.3');
                            fprintf('Checkpoint COARSE salvo (%d pontos).\n', coarse_done_points);
                        end
                    end
                end
            end
        end
    end

    Tcoarse = array2table(Rows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm','S_deg_per_RIU','S_sign'});

    % guarda Tcoarse no checkpoint final do estágio
    CKPT.stage            = 'COARSE';
    CKPT.done_points      = coarse_done_points;
    CKPT.runs_done_global = runs_done_global;
    CKPT.Tcoarse          = Tcoarse;
    CKPT.Tfine            = Tfine;
    CKPT.bestStruct       = bestStruct;
    save(CKPT_FILE,'CKPT','-v7.3');
    fprintf('COARSE finalizado: %d pontos.\n', coarse_done_points);
end

% =====================================================================
%                  Seleção de seeds (TOP-K) para FINE
% =====================================================================
if isempty(Tcoarse)
    error('Tcoarse vazio (nada calculado em COARSE).');
end

[~, ordC] = sort(abs(Tcoarse.S_deg_per_RIU), 'descend');
kkeepC    = min(TOPK_COARSE, numel(ordC));
seeds     = Tcoarse(ordC(1:kkeepC), :);

fprintf('\nSeeds para FINE (TOP-%d por |S|):\n', kkeepC);
disp(seeds);

% =====================================================================
%                             FINE
% =====================================================================
fine_done_points = 0;
if RESUME && any(RES_STAGE == ["VALID"])
    fprintf('Pulando FINE (já feito).\n');
else
    fprintf('\n===== STAGE FINE =====\n');
    FineRows = [];   % Ldom, Lden, h_si, h_cey, h_au, S, S_sign

    if RESUME && RES_STAGE == "FINE"
        runs_done_global   = CKPT.runs_done_global;
        fine_done_points   = CKPT.done_points;
        if isfield(CKPT,'FineRows'), FineRows = CKPT.FineRows; end
        if isfield(CKPT,'seeds'),    seeds    = CKPT.seeds;    end
        fprintf('>>> RESUME FINE: já tinha %d pontos.\n', fine_done_points);
    end

    flat_idx = 0;

    for s = 1:height(seeds)
        % listas ao redor do seed (clampadas ao grid COARSE)
        Hlist = max(min(hau_grid_nm),  seeds.h_au_nm(s)   - fine_hau_delta)  : fine_hau_step  : min(max(hau_grid_nm),  seeds.h_au_nm(s)   + fine_hau_delta);
        Glist = max(min(hcey_grid_nm), seeds.h_ceyig_nm(s)- fine_hcey_delta) : fine_hcey_step : min(max(hcey_grid_nm), seeds.h_ceyig_nm(s)+ fine_hcey_delta);
        Llist = max(min(Ldomain_grid_nm), seeds.L_domain_nm(s)- fine_Ldom_delta) : fine_Ldom_step : min(max(Ldomain_grid_nm), seeds.L_domain_nm(s)+ fine_Ldom_delta);
        Dlist = max(min(ldente_grid_nm),  seeds.l_dente_nm(s) - fine_Lden_delta) : fine_Lden_step : min(max(ldente_grid_nm),  seeds.l_dente_nm(s) + fine_Lden_delta);
        Slist = max(min(hsi_grid_nm),     seeds.h_si_nm(s)    - fine_hsi_delta)  : fine_hsi_step  : min(max(hsi_grid_nm),     seeds.h_si_nm(s)    + fine_hsi_delta);

        for iL = 1:numel(Llist)
            setParamNm(model, PARAM_LDOM, Llist(iL));
            for iD = 1:numel(Dlist)
                setParamNm(model, PARAM_LDEN, Dlist(iD));
                for iS = 1:numel(Slist)
                    setParamNm(model, PARAM_HSI, Slist(iS));
                    for ig = 1:numel(Glist)
                        setParamNm(model, PARAM_HCEY, Glist(ig));
                        for ih = 1:numel(Hlist)
                            setParamNm(model, PARAM_HAU, Hlist(ih));

                            flat_idx = flat_idx + 1;
                            if flat_idx <= fine_done_points
                                continue;
                            end

                            % FINE: alpha mais fino
                            [Sval, alpha_peaks] = evaluate_sensitivity( ...
                                model, STUDY_TAG, PARAM_N, na, ...
                                ALPHA_NAME, MSIGN_NAME, ...
                                alpha_start, alpha_step_fine, alpha_stop, ...
                                M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                            FineRows = [FineRows; ...
                                Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), Sval, sign(Sval)]; %#ok<AGROW>

                            runs_done_global = runs_done_global + numel(na)*2;
                            fine_done_points = fine_done_points + 1;
                            done_points_global = done_points_global + 1;

                            fprintf(['FINE   | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                     ' | S=%.6f deg/RIU (sign=%+d) | ponto=%d\n'], ...
                                Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), Sval, sign(Sval), fine_done_points);

                            % checkpoint
                            if mod(fine_done_points, CKPT_EVERY_POINTS) == 0
                                CKPT.stage            = 'FINE';
                                CKPT.done_points      = fine_done_points;
                                CKPT.runs_done_global = runs_done_global;
                                CKPT.FineRows         = FineRows;
                                CKPT.seeds            = seeds;
                                CKPT.Tcoarse          = Tcoarse;
                                CKPT.Tfine            = [];
                                CKPT.bestStruct       = bestStruct;
                                save(CKPT_FILE,'CKPT','-v7.3');
                                fprintf('Checkpoint FINE salvo (%d pontos).\n', fine_done_points);
                            end
                        end
                    end
                end
            end
        end
    end

    Tfine = array2table(FineRows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm','S_deg_per_RIU','S_sign'});

    CKPT.stage            = 'FINE';
    CKPT.done_points      = fine_done_points;
    CKPT.runs_done_global = runs_done_global;
    CKPT.Tcoarse          = Tcoarse;
    CKPT.Tfine            = Tfine;
    CKPT.bestStruct       = bestStruct;
    save(CKPT_FILE,'CKPT','-v7.3');
    fprintf('FINE finalizado: %d pontos.\n', fine_done_points);
end

% =====================================================================
%                       Escolhe melhor (global) por |S|
% =====================================================================
AllCand = Tcoarse;
if ~isempty(Tfine)
    AllCand = [AllCand; Tfine]; %#ok<AGROW>
end

[~, ordAll] = sort(abs(AllCand.S_deg_per_RIU), 'descend');
bestRow     = AllCand(ordAll(1), :);

bestStruct.Ldom_best  = bestRow.L_domain_nm;
bestStruct.Lden_best  = bestRow.l_dente_nm;
bestStruct.hsi_best   = bestRow.h_si_nm;
bestStruct.hcey_best  = bestRow.h_ceyig_nm;
bestStruct.hau_best   = bestRow.h_au_nm;
bestStruct.S_best     = bestRow.S_deg_per_RIU;

fprintf('\n===== MELHOR CONFIGURAÇÃO POR |S| =====\n');
disp(bestRow);

CKPT.stage            = 'VALID';
CKPT.done_points      = 0;
CKPT.runs_done_global = runs_done_global;
CKPT.Tcoarse          = Tcoarse;
CKPT.Tfine            = Tfine;
CKPT.bestStruct       = bestStruct;
save(CKPT_FILE,'CKPT','-v7.3');

% =====================================================================
%                 VALID: curva densa + plot de sensibilidade
% =====================================================================
fprintf('\n===== STAGE VALID (reavalia melhor ponto) =====\n');

setParamNm(model, PARAM_LDOM, bestStruct.Ldom_best);
setParamNm(model, PARAM_LDEN, bestStruct.Lden_best);
setParamNm(model, PARAM_HSI,  bestStruct.hsi_best);
setParamNm(model, PARAM_HCEY, bestStruct.hcey_best);
setParamNm(model, PARAM_HAU,  bestStruct.hau_best);

alpha_step_dense = 0.01;

S_dense      = 0;
alpha_peaks  = zeros(size(na));
TM_all       = cell(size(na));
alpha_all    = cell(size(na));

for i = 1:numel(na)
    n_val = na(i);
    model.param.set(PARAM_N, sprintf('%.12g', n_val));

    [alpha_deg, Rplus, Rminus, TM] = solveAndGetRplusRminus( ...
        model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
        alpha_start, alpha_step_dense, alpha_stop, ...
        M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

    [~, k] = max(abs(TM));
    alpha_peaks(i) = alpha_deg(k);
    TM_all{i}      = TM;
    alpha_all{i}   = alpha_deg;

    runs_done_global = runs_done_global + 2;
end

p = polyfit(na, alpha_peaks, 1);
S_dense = p(1);

fprintf('alpha_peak (deg) = [');
fprintf(' %.4f', alpha_peaks);
fprintf(' ]\n');
fprintf('Sensibilidade (VALID) S ≈ %.6f deg/RIU\n', S_dense);

% ------------------------- Plots finais ------------------------------
figure('Name','Sensibilidade TMOKE vs n','NumberTitle','off');
hold on; grid on;
plot(na, alpha_peaks, 'o','LineWidth',1.5,'DisplayName','dados (\alpha_{peak})');
na_fit    = linspace(min(na), max(na), 100);
alpha_fit = polyval(p, na_fit);
plot(na_fit, alpha_fit, '-','LineWidth',1.2, ...
    'DisplayName', sprintf('ajuste linear (S = %.3f deg/RIU)', S_dense));
xlabel('Índice de refração n');
ylabel('\alpha_{peak} [deg]');
title('Sensibilidade de \alpha_{peak} em função de n (melhor 5D)');
legend('Location','best');

figure('Name','TMOKE vs \alpha para diferentes n (best 5D)','NumberTitle','off');
hold on; grid on;
for i = 1:numel(na)
    plot(alpha_all{i}, TM_all{i}, 'LineWidth',1.2, ...
        'DisplayName', sprintf('n = %.2f', na(i)));
end
xlabel('\alpha [deg]');
ylabel('TMOKE(\alpha)');
title('TMOKE(\alpha) para n = 1.33, 1.36, 1.39 (melhor 5D)');
legend('Location','best');

% --------------------------- CSVs ------------------------------------
writetable(Tcoarse, fullfile(homeDir,'tmoke_sens_5D_coarse.csv'));
if ~isempty(Tfine)
    writetable(Tfine, fullfile(homeDir,'tmoke_sens_5D_fine.csv'));
end

BestSensRow = bestRow;
BestSensRow.alpha_peaks_deg = {alpha_peaks};
BestSensRow.S_dense_deg_per_RIU = S_dense;
writetable(BestSensRow, fullfile(homeDir,'tmoke_sens_5D_best.csv'));

% -------------------------- Summary ----------------------------------
elapsed_total = toc(t0);
fprintf('\n===== SUMMARY =====\n');
fprintf('Best 5D:\n');
fprintf('  L_domain = %4.0f nm\n', bestStruct.Ldom_best);
fprintf('  l_dente  = %4.0f nm\n', bestStruct.Lden_best);
fprintf('  h_si     = %3.0f nm\n', bestStruct.hsi_best);
fprintf('  h_ceyig  = %.3f nm\n', bestStruct.hcey_best);
fprintf('  h_au     = %.3f nm\n', bestStruct.hau_best);
fprintf('  S (FINE) ≈ %.6f deg/RIU\n', bestStruct.S_best);
fprintf('  S (VALID)≈ %.6f deg/RIU\n', S_dense);
fprintf('Runs totais ≈ %d\n', runs_done_global);
fprintf('Tempo total: %.1f s\n', elapsed_total);

% =====================================================================
%                           FUNÇÕES LOCAIS
% =====================================================================
function setParamNm(mdl, name, val_nm)
    mdl.param.set(name, sprintf('%.12g[nm]', val_nm));
end

function [Sval, alpha_peaks] = evaluate_sensitivity( ...
        mdl, studyTag, PARAM_N, na, ...
        alphaName, mName, aStartDeg, aStepDeg, aStopDeg, ...
        mPlusStr, mMinusStr, ttagPlus, ttagMinus)
% Calcula S = d(alpha_peak)/dn (deg/RIU) para uma geometria fixa.
    alpha_peaks = zeros(size(na));
    for i = 1:numel(na)
        n_val = na(i);
        mdl.param.set(PARAM_N, sprintf('%.12g', n_val));

        [alpha_deg, Rplus, Rminus, TM] = solveAndGetRplusRminus( ...
            mdl, studyTag, alphaName, mName, ...
            aStartDeg, aStepDeg, aStopDeg, ...
            mPlusStr, mMinusStr, ttagPlus, ttagMinus);

        if isempty(TM)
            error('TMOKE vazio na avaliação de sensibilidade.');
        end

        [~, k] = max(abs(TM));
        alpha_peaks(i) = alpha_deg(k);
    end

    p = polyfit(na, alpha_peaks, 1);
    Sval = p(1);   % slope [deg/RIU]
end

function [alpha_deg, Rplus, Rminus, TMOKE] = solveAndGetRplusRminus( ...
    mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, ...
    mPlusStr, mMinusStr, ttagPlus, ttagMinus)

    ensureTable(mdl, ttagPlus);
    ensureTable(mdl, ttagMinus);

    Npts = 1 + floor((aStopDeg - aStartDeg)/aStepDeg + 1e-9);

    % m = +1  -> tabela ttagPlus
    redirectAllNumericalsToTable(mdl, ttagPlus);
    clearTable(mdl, ttagPlus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, ...
                     aStartDeg, aStepDeg, aStopDeg, mPlusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha1_deg, R1] = readAlphaAndRFromNamedTable(mdl, ttagPlus, Npts);

    % m = -1  -> tabela ttagMinus
    redirectAllNumericalsToTable(mdl, ttagMinus);
    clearTable(mdl, ttagMinus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, ...
                     aStartDeg, aStepDeg, aStopDeg, mMinusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha2_deg, R2] = readAlphaAndRFromNamedTable(mdl, ttagMinus, Npts);

    assert(numel(alpha1_deg)==Npts && numel(alpha2_deg)==Npts, ...
        'Unexpected α sweep length.');
    assert(max(abs(alpha1_deg - alpha2_deg)) < 1e-9, ...
        'α grids differ between m=+1 and m=-1.');

    alpha_deg = alpha1_deg;
    Rplus  = R1;
    Rminus = R2;

    denom = Rplus + Rminus;
    denom(abs(denom) < 1e-9) = 1e-9;
    TMOKE = (Rplus - Rminus) ./ denom;
end

function setTwoParamSweep(mdl, studyTag, alphaName, mName, ...
                          aStartDeg, aStepDeg, aStopDeg, mStr)
    ptag = getParametricTag(mdl, studyTag);
    mdl.study(studyTag).feature(ptag).set('pname', {alphaName, mName});
    mdl.study(studyTag).feature(ptag).set('punit', {'deg','1'});
    alist = sprintf('range(%.12g[deg], %.12g[deg], %.12g[deg])', ...
        aStartDeg, aStepDeg, aStopDeg);
    mdl.study(studyTag).feature(ptag).set('plistarr', {alist, mStr});
end

function ptag = getParametricTag(mdl, studyTag)
    ptag = 'param';
    try
        fts = cell(mdl.study(studyTag).feature().tags());
        for i = 1:numel(fts)
            typ = char(mdl.study(studyTag).feature(fts{i}).featureType());
            if contains(lower(typ), 'param')
                ptag = fts{i};
                break;
            end
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
    if isempty(heads) && isfield(S,'head')   && ~isempty(S.head),   heads = string(S.head);   end
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
        rIdx = find(contains(hlow,'reflectance') | contains(hlow,'total reflectance') | ...
                    contains(hlow,' total r'), 1, 'first');
        if isempty(aIdx) || isempty(rIdx)
            if ~isfield(S,'data') || isempty(S.data) || size(S.data,2) < 2
                error('Alpha/Reflectance not found in %s and not enough columns.', ttag);
            end
            aIdx = 1; rIdx = 2;
        end
        a = S.data(:, aIdx); 
        R = S.data(:, rIdx);
    end

    assert(numel(a) >= Npts && numel(R) >= Npts, ...
        'Table %s has %d rows; expected >= %d.', ttag, numel(a), Npts);
    a = a(end-Npts+1:end);
    R = R(end-Npts+1:end);

    if ~isempty(heads)
        hlow = lower(heads);
        if any(contains(hlow, '(rad)')) || any(contains(hlow, '[rad]'))
            a = a * 180/pi;
        end
    end

    alpha_deg = a;
    Rcol      = R;
end
