% =====================================================================
% TMOKE 5D — Otimização de Sensibilidade vs n (deg/RIU)
% VERSÃO "TESTE RÁPIDO":
%   - poucos steps de parâmetros (rodar rápido só pra validar o fluxo)
%   - a cada iteração (ponto de geometria) plota:
%       (1) TMOKE(α) para n = [1.33 1.36 1.39]
%       (5) Heatmap TMOKE(α,n)
%   - no final plota TODOS os gráficos úteis (lista completa)
%
% Checkpoint:
%   - /checkpoints/tmoke_sens_5d_checkpoint.mat
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
na = [1.33 1.36 1.39];       % 3 valores de n (pra sensibilidade)
alpha_start        = 0;
alpha_stop         = 89;

% ======= TESTE RÁPIDO (alpha mais grosso) =======
alpha_step_coarse  = 5.0;     % antes: 1.0
alpha_step_fine    = 2.0;     % antes: 0.1
alpha_step_dense   = 0.5;     % antes: 0.01 (VALID)

% ------------------------ Grade COARSE (TESTE RÁPIDO) -----------------
% Pouquíssimos pontos, só pra validar que o pipeline está ok.
% (Depois você aumenta os grids.)
Ldomain_grid_nm = [800 850];     % 2
ldente_grid_nm  = [500 600];     % 2
hsi_grid_nm     = [220 260];     % 2
hcey_grid_nm    = [100 140];     % 2
hau_grid_nm     = [20 40];       % 2
% Total COARSE = 2*2*2*2*2 = 32 pontos

% ------------------------ Janela FINE (TESTE RÁPIDO) ------------------
fine_hau_delta   = 2;    fine_hau_step   = 2;   % poucos pontos
fine_hcey_delta  = 5;    fine_hcey_step  = 5;
fine_Ldom_delta  = 10;   fine_Ldom_step  = 10;
fine_Lden_delta  = 10;   fine_Lden_step  = 10;
fine_hsi_delta   = 10;   fine_hsi_step   = 10;

TOPK_COARSE = 1;   % seeds do COARSE que entram no FINE (deixa 1 no teste)

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
CKPT_EVERY_POINTS  = 5;  % salva a cada 5 pontos no teste

% --------------------------- Carrega modelo --------------------------
ModelUtil.clear;
model = mphload(mph_file);
ModelUtil.showProgress(true);

t0 = tic;

% =====================================================================
%                 ESTIMATIVA DE NÚMERO DE ITERAÇÕES
% =====================================================================
nL = numel(Ldomain_grid_nm);
nD = numel(ldente_grid_nm);
nS = numel(hsi_grid_nm);
nG = numel(hcey_grid_nm);
nH = numel(hau_grid_nm);

coarse_total_pts  = nL*nD*nS*nG*nH;
coarse_total_runs = coarse_total_pts * numel(na) * 2;

fine_H = 1 + floor(2*fine_hau_delta   / fine_hau_step);
fine_G = 1 + floor(2*fine_hcey_delta  / fine_hcey_step);
fine_L = 1 + floor(2*fine_Ldom_delta  / fine_Ldom_step);
fine_D = 1 + floor(2*fine_Lden_delta  / fine_Lden_step);
fine_S = 1 + floor(2*fine_hsi_delta   / fine_hsi_step);
fine_pts_per_seed_est = fine_L*fine_D*fine_S*fine_G*fine_H;
fine_total_pts_est    = TOPK_COARSE * fine_pts_per_seed_est;
fine_total_runs_est   = fine_total_pts_est * numel(na) * 2;

total_pts_est_global  = coarse_total_pts + fine_total_pts_est;
total_runs_est_global = coarse_total_runs + fine_total_runs_est;

fprintf('\n===== ESTIMATIVA DE REPETIÇÕES (TESTE) =====\n');
fprintf('COARSE: %d pontos (~%d runs COMSOL)\n', coarse_total_pts, coarse_total_runs);
fprintf('FINE (estimado): ~%d pontos (~%d runs COMSOL)\n', fine_total_pts_est, fine_total_runs_est);
fprintf('TOTAL (estimado): ~%d pontos (~%d runs COMSOL)\n\n', total_pts_est_global, total_runs_est_global);

iter_global = 0;

% =====================================================================
%                             COARSE
% =====================================================================
coarse_done_points = 0;
if RESUME && any(RES_STAGE == ["FINE","VALID"])
    fprintf('Pulando COARSE (já feito).\n');
else
    fprintf('\n===== STAGE COARSE =====\n');
    Rows = [];   % Ldom, Lden, h_si, h_cey, h_au, S, S_sign

    if RESUME && RES_STAGE == "COARSE"
        runs_done_global   = CKPT.runs_done_global;
        coarse_done_points = CKPT.done_points;
        if isfield(CKPT,'Rows'), Rows = CKPT.Rows; end
        if isfield(CKPT,'iter_global'), iter_global = CKPT.iter_global; end
        fprintf('>>> RESUME COARSE: já tinha %d pontos.\n', coarse_done_points);
    end

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
                            continue;
                        end

                        % ---- avalia S e (IMPORTANTE) retorna TMOKE(α) pra plotar 1 e 5 ----
                        [Sval, alpha_peaks, alpha_grid, TM_mat] = evaluate_sensitivity_with_tmoke( ...
                            model, STUDY_TAG, PARAM_N, na, ...
                            ALPHA_NAME, MSIGN_NAME, ...
                            alpha_start, alpha_step_coarse, alpha_stop, ...
                            M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                        Rows = [Rows; ...
                            Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), ...
                            hcey_grid_nm(ig), hau_grid_nm(ih), Sval, sign(Sval)]; %#ok<AGROW>

                        runs_done_global   = runs_done_global + numel(na)*2;
                        coarse_done_points = coarse_done_points + 1;
                        done_points_global = done_points_global + 1;

                        iter_global = iter_global + 1;
                        frac_global = iter_global / max(1,total_pts_est_global);

                        fprintf(['COARSE | it=%4d/%4d (%.1f%%) | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                 ' | S=%.6f deg/RIU (sign=%+d)\n'], ...
                            iter_global, total_pts_est_global, 100*frac_global, ...
                            Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), ...
                            hcey_grid_nm(ig), hau_grid_nm(ih), Sval, sign(Sval));

                        % ---- plot live: (1) TMOKE(α) e (5) Heatmap TMOKE(α,n) ----
                        geomLabel = sprintf('COARSE it=%d | Ldom=%.0f Lden=%.0f hsi=%.0f hcey=%.1f hau=%.1f', ...
                            iter_global, Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), hcey_grid_nm(ig), hau_grid_nm(ih));
                        livePlot_TMOKE_and_Heatmap(alpha_grid, na, TM_mat, geomLabel);

                        % ---- plot de S vs iteração (se quiser manter) ----
                        updateSensitivityPlot(iter_global, Sval, 'COARSE');

                        if mod(coarse_done_points, CKPT_EVERY_POINTS) == 0
                            CKPT.stage            = 'COARSE';
                            CKPT.done_points      = coarse_done_points;
                            CKPT.runs_done_global = runs_done_global;
                            CKPT.Rows             = Rows;
                            CKPT.Tcoarse          = [];
                            CKPT.Tfine            = Tfine;
                            CKPT.bestStruct       = bestStruct;
                            CKPT.iter_global      = iter_global;
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

    CKPT.stage            = 'COARSE';
    CKPT.done_points      = coarse_done_points;
    CKPT.runs_done_global = runs_done_global;
    CKPT.Tcoarse          = Tcoarse;
    CKPT.Tfine            = Tfine;
    CKPT.bestStruct       = bestStruct;
    CKPT.iter_global      = iter_global;
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
        if isfield(CKPT,'iter_global'), iter_global = CKPT.iter_global; end
        fprintf('>>> RESUME FINE: já tinha %d pontos.\n', fine_done_points);
    end

    flat_idx = 0;

    for s = 1:height(seeds)
        Hlist = (seeds.h_au_nm(s)    - fine_hau_delta)  : fine_hau_step  : (seeds.h_au_nm(s)    + fine_hau_delta);
        Glist = (seeds.h_ceyig_nm(s) - fine_hcey_delta) : fine_hcey_step : (seeds.h_ceyig_nm(s) + fine_hcey_delta);
        Llist = (seeds.L_domain_nm(s)- fine_Ldom_delta) : fine_Ldom_step : (seeds.L_domain_nm(s)+ fine_Ldom_delta);
        Dlist = (seeds.l_dente_nm(s) - fine_Lden_delta) : fine_Lden_step : (seeds.l_dente_nm(s) + fine_Lden_delta);
        Slist = (seeds.h_si_nm(s)    - fine_hsi_delta)  : fine_hsi_step  : (seeds.h_si_nm(s)    + fine_hsi_delta);

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

                            [Sval, alpha_peaks, alpha_grid, TM_mat] = evaluate_sensitivity_with_tmoke( ...
                                model, STUDY_TAG, PARAM_N, na, ...
                                ALPHA_NAME, MSIGN_NAME, ...
                                alpha_start, alpha_step_fine, alpha_stop, ...
                                M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                            FineRows = [FineRows; ...
                                Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), Sval, sign(Sval)]; %#ok<AGROW>

                            runs_done_global = runs_done_global + numel(na)*2;
                            fine_done_points = fine_done_points + 1;
                            done_points_global = done_points_global + 1;

                            iter_global = iter_global + 1;
                            frac_global = iter_global / max(1,total_pts_est_global);

                            fprintf(['FINE   | it=%4d/%4d (%.1f%%) | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                     ' | S=%.6f deg/RIU (sign=%+d)\n'], ...
                                iter_global, total_pts_est_global, 100*frac_global, ...
                                Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), Sval, sign(Sval));

                            geomLabel = sprintf('FINE it=%d | Ldom=%.0f Lden=%.0f hsi=%.0f hcey=%.1f hau=%.1f', ...
                                iter_global, Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih));
                            livePlot_TMOKE_and_Heatmap(alpha_grid, na, TM_mat, geomLabel);

                            updateSensitivityPlot(iter_global, Sval, 'FINE');

                            if mod(fine_done_points, CKPT_EVERY_POINTS) == 0
                                CKPT.stage            = 'FINE';
                                CKPT.done_points      = fine_done_points;
                                CKPT.runs_done_global = runs_done_global;
                                CKPT.FineRows         = FineRows;
                                CKPT.seeds            = seeds;
                                CKPT.Tcoarse          = Tcoarse;
                                CKPT.Tfine            = [];
                                CKPT.bestStruct       = bestStruct;
                                CKPT.iter_global      = iter_global;
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
    CKPT.iter_global      = iter_global;
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

fprintf('\n===== MELHOR CONFIGURAÇÃO POR |S| (COARSE+FINE) =====\n');
disp(bestRow);

CKPT.stage            = 'VALID';
CKPT.done_points      = 0;
CKPT.runs_done_global = runs_done_global;
CKPT.Tcoarse          = Tcoarse;
CKPT.Tfine            = Tfine;
CKPT.bestStruct       = bestStruct;
CKPT.iter_global      = iter_global;
save(CKPT_FILE,'CKPT','-v7.3');

% =====================================================================
%                 VALID: curva densa + guarda dados p/ gráficos
% =====================================================================
fprintf('\n===== STAGE VALID (reavalia melhor ponto) =====\n');

setParamNm(model, PARAM_LDOM, bestStruct.Ldom_best);
setParamNm(model, PARAM_LDEN, bestStruct.Lden_best);
setParamNm(model, PARAM_HSI,  bestStruct.hsi_best);
setParamNm(model, PARAM_HCEY, bestStruct.hcey_best);
setParamNm(model, PARAM_HAU,  bestStruct.hau_best);

alpha_peaks  = zeros(size(na));
TM_all       = cell(size(na));
alpha_all    = cell(size(na));
Rplus_all    = cell(size(na));
Rminus_all   = cell(size(na));
TMmax_abs    = zeros(size(na));

for i = 1:numel(na)
    n_val = na(i);
    model.param.set(PARAM_N, sprintf('%.12g', n_val));

    [alpha_deg, Rplus, Rminus, TM] = solveAndGetRplusRminus( ...
        model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
        alpha_start, alpha_step_dense, alpha_stop, ...
        M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

    [tmmax, k] = max(abs(TM));
    alpha_peaks(i) = alpha_deg(k);

    TM_all{i}     = TM;
    alpha_all{i}  = alpha_deg;
    Rplus_all{i}  = Rplus;
    Rminus_all{i} = Rminus;
    TMmax_abs(i)  = tmmax;

    runs_done_global = runs_done_global + 2;
end

p = polyfit(na, alpha_peaks, 1);
S_dense = p(1);

fprintf('alpha_peak (deg) = ['); fprintf(' %.4f', alpha_peaks); fprintf(' ]\n');
fprintf('Sensibilidade (VALID) S ≈ %.6f deg/RIU\n', S_dense);

% =====================================================================
%                 GRÁFICOS FINAIS (TODOS)
% =====================================================================

% (A) (1) TMOKE(α) por n (melhor ponto)
figure('Name','(1) TMOKE vs alpha (VALID)','NumberTitle','off');
hold on; grid on;
for i = 1:numel(na)
    plot(alpha_all{i}, TM_all{i}, 'LineWidth',1.2, ...
        'DisplayName', sprintf('n = %.2f', na(i)));
end
xlabel('\alpha [deg]');
ylabel('TMOKE(\alpha)');
title('TMOKE(\alpha) no melhor ponto (VALID)');
legend('Location','best');

% (B) (2) alpha_peak vs n (com fit)
figure('Name','(2) alpha_{peak} vs n (VALID)','NumberTitle','off');
hold on; grid on;
plot(na, alpha_peaks, 'o','LineWidth',1.5,'DisplayName','dados (\alpha_{peak})');
na_fit    = linspace(min(na), max(na), 100);
alpha_fit = polyval(p, na_fit);
plot(na_fit, alpha_fit, '-','LineWidth',1.2, ...
    'DisplayName', sprintf('ajuste linear (S = %.3f deg/RIU)', S_dense));
xlabel('Índice de refração n');
ylabel('\alpha_{peak} [deg]');
title('\alpha_{peak} em função de n (VALID)');
legend('Location','best');

% (C) (3) |TMOKE|max vs n
figure('Name','(3) |TMOKE|_{max} vs n (VALID)','NumberTitle','off');
grid on; hold on;
plot(na, TMmax_abs, 'o-','LineWidth',1.5);
xlabel('Índice de refração n');
ylabel('|TMOKE|_{max}');
title('Magnitude do pico |TMOKE|_{max} vs n (VALID)');

% (D) (5) Heatmap TMOKE(α,n) (VALID)
% Monta matriz TM_mat_valid: linhas = n, colunas = alpha
alpha_grid_valid = alpha_all{1};
TM_mat_valid = zeros(numel(na), numel(alpha_grid_valid));
for i = 1:numel(na)
    % assume mesma grade de alpha no VALID
    TM_mat_valid(i,:) = TM_all{i}(:).';
end

figure('Name','(5) Heatmap TMOKE(alpha,n) (VALID)','NumberTitle','off');
imagesc(alpha_grid_valid, na, TM_mat_valid);
axis xy; colorbar;
xlabel('\alpha [deg]');
ylabel('n');
title('Heatmap TMOKE(\alpha,n) no melhor ponto (VALID)');

% (E) (6) |S| histogram (COARSE e FINE)
figure('Name','(6) Histograma |S| (COARSE/FINE)','NumberTitle','off');
hold on; grid on;
histogram(abs(Tcoarse.S_deg_per_RIU), 'DisplayName','COARSE |S|');
if ~isempty(Tfine)
    histogram(abs(Tfine.S_deg_per_RIU), 'DisplayName','FINE |S|');
end
xlabel('|S| [deg/RIU]');
ylabel('contagem');
title('Distribuição de |S|');
legend('Location','best');

% (F) (7) Scatter |S| vs cada parâmetro (COARSE+FINE)
All = AllCand;
figure('Name','(7) |S| vs parâmetros (COARSE+FINE)','NumberTitle','off');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

nexttile; grid on; scatter(All.L_domain_nm, abs(All.S_deg_per_RIU), 30, 'filled');
xlabel('L\_domain [nm]'); ylabel('|S|'); title('|S| vs L\_domain');

nexttile; grid on; scatter(All.l_dente_nm,  abs(All.S_deg_per_RIU), 30, 'filled');
xlabel('l\_dente [nm]'); ylabel('|S|'); title('|S| vs l\_dente');

nexttile; grid on; scatter(All.h_si_nm,     abs(All.S_deg_per_RIU), 30, 'filled');
xlabel('h\_si [nm]'); ylabel('|S|'); title('|S| vs h\_si');

nexttile; grid on; scatter(All.h_ceyig_nm,  abs(All.S_deg_per_RIU), 30, 'filled');
xlabel('h\_ceyig [nm]'); ylabel('|S|'); title('|S| vs h\_ceyig');

nexttile; grid on; scatter(All.h_au_nm,     abs(All.S_deg_per_RIU), 30, 'filled');
xlabel('h\_au [nm]'); ylabel('|S|'); title('|S| vs h\_au');

nexttile; axis off;
text(0,0.7,sprintf('Best |S| (FINE): %.3f deg/RIU', abs(bestStruct.S_best)),'FontSize',12);
text(0,0.5,sprintf('Best S (VALID):  %.3f deg/RIU', abs(S_dense)),'FontSize',12);

% --------------------------- CSVs ------------------------------------
writetable(Tcoarse, fullfile(homeDir,'tmoke_sens_5D_coarse_TEST.csv'));
if ~isempty(Tfine)
    writetable(Tfine, fullfile(homeDir,'tmoke_sens_5D_fine_TEST.csv'));
end

BestSensRow = bestRow;
BestSensRow.alpha_peaks_deg         = {alpha_peaks};
BestSensRow.S_dense_deg_per_RIU     = S_dense;
writetable(BestSensRow, fullfile(homeDir,'tmoke_sens_5D_best_TEST.csv'));

elapsed_total = toc(t0);
fprintf('\n===== SUMMARY (TESTE) =====\n');
fprintf('Best 5D:\n');
fprintf('  L_domain = %4.0f nm\n', bestStruct.Ldom_best);
fprintf('  l_dente  = %4.0f nm\n', bestStruct.Lden_best);
fprintf('  h_si     = %3.0f nm\n', bestStruct.hsi_best);
fprintf('  h_ceyig  = %.3f nm\n', bestStruct.hcey_best);
fprintf('  h_au     = %.3f nm\n', bestStruct.hau_best);
fprintf('  S (FINE)  ≈ %.6f deg/RIU\n', bestStruct.S_best);
fprintf('  S (VALID) ≈ %.6f deg/RIU\n', S_dense);
fprintf('Runs totais (aprox) ≈ %d\n', runs_done_global);
fprintf('Tempo total: %.1f s\n', elapsed_total);

% =====================================================================
%                           FUNÇÕES LOCAIS
% =====================================================================

function setParamNm(mdl, name, val_nm)
    mdl.param.set(name, sprintf('%.12g[nm]', val_nm));
end

function [Sval, alpha_peaks, alpha_grid, TM_mat] = evaluate_sensitivity_with_tmoke( ...
        mdl, studyTag, PARAM_N, na, ...
        alphaName, mName, aStartDeg, aStepDeg, aStopDeg, ...
        mPlusStr, mMinusStr, ttagPlus, ttagMinus)
% Calcula S e também retorna TMOKE(α) pra cada n:
%   - alpha_grid: vetor de alpha (deg)
%   - TM_mat: matriz [numel(na) x numel(alpha_grid)] com TMOKE
    alpha_peaks = zeros(size(na));
    alpha_grid = [];
    TM_mat = [];

    for i = 1:numel(na)
        n_val = na(i);
        mdl.param.set(PARAM_N, sprintf('%.12g', n_val));

        [alpha_deg, ~, ~, TM] = solveAndGetRplusRminus( ...
            mdl, studyTag, alphaName, mName, ...
            aStartDeg, aStepDeg, aStopDeg, ...
            mPlusStr, mMinusStr, ttagPlus, ttagMinus);

        if isempty(TM)
            error('TMOKE vazio na avaliação de sensibilidade.');
        end

        if isempty(alpha_grid)
            alpha_grid = alpha_deg(:).';
            TM_mat = zeros(numel(na), numel(alpha_grid));
        end

        TM_mat(i,:) = TM(:).';

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

    % m = +1
    redirectAllNumericalsToTable(mdl, ttagPlus);
    clearTable(mdl, ttagPlus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, ...
                     aStartDeg, aStepDeg, aStopDeg, mPlusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha1_deg, R1] = readAlphaAndRFromNamedTable(mdl, ttagPlus, Npts);

    % m = -1
    redirectAllNumericalsToTable(mdl, ttagMinus);
    clearTable(mdl, ttagMinus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, ...
                     aStartDeg, aStepDeg, aStopDeg, mMinusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha2_deg, R2] = readAlphaAndRFromNamedTable(mdl, ttagMinus, Npts);

    assert(numel(alpha1_deg)==Npts && numel(alpha2_deg)==Npts, 'Unexpected α sweep length.');
    assert(max(abs(alpha1_deg - alpha2_deg)) < 1e-9, 'α grids differ between m=+1 and m=-1.');

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

function updateSensitivityPlot(iter_idx, Sval, stage)
% Plot online: S vs iteração global, separando COARSE/FINE
    persistent fig ax hCoarse hFine inited
    if isempty(inited) || ~isgraphics(fig)
        fig = figure('Name','Sensibilidade vs Iteração','NumberTitle','off');
        ax  = axes('Parent',fig); hold(ax,'on'); grid(ax,'on');
        xlabel(ax,'Iteração global (ponto de geometria)');
        ylabel(ax,'S [deg/RIU]');
        title(ax,'Evolução da sensibilidade ao longo da busca');
        hCoarse = plot(ax, nan, nan, 'o-','DisplayName','COARSE');
        hFine   = plot(ax, nan, nan, 'x-','DisplayName','FINE');
        legend(ax,'Location','best');
        inited = true;
    end

    switch upper(stage)
        case 'COARSE'
            x = get(hCoarse,'XData'); y = get(hCoarse,'YData');
            x = [x, iter_idx];        y = [y, Sval];
            set(hCoarse,'XData',x,'YData',y);
        case 'FINE'
            x = get(hFine,'XData'); y = get(hFine,'YData');
            x = [x, iter_idx];      y = [y, Sval];
            set(hFine,'XData',x,'YData',y);
    end
    drawnow limitrate;
end

function livePlot_TMOKE_and_Heatmap(alpha_deg, na, TM_mat, geomLabel)
% A cada iteração plota:
% (1) TMOKE(alpha) para cada n
% (5) Heatmap TMOKE(alpha,n)
    persistent fig ax1 ax2 hLines hImg inited

    if isempty(inited) || ~isgraphics(fig)
        fig = figure('Name','LIVE: (1) TMOKE(alpha) + (5) Heatmap','NumberTitle','off');
        ax1 = subplot(1,2,1,'Parent',fig); grid(ax1,'on'); hold(ax1,'on');
        ax2 = subplot(1,2,2,'Parent',fig);

        hLines = gobjects(numel(na),1);
        for i = 1:numel(na)
            hLines(i) = plot(ax1, alpha_deg, TM_mat(i,:), 'LineWidth', 1.2, ...
                'DisplayName', sprintf('n=%.2f', na(i)));
        end
        xlabel(ax1,'\alpha [deg]'); ylabel(ax1,'TMOKE(\alpha)');
        title(ax1,'(1) TMOKE(\alpha) por n');
        legend(ax1,'Location','best');

        hImg = imagesc(ax2, alpha_deg, na, TM_mat);
        axis(ax2,'xy'); colorbar(ax2);
        xlabel(ax2,'\alpha [deg]'); ylabel(ax2,'n');
        title(ax2,'(5) Heatmap TMOKE(\alpha,n)');

        inited = true;
    else
        for i = 1:numel(na)
            set(hLines(i),'XData',alpha_deg,'YData',TM_mat(i,:));
        end
        set(hImg,'XData',alpha_deg,'YData',na,'CData',TM_mat);
        axis(ax2,'xy');
    end

    sgtitle(fig, geomLabel, 'Interpreter','none');
    drawnow limitrate;
end
