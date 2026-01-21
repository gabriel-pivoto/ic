%% =====================================================================
% MATLAB + COMSOL — 5D Hierarchical Search (TMOKE + Sensibilidade S)
% ---------------------------------------------------------------------
% OBJETIVO (em 1 frase):
%   Encontrar a melhor geometria 5D que seja boa ao mesmo tempo em:
%     (A) |TMOKE|max  e  (B) Sensibilidade angular S = d(alpha_peak)/dn [deg/RIU]
%
% COMO ESTE SCRIPT “JUNTA” OS DOIS CÓDIGOS:
% ---------------------------------------------------------------------
% Antes você tinha 2 pipelines separados:
%   1) Otimização por TMOKE (max |TMOKE|)   -> barato (2 runs por ponto)
%   2) Otimização por sensibilidade S(n)    -> caro (>= 2*#n runs por ponto)
%
% Aqui, em CADA ponto de geometria (COARSE/FINE/SUPER), nós calculamos:
%   - Métrica TMOKE:  maxAbsTMOKE_base = max(|TMOKE(alpha)|) para n_base
%   - Métrica Sens.:  S_est = (alpha_peak(n2) - alpha_peak(n1)) / (n2 - n1)
%                    usando APENAS 2 valores de n (rápido) durante a busca.
%
% Depois (VALID), no melhor ponto final, calculamos S “de verdade” com 3 n:
%   - n_valid = [1.33 1.36 1.39] (ou o que você quiser)
%   - S_dense = slope do fit linear alpha_peak vs n (deg/RIU)
%
% SELEÇÃO “MELHOR DOS DOIS” (trade-off):
% ---------------------------------------------------------------------
% Em cada estágio (COARSE/FINE/SUPER), nós ranqueamos os candidatos por:
%   - rank_TM = rank de |TMOKE|max (maior é melhor)
%   - rank_S  = rank de |S_est|     (maior é melhor)
%   - score_comb = rank_TM + rank_S   (menor score = melhor trade-off)
%
% Assim você obtém:
%   - best_tradeoff: melhor compromisso (TMOKE e S juntos)
%   - best_TMOKE: melhor só em TMOKE
%   - best_S: melhor só em sensibilidade
%
% FLUXOGRAMA (alto nível)
% ---------------------------------------------------------------------
% [Start] -> Load .mph -> (Resume checkpoint?)
%   |
%   v
% [COARSE]  (varre grade ampla; para cada ponto calcula TMOKE e S_est)
%   |
%   v
% Seleciona seeds TOP-K (pelo score_comb)
%   |
%   v
% [FINE]    (refina parâmetros + janela alpha; recalcula TMOKE e S_est)
%   |
%   v
% Seleciona seeds TOP-K (pelo score_comb)
%   |
%   v
% [SUPER]   (refino final; recalcula TMOKE e S_est)
%   |
%   v
% Escolhe best_tradeoff + best_TMOKE + best_S
%   |
%   v
% [VALID]   (no best_tradeoff: curva densa para n_valid e calcula S_dense)
%   |
%   v
% [FULL]    (curva final (baseline) e exporta CSV/XLSX)
%   |
%   v
% [SAVE FIGURES] (exporta todas as figures abertas)
%   |
%   v
% [End]
% =====================================================================

clear; clc; close all; format long; tic
import com.comsol.model.*;
import com.comsol.model.util.*;

%% --------------------------- Paths/Project -------------------------
homeDir  = 'C:\Users\gabri\Documents\projetoIC';
mph_file = fullfile(homeDir,'usandoMatlab.mph');
addpath(genpath(homeDir));

%% ------------------------ Run budget (opcional) ---------------------
% Mantém o “espírito” do script TMOKE (< 20k). Agora, como calculamos S_est
% (2 valores de n) em cada ponto, o custo por ponto dobra.
MAX_RUNS = 20000;

%% ------------------------ Checkpoint config -------------------------
CKPT_DIR          = fullfile(homeDir,'checkpoints');
if ~exist(CKPT_DIR,'dir'), mkdir(CKPT_DIR); end
CKPT_FILE         = fullfile(CKPT_DIR,'tmoke_sens_5d_checkpoint.mat');
PROGRESS_XLSX     = fullfile(CKPT_DIR,'tmoke_sens_progress.xlsx');

CKPT_EVERY_POINTS = 10;      % salva a cada 10 pontos (não a cada run)
points_since_ckpt = 0;

%% ------------------------ Resume / Runlog ---------------------------
RESUME = false; CKPT = struct; RESUME_STAGE = "";
if isfile(CKPT_FILE)
    try
        S = load(CKPT_FILE,'CKPT');
        CKPT = S.CKPT;
        RESUME = true;
        RESUME_STAGE = string(CKPT.stage);
        fprintf('>>> RESUME detected | stage=%s | done_points=%d | runs_done_global=%d\n', ...
            RESUME_STAGE, CKPT.done_points, CKPT.runs_done_global);
    catch
        warning('Checkpoint exists but could not be loaded. Starting fresh.');
    end
end

diary(fullfile(CKPT_DIR, sprintf('runlog_%s.txt', datestr(now,'yyyymmdd_HHMMss'))));

%% ------------------------ Model Tags/Params -------------------------
STUDY_TAG   = 'std1';

PARAM_HAU   = 'h_au';
PARAM_HCEY  = 'h_ceyig';
PARAM_LDOM  = 'L_domain';
PARAM_LDEN  = 'l_dente';
PARAM_HSI   = 'h_si';
PARAM_N     = 'n';               % <<< vem do seu script de sensibilidade

ALPHA_NAME  = 'alpha';
MSIGN_NAME  = 'm';
M_PLUS      = '1';
M_MINUS     = '-1';

RPLUS_TABLE_TAG  = 'tblRplus';
RMINUS_TABLE_TAG = 'tblRminus';

%% --------------------- Índices de refração (n) ----------------------
% Durante a busca (COARSE/FINE/SUPER):
% - usamos 2 pontos de n para estimar S rápido (S_est)
% - n_base também é o usado para a métrica de TMOKE (maxAbsTMOKE_base)
n_fast   = [1.33, 1.39];         % <<< S_est rápido (2 pontos)
n_base   = n_fast(1);            % <<< TMOKE “base” (pode trocar se quiser)

% Na validação final (VALID):
% - usamos 3 pontos (ou mais) para S_dense mais confiável
n_valid  = [1.33, 1.36, 1.39];

%% ------------------------ COARSE grids (exact) ----------------------
Ldomain_grid_nm = 800:50:850;         % 3
ldente_grid_nm  = 500:50:600;         % 3
hsi_grid_nm     = [220, 240, 260];    % 3
hcey_grid_nm    = [100, 140];         % 2
hau_grid_nm     = 20:10:60;           % 5
% COARSE points = 270

%% ---------------------- Alpha ranges/steps --------------------------
alpha_coarse_rng = [0, 1.0, 89];  % COARSE sweep completo
alpha_fine_step  = 0.1;           % FINE step
alpha_super_step = 0.01;          % SUPERFINE step
alpha_dense_step = 0.01;          % VALID step
alpha_full_step  = 0.01;          % FULL step

%% ----------------------- TOP-K strategy ----------------------------
TOPK_COARSE = 1;
TOPK_FINE   = 1;

%% --------------------- Refinement windows --------------------------
fine_hau_delta   = 2;    fine_hau_step   = 1;
fine_hcey_delta  = 5;    fine_hcey_step  = 1;
fine_Ldom_delta  = 10;   fine_Ldom_step  = 5;
fine_Lden_delta  = 10;   fine_Lden_step  = 5;
fine_hsi_delta   = 5;    fine_hsi_step   = 5;

super_hau_delta  = 1;    super_hau_step  = 0.5;
super_hcey_delta = 2;    super_hcey_step = 1;
super_Ldom_delta = 4;    super_Ldom_step = 2;
super_Lden_delta = 4;    super_Lden_step = 2;
super_hsi_delta  = 2;    super_hsi_step  = 1;

%% ------------------ Alpha windows (centered at last peak) ----------
FINE_ALPHA_HALFSPAN   = 5;     % TMOKE window nominal no FINE
SUPER_ALPHA_HALFSPAN  = 2;     % TMOKE window nominal no SUPER

% Para sensibilidade, é comum o alpha_peak deslocar um pouco com n.
% Então usamos uma janela levemente maior para calcular S_est com segurança.
FINE_ALPHA_HALFSPAN_SENS  = FINE_ALPHA_HALFSPAN  + 1;   % ±6°
SUPER_ALPHA_HALFSPAN_SENS = SUPER_ALPHA_HALFSPAN + 2;   % ±4°

%% ----------------------- Outputs / Plots ----------------------------
SAVE_SNAPSHOT = true;
MAKE_PLOTS    = true;
PLOT_LIVE     = true;

%% ----------------------- Save figures (FINAL) -----------------------
SAVE_FIGS   = true;
FIG_FORMATS = {'png','pdf'};
RUN_STAMP   = datestr(now,'yyyy-mm-dd_HHMMSS');
FIG_DIR     = fullfile(homeDir,'plots', ['tmoke_sens_5d_' RUN_STAMP]);
if SAVE_FIGS && ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

%% ------------------------- Load model ------------------------------
ModelUtil.clear;
model = mphload(mph_file);
ModelUtil.showProgress(true);

%% ==================================================================
%                    GLOBAL PLANNING (runs / ETA)
% ==================================================================
t0 = tic;
runs_done_global = 0;
is_total_exact   = false;

nL = numel(Ldomain_grid_nm);
nD = numel(ldente_grid_nm);
nS = numel(hsi_grid_nm);
nG = numel(hcey_grid_nm);
nH = numel(hau_grid_nm);

% Cada n precisa de 2 runs (m=±1). Durante a busca usamos n_fast (2 valores).
runs_per_point_search = 2 * numel(n_fast);     % 4 runs por ponto

coarse_total_pts  = nL*nD*nS*nG*nH;
coarse_total_runs = runs_per_point_search * coarse_total_pts;

% Estimativas do FINE (antes do clamp exato)
fine_H = 1 + floor(2*fine_hau_delta   / fine_hau_step);
fine_G = 1 + floor(2*fine_hcey_delta  / fine_hcey_step);
fine_L = 1 + floor(2*fine_Ldom_delta  / fine_Ldom_step);
fine_D = 1 + floor(2*fine_Lden_delta  / fine_Lden_step);
fine_S = 1 + floor(2*fine_hsi_delta   / fine_hsi_step);
fine_pts_per_seed_est = fine_L*fine_D*fine_S*fine_G*fine_H;
fine_total_pts_est    = TOPK_COARSE * fine_pts_per_seed_est;
fine_total_runs_est   = runs_per_point_search * fine_total_pts_est;

% SUPERFINE (exato upfront)
super_H = 1 + floor(2*super_hau_delta   / super_hau_step);
super_G = 1 + floor(2*super_hcey_delta  / super_hcey_step);
super_L = 1 + floor(2*super_Ldom_delta  / super_Ldom_step);
super_D = 1 + floor(2*super_Lden_delta  / super_Lden_step);
super_S = 1 + floor(2*super_hsi_delta   / super_hsi_step);
super_pts_per_seed = super_L*super_D*super_S*super_G*super_H;
super_total_pts    = TOPK_FINE * super_pts_per_seed;
super_total_runs   = runs_per_point_search * super_total_pts;

% EXTRAS:
% - VALID: 2 runs por n_valid (m=±1) => 2*len(n_valid)
% - snapshot: 2 runs (m=±1) num alpha fixo no n_base
% - FULL: 2 runs (m=±1) sweep completo no n_base
extras_runs_fixed = 2*numel(n_valid) + (SAVE_SNAPSHOT*2) + 2;

grand_total_target = coarse_total_runs + fine_total_runs_est + super_total_runs + extras_runs_fixed;
is_total_exact = false;

fprintf('START\n');
fprintf('  COARSE (exact):    %d runs\n', coarse_total_runs);
fprintf('  FINE   (estimate): %d runs (TOP-%d)\n', fine_total_runs_est, TOPK_COARSE);
fprintf('  SUPER  (exact):    %d runs (TOP-%d)\n', super_total_runs, TOPK_FINE);
fprintf('  EXTRAS (exact):    %d runs\n', extras_runs_fixed);
fprintf('  TOTAL  (estimate): %d runs\n\n', grand_total_target);

%% ==================================================================
%                               COARSE
% ==================================================================
Tcoarse = [];
seeds   = [];

skip_COARSE = RESUME && any(strcmp(RESUME_STAGE, ["FINE","SUPER","VALID","FULL"]));
if skip_COARSE
    if isfield(CKPT.payload,'Tcoarse'), Tcoarse = CKPT.payload.Tcoarse; end
    if isfield(CKPT.payload,'seeds'),   seeds   = CKPT.payload.seeds;   end
    fprintf('SKIP COARSE -> restored Tcoarse (rows=%d), seeds (rows=%d)\n', ...
        size(Tcoarse,1), height(seeds));
else
    fprintf('STAGE COARSE — EXACT: %d runs\n', coarse_total_runs);
    stage_runs_start = runs_done_global;
    stage_total_runs = coarse_total_runs;
    stage_t0 = tic;

    % Colunas:
    % [Ldom, Lden, hsi, hcey, hau,
    %  maxAbsTMOKE_base, alpha_peak_base, TMOKE_at_peak_base,
    %  alpha_peak_n2, S_est]
    Rows = [];
    coarse_point_idx = 0;

    if RESUME && RESUME_STAGE=="COARSE"
        runs_done_global = CKPT.runs_done_global;
        if isfield(CKPT.payload,'Rows'), Rows = CKPT.payload.Rows; end
        coarse_point_idx = CKPT.done_points;
        fprintf('>>> COARSE resume: skipping %d points already computed.\n', coarse_point_idx);
    end

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

                        coarse_point_idx = coarse_point_idx + 1;
                        if RESUME && RESUME_STAGE=="COARSE" && coarse_point_idx <= CKPT.done_points
                            continue;
                        end

                        % =======================
                        % (1) n = n_fast(1) (base)
                        % =======================
                        setParamScalar(model, PARAM_N, n_fast(1)); % n sem unidade
                        [a1, Rp1, Rm1, TM1] = solveAndGetRplusRminus( ...
                            model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                            alpha_coarse_rng(1), alpha_coarse_rng(2), alpha_coarse_rng(3), ...
                            M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                        [pk1, k1] = max(abs(TM1));
                        alpha1 = a1(k1);
                        tm_at_pk1 = TM1(k1);

                        if k1==1 || k1==numel(TM1)
                            logline('WARN COARSE (n=%.2f) peak at α-edge (α*=%.3f in [%.3f, %.3f])\n', ...
                                n_fast(1), alpha1, a1(1), a1(end));
                        end

                        if PLOT_LIVE
                            updateLivePlot('COARSE', ...
                                Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), ...
                                hcey_grid_nm(ig), hau_grid_nm(ih), ...
                                a1, TM1, alpha1, tm_at_pk1, Rp1, Rm1);
                        end

                        % ======================
                        % (2) n = n_fast(2) (p/ S_est)
                        % ======================
                        setParamScalar(model, PARAM_N, n_fast(2));
                        [a2, ~, ~, TM2] = solveAndGetRplusRminus( ...
                            model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                            alpha_coarse_rng(1), alpha_coarse_rng(2), alpha_coarse_rng(3), ...
                            M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                        [~, k2] = max(abs(TM2));
                        alpha2 = a2(k2);

                        % Sensibilidade rápida (2 pontos)
                        S_est = (alpha2 - alpha1) / (n_fast(2) - n_fast(1));

                        % Guarda ponto
                        Rows = [Rows; ...
                            Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), ...
                            hcey_grid_nm(ig), hau_grid_nm(ih), ...
                            pk1, alpha1, tm_at_pk1, ...
                            alpha2, S_est]; %#ok<AGROW>

                        % Custo do ponto: 2*n_fast runs (m=±1 para cada n)
                        runs_done_global = runs_done_global + runs_per_point_search;

                        % ETA
                        stage_done = runs_done_global - stage_runs_start;
                        [frac_stage, eta_stage]   = frac_eta(stage_done, stage_total_runs, toc(stage_t0));
                        [frac_global, eta_global] = frac_eta(runs_done_global, grand_total_target, toc(t0));

                        logline(['COARSE | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                 ' | |TM|=%.5f @ α=%.2f° | S_est=%.4f deg/RIU' ...
                                 ' [Stage %5.1f%% | ETA %s | Global %s%% | ETA %s]\n'], ...
                            Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), ...
                            hcey_grid_nm(ig), hau_grid_nm(ih), ...
                            pk1, alpha1, S_est, ...
                            100*frac_stage, fmt_time_long(eta_stage), fmt_pct(frac_global), fmt_eta_flag(eta_global, is_total_exact));

                        % Checkpoint
                        points_since_ckpt = points_since_ckpt + 1;
                        payload = struct('Rows',Rows);
                        points_since_ckpt = maybe_checkpoint( ...
                            'COARSE', CKPT_EVERY_POINTS, CKPT_FILE, PROGRESS_XLSX, ...
                            runs_done_global, points_since_ckpt, coarse_point_idx, payload, ...
                            [], [], [], [], []);
                    end
                end
            end
        end
    end

    Tcoarse = array2table(Rows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm', ...
         'maxAbsTMOKE_base','alpha_peak_base_deg','TMOKE_at_peak_base', ...
         'alpha_peak_n2_deg','S_est_deg_per_RIU'});

    % Seleção TOP-K por score combinado (rank TMOKE + rank |S|)
    seeds = selectTopK_tradeoff(Tcoarse, TOPK_COARSE, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');

    save_checkpoint(CKPT_FILE, 'FINE', runs_done_global, 0, struct('Tcoarse', Tcoarse, 'seeds', seeds));
    write_progress_xlsx(PROGRESS_XLSX, 'coarse', Tcoarse);
end

%% ==================================================================
%                        FINE PLANNING (EXACT NOW)
% ==================================================================
if isempty(seeds) && RESUME && RESUME_STAGE=="FINE"
    if isfield(CKPT.payload,'seeds')
        seeds = CKPT.payload.seeds;
    else
        error('Checkpoint at FINE stage has no "seeds". Rerun COARSE.');
    end
end

% Conta exata do FINE (com clamp)
fine_total_pts = 0;
for s = 1:height(seeds)
    Hlist = max(min(hau_grid_nm),  seeds.h_au_nm(s)   - fine_hau_delta)  : fine_hau_step  : min(max(hau_grid_nm),  seeds.h_au_nm(s)   + fine_hau_delta);
    Glist = max(min(hcey_grid_nm), seeds.h_ceyig_nm(s)- fine_hcey_delta) : fine_hcey_step : min(max(hcey_grid_nm), seeds.h_ceyig_nm(s)+ fine_hcey_delta);
    Llist = max(min(Ldomain_grid_nm), seeds.L_domain_nm(s)- fine_Ldom_delta) : fine_Ldom_step : min(max(Ldomain_grid_nm), seeds.L_domain_nm(s)+ fine_Ldom_delta);
    Dlist = max(min(ldente_grid_nm),  seeds.l_dente_nm(s) - fine_Lden_delta) : fine_Lden_step : min(max(ldente_grid_nm),  seeds.l_dente_nm(s) + fine_Lden_delta);
    Slist = max(min(hsi_grid_nm),     seeds.h_si_nm(s)    - fine_hsi_delta)  : fine_hsi_step  : min(max(hsi_grid_nm),     seeds.h_si_nm(s)    + fine_hsi_delta);
    fine_total_pts = fine_total_pts + numel(Llist)*numel(Dlist)*numel(Slist)*numel(Glist)*numel(Hlist);
end
fine_total_runs = runs_per_point_search * fine_total_pts;

% Total global vira EXATO agora
grand_total_target = coarse_total_runs + fine_total_runs + super_total_runs + extras_runs_fixed;
is_total_exact = true;

fprintf('\nSTAGE FINE — EXACT: %d runs (TOP-%d)\n', fine_total_runs, height(seeds));
fprintf('GLOBAL TOTAL (EXACT): %d runs\n\n', grand_total_target);
if grand_total_target > MAX_RUNS
    error('Planned total exceeds %d runs (%d). Ajuste deltas/steps/TOPK.', MAX_RUNS, grand_total_target);
end

%% ==================================================================
%                                 FINE
% ==================================================================
Tfine = [];
seedsSuper = [];

if RESUME && any(strcmp(RESUME_STAGE, ["SUPER","VALID","FULL"]))
    if isfield(CKPT.payload,'Tfine'),      Tfine = CKPT.payload.Tfine; end
    if isfield(CKPT.payload,'seedsSuper'), seedsSuper = CKPT.payload.seedsSuper; end
    fprintf('SKIP FINE -> restored Tfine (rows=%d), seedsSuper (rows=%d)\n', ...
        size(Tfine,1), height(seedsSuper));
else
    FineRows = [];
    stage_runs_start = runs_done_global;
    stage_total_runs = fine_total_runs;
    stage_t0 = tic;

    fine_point_idx = 0;
    if RESUME && RESUME_STAGE=="FINE" && CKPT.done_points > 0
        runs_done_global = CKPT.runs_done_global;
        fine_point_idx   = CKPT.done_points;
        if isfield(CKPT.payload,'FineRows'), FineRows = CKPT.payload.FineRows; end
        if isfield(CKPT.payload,'seeds'),    seeds    = CKPT.payload.seeds;    end
        fprintf('>>> FINE resume: skipping %d points already computed.\n', fine_point_idx);
    end

    for s = 1:height(seeds)
        % Centro da janela em alpha vem do pico no n_base do estágio anterior
        aCenter = seeds.alpha_peak_base_deg(s);

        % Clamp dos parâmetros
        Hlist = max(min(hau_grid_nm),  seeds.h_au_nm(s)   - fine_hau_delta)  : fine_hau_step  : min(max(hau_grid_nm),  seeds.h_au_nm(s)   + fine_hau_delta);
        Glist = max(min(hcey_grid_nm), seeds.h_ceyig_nm(s)- fine_hcey_delta) : fine_hcey_step : min(max(hcey_grid_nm), seeds.h_ceyig_nm(s)+ fine_hcey_delta);
        Llist = max(min(Ldomain_grid_nm), seeds.L_domain_nm(s)- fine_Ldom_delta) : fine_Ldom_step : min(max(Ldomain_grid_nm), seeds.L_domain_nm(s)+ fine_Ldom_delta);
        Dlist = max(min(ldente_grid_nm),  seeds.l_dente_nm(s) - fine_Lden_delta) : fine_Lden_step : min(max(ldente_grid_nm),  seeds.l_dente_nm(s) + fine_Lden_delta);
        Slist = max(min(hsi_grid_nm),     seeds.h_si_nm(s)    - fine_hsi_delta)  : fine_hsi_step  : min(max(hsi_grid_nm),     seeds.h_si_nm(s)    + fine_hsi_delta);

        % Janela de alpha (um pouco maior para S_est)
        aStart = max(0,  aCenter - FINE_ALPHA_HALFSPAN_SENS);
        aStop  = min(89, aCenter + FINE_ALPHA_HALFSPAN_SENS);

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

                            fine_point_idx = fine_point_idx + 1;
                            if RESUME && RESUME_STAGE=="FINE" && fine_point_idx <= CKPT.done_points
                                continue;
                            end

                            % n = n_fast(1)
                            setParamScalar(model, PARAM_N, n_fast(1));
                            [a1, Rp1, Rm1, TM1] = solveAndGetRplusRminus( ...
                                model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                                aStart, alpha_fine_step, aStop, ...
                                M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                            [pk1, k1] = max(abs(TM1));
                            alpha1 = a1(k1);
                            tm_at_pk1 = TM1(k1);

                            if PLOT_LIVE
                                updateLivePlot('FINE', Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), ...
                                    a1, TM1, alpha1, tm_at_pk1, Rp1, Rm1);
                            end

                            % n = n_fast(2) (para S_est)
                            setParamScalar(model, PARAM_N, n_fast(2));
                            [a2, ~, ~, TM2] = solveAndGetRplusRminus( ...
                                model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                                aStart, alpha_fine_step, aStop, ...
                                M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                            [~, k2] = max(abs(TM2));
                            alpha2 = a2(k2);

                            S_est = (alpha2 - alpha1) / (n_fast(2) - n_fast(1));

                            FineRows = [FineRows; ...
                                Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), ...
                                pk1, alpha1, tm_at_pk1, alpha2, S_est]; %#ok<AGROW>

                            runs_done_global = runs_done_global + runs_per_point_search;

                            stage_done = runs_done_global - stage_runs_start;
                            [frac_stage, eta_stage]   = frac_eta(stage_done, stage_total_runs, toc(stage_t0));
                            [frac_global, eta_global] = frac_eta(runs_done_global, grand_total_target, toc(t0));
                            logline(['FINE   | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                     ' | |TM|=%.5f @ α=%.3f° | S_est=%.4f' ...
                                     ' [Stage %5.1f%% | ETA %s | Global %5.1f%% | ETA %s]\n'], ...
                                Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), ...
                                pk1, alpha1, S_est, ...
                                100*frac_stage, fmt_time_long(eta_stage), 100*frac_global, fmt_time_long(eta_global));

                            points_since_ckpt = points_since_ckpt + 1;
                            payload = struct('FineRows',FineRows, 'seeds', seeds);
                            points_since_ckpt = maybe_checkpoint( ...
                                'FINE', CKPT_EVERY_POINTS, CKPT_FILE, PROGRESS_XLSX, ...
                                runs_done_global, points_since_ckpt, fine_point_idx, payload, ...
                                Tcoarse, [], [], [], []);
                        end
                    end
                end
            end
        end
    end

    Tfine = array2table(FineRows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm', ...
         'maxAbsTMOKE_base','alpha_peak_base_deg','TMOKE_at_peak_base', ...
         'alpha_peak_n2_deg','S_est_deg_per_RIU'});

    seedsSuper = selectTopK_tradeoff(Tfine, TOPK_FINE, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');

    save_checkpoint(CKPT_FILE, 'SUPER', runs_done_global, 0, struct( ...
        'Tcoarse', Tcoarse, 'Tfine', Tfine, 'seedsSuper', seedsSuper));
    write_progress_xlsx(PROGRESS_XLSX, 'fine', Tfine);
end

%% ==================================================================
%                     SUPERFINE PLANNING (EXACT)
% ==================================================================
if isempty(seedsSuper) && RESUME && RESUME_STAGE=="SUPER"
    if isfield(CKPT.payload,'seedsSuper')
        seedsSuper = CKPT.payload.seedsSuper;
    else
        error('Checkpoint at SUPER stage has no "seedsSuper". Rerun FINE.');
    end
end

super_total_pts = 0;
for s = 1:height(seedsSuper)
    Hlist = (seedsSuper.h_au_nm(s)    - super_hau_delta)  : super_hau_step  : (seedsSuper.h_au_nm(s)    + super_hau_delta);
    Glist = (seedsSuper.h_ceyig_nm(s) - super_hcey_delta) : super_hcey_step : (seedsSuper.h_ceyig_nm(s) + super_hcey_delta);
    Llist = (seedsSuper.L_domain_nm(s)- super_Ldom_delta) : super_Ldom_step : (seedsSuper.L_domain_nm(s)+ super_Ldom_delta);
    Dlist = (seedsSuper.l_dente_nm(s) - super_Lden_delta) : super_Lden_step : (seedsSuper.l_dente_nm(s) + super_Lden_delta);
    Slist = (seedsSuper.h_si_nm(s)    - super_hsi_delta)  : super_hsi_step  : (seedsSuper.h_si_nm(s)    + super_hsi_delta);
    super_total_pts = super_total_pts + numel(Llist)*numel(Dlist)*numel(Slist)*numel(Glist)*numel(Hlist);
end
super_total_runs = runs_per_point_search * super_total_pts;

fprintf('\nSTAGE SUPER — EXACT: %d runs (TOP-%d)\n', super_total_runs, height(seedsSuper));
fprintf('GLOBAL TOTAL (EXACT): %d runs\n\n', grand_total_target);

%% ==================================================================
%                             SUPERFINE
% ==================================================================
Tsuper = [];
best_tradeoff = table(); best_TMOKE = table(); best_S = table();

if RESUME && any(strcmp(RESUME_STAGE, ["VALID","FULL"]))
    if isfield(CKPT.payload,'Tsuper'), Tsuper = CKPT.payload.Tsuper; end
    if isfield(CKPT.payload,'best_tradeoff'), best_tradeoff = CKPT.payload.best_tradeoff; end
    if isfield(CKPT.payload,'best_TMOKE'),    best_TMOKE    = CKPT.payload.best_TMOKE;    end
    if isfield(CKPT.payload,'best_S'),        best_S        = CKPT.payload.best_S;        end
    fprintf('SKIP SUPER -> restored Tsuper (rows=%d) and best selections.\n', size(Tsuper,1));
else
    SuperRows = [];
    stage_runs_start = runs_done_global;
    stage_total_runs = super_total_runs;
    stage_t0 = tic;

    super_point_idx = 0;
    if RESUME && RESUME_STAGE=="SUPER" && CKPT.done_points > 0
        runs_done_global = CKPT.runs_done_global;
        super_point_idx  = CKPT.done_points;
        if isfield(CKPT.payload,'SuperRows'),  SuperRows  = CKPT.payload.SuperRows;  end
        if isfield(CKPT.payload,'seedsSuper'), seedsSuper = CKPT.payload.seedsSuper; end
        fprintf('>>> SUPER resume: skipping %d points already computed.\n', super_point_idx);
    end

    for s = 1:height(seedsSuper)
        aCenter = seedsSuper.alpha_peak_base_deg(s);

        Hlist = (seedsSuper.h_au_nm(s)    - super_hau_delta)  : super_hau_step  : (seedsSuper.h_au_nm(s)    + super_hau_delta);
        Glist = (seedsSuper.h_ceyig_nm(s) - super_hcey_delta) : super_hcey_step : (seedsSuper.h_ceyig_nm(s) + super_hcey_delta);
        Llist = (seedsSuper.L_domain_nm(s)- super_Ldom_delta) : super_Ldom_step : (seedsSuper.L_domain_nm(s)+ super_Ldom_delta);
        Dlist = (seedsSuper.l_dente_nm(s) - super_Lden_delta) : super_Lden_step : (seedsSuper.l_dente_nm(s) + super_Lden_delta);
        Slist = (seedsSuper.h_si_nm(s)    - super_hsi_delta)  : super_hsi_step  : (seedsSuper.h_si_nm(s)    + super_hsi_delta);

        aStart = max(0,  aCenter - SUPER_ALPHA_HALFSPAN_SENS);
        aStop  = min(89, aCenter + SUPER_ALPHA_HALFSPAN_SENS);

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

                            super_point_idx = super_point_idx + 1;
                            if RESUME && RESUME_STAGE=="SUPER" && super_point_idx <= CKPT.done_points
                                continue;
                            end

                            setParamScalar(model, PARAM_N, n_fast(1));
                            [a1, Rp1, Rm1, TM1] = solveAndGetRplusRminus( ...
                                model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                                aStart, alpha_super_step, aStop, ...
                                M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                            [pk1, k1] = max(abs(TM1));
                            alpha1 = a1(k1);
                            tm_at_pk1 = TM1(k1);

                            if PLOT_LIVE
                                updateLivePlot('SUPER', Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), ...
                                    a1, TM1, alpha1, tm_at_pk1, Rp1, Rm1);
                            end

                            setParamScalar(model, PARAM_N, n_fast(2));
                            [a2, ~, ~, TM2] = solveAndGetRplusRminus( ...
                                model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                                aStart, alpha_super_step, aStop, ...
                                M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                            [~, k2] = max(abs(TM2));
                            alpha2 = a2(k2);

                            S_est = (alpha2 - alpha1) / (n_fast(2) - n_fast(1));

                            SuperRows = [SuperRows; ...
                                Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), ...
                                pk1, alpha1, tm_at_pk1, alpha2, S_est]; %#ok<AGROW>

                            runs_done_global = runs_done_global + runs_per_point_search;

                            stage_done = runs_done_global - stage_runs_start;
                            [frac_stage, eta_stage]   = frac_eta(stage_done, stage_total_runs, toc(stage_t0));
                            [frac_global, eta_global] = frac_eta(runs_done_global, grand_total_target, toc(t0));
                            logline(['SUPER  | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                     ' | |TM|=%.5f @ α=%.4f° | S_est=%.4f' ...
                                     ' [Stage %5.1f%% | ETA %s | Global %5.1f%% | ETA %s]\n'], ...
                                Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), ...
                                pk1, alpha1, S_est, ...
                                100*frac_stage, fmt_time_long(eta_stage), 100*frac_global, fmt_time_long(eta_global));

                            points_since_ckpt = points_since_ckpt + 1;
                            payload = struct('SuperRows',SuperRows, 'seedsSuper', seedsSuper);
                            points_since_ckpt = maybe_checkpoint( ...
                                'SUPER', CKPT_EVERY_POINTS, CKPT_FILE, PROGRESS_XLSX, ...
                                runs_done_global, points_since_ckpt, super_point_idx, payload, ...
                                Tcoarse, Tfine, [], [], []);
                        end
                    end
                end
            end
        end
    end

    Tsuper = array2table(SuperRows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm', ...
         'maxAbsTMOKE_base','alpha_peak_base_deg','TMOKE_at_peak_base', ...
         'alpha_peak_n2_deg','S_est_deg_per_RIU'});

    % Escolhas finais dentro do SUPER:
    best_tradeoff = selectTopK_tradeoff(Tsuper, 1, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');
    best_TMOKE    = selectTopK_single(Tsuper, 1, 'maxAbsTMOKE_base');      % max |TMOKE|
    best_S        = selectTopK_single_abs(Tsuper, 1, 'S_est_deg_per_RIU'); % max |S|

    fprintf('\n===== BEST (SUPERFINE) =====\n');
    fprintf('Trade-off (rankTM+rankS):\n'); disp(best_tradeoff);
    fprintf('Best TMOKE only:\n');         disp(best_TMOKE);
    fprintf('Best |S| only:\n');          disp(best_S);

    save_checkpoint(CKPT_FILE, 'VALID', runs_done_global, 0, struct( ...
        'Tcoarse', Tcoarse, 'Tfine', Tfine, 'Tsuper', Tsuper, ...
        'best_tradeoff', best_tradeoff, 'best_TMOKE', best_TMOKE, 'best_S', best_S));
    write_progress_xlsx(PROGRESS_XLSX, 'super', Tsuper);
end

%% ==================================================================
%                 VALID: Curvas densas + S_dense (3 n)
% ==================================================================
% A partir daqui, a geometria fixa é a best_tradeoff (melhor compromisso).
if isempty(best_tradeoff) && RESUME && any(strcmp(RESUME_STAGE, ["VALID","FULL"]))
    best_tradeoff = CKPT.payload.best_tradeoff;
end

Ldom_best = best_tradeoff.L_domain_nm(1);
Lden_best = best_tradeoff.l_dente_nm(1);
hsi_best  = best_tradeoff.h_si_nm(1);
hcey_best = best_tradeoff.h_ceyig_nm(1);
hau_best  = best_tradeoff.h_au_nm(1);

setParamNm(model, PARAM_LDOM, Ldom_best);
setParamNm(model, PARAM_LDEN, Lden_best);
setParamNm(model, PARAM_HSI,  hsi_best);
setParamNm(model, PARAM_HCEY, hcey_best);
setParamNm(model, PARAM_HAU,  hau_best);

alpha_peaks = zeros(size(n_valid));
TMmax_abs   = zeros(size(n_valid));
TM_all      = cell(size(n_valid));
alpha_all   = cell(size(n_valid));
Rp_all      = cell(size(n_valid));
Rm_all      = cell(size(n_valid));

for i = 1:numel(n_valid)
    setParamScalar(model, PARAM_N, n_valid(i));

    [a, Rp, Rm, TM] = solveAndGetRplusRminus( ...
        model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
        0, alpha_dense_step, 89, ...
        M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

    [tmmax, k] = max(abs(TM));
    alpha_peaks(i) = a(k);
    TMmax_abs(i)   = tmmax;

    alpha_all{i} = a;
    TM_all{i}    = TM;
    Rp_all{i}    = Rp;
    Rm_all{i}    = Rm;

    runs_done_global = runs_done_global + 2;
end

p = polyfit(n_valid, alpha_peaks, 1);
S_dense = p(1);

fprintf('\n===== VALID (best_tradeoff) =====\n');
fprintf('Best geom: Ldom=%g | Lden=%g | hsi=%g | hcey=%g | hau=%g\n', ...
    Ldom_best, Lden_best, hsi_best, hcey_best, hau_best);
fprintf('alpha_peak(n) = ['); fprintf(' %.4f', alpha_peaks); fprintf(' ]\n');
fprintf('S_dense ≈ %.6f deg/RIU (fit linear)\n', S_dense);

% Define “alpha_best / tmoke_best” como o pico no n_base (se estiver em n_valid)
[~, idxBase] = min(abs(n_valid - n_base));
alpha_best = alpha_peaks(idxBase);

% Recalcula tmoke_best no n_base (a partir dos dados já coletados)
[tmmax_base, kbase] = max(abs(TM_all{idxBase}));
tmoke_best = TM_all{idxBase}(kbase);

% Tabela densa consolidada (todas as linhas de alpha para cada n)
BestDense = table();
for i = 1:numel(n_valid)
    ncol = repmat(n_valid(i), numel(alpha_all{i}), 1);
    Ttmp = table( ...
        repmat(Ldom_best,numel(alpha_all{i}),1), repmat(Lden_best,numel(alpha_all{i}),1), repmat(hsi_best,numel(alpha_all{i}),1), ...
        repmat(hcey_best,numel(alpha_all{i}),1), repmat(hau_best,numel(alpha_all{i}),1), ...
        ncol, alpha_all{i}(:), Rp_all{i}(:), Rm_all{i}(:), TM_all{i}(:), ...
        'VariableNames', {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm', ...
                          'n','alpha_deg','Rplus','Rminus','TMOKE'});
    BestDense = [BestDense; Ttmp]; %#ok<AGROW>
end

save_checkpoint(CKPT_FILE, 'FULL', runs_done_global, 0, struct( ...
    'Tcoarse', Tcoarse, 'Tfine', Tfine, 'Tsuper', Tsuper, ...
    'best_tradeoff', best_tradeoff, 'best_TMOKE', best_TMOKE, 'best_S', best_S, ...
    'BestDense', BestDense, 'S_dense', S_dense, ...
    'alpha_best', alpha_best, 'tmoke_best', tmoke_best));
write_progress_xlsx(PROGRESS_XLSX, 'dense', BestDense);

%% ==================================================================
%                           Snapshot (optional)
% ==================================================================
% Snapshot do melhor ponto (best_tradeoff) no alpha_best e n_base
if SAVE_SNAPSHOT
    setParamScalar(model, PARAM_N, n_base);

    setAlphaMSweep(model, STUDY_TAG, ALPHA_NAME, alpha_best, 0.01, alpha_best, MSIGN_NAME, sprintf('%s %s', M_PLUS, M_MINUS));
    model.study(STUDY_TAG).run;
    refreshDerivedValues(model);

    runs_done_global = runs_done_global + 2;

    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    snap_name = fullfile(homeDir, sprintf( ...
        'snapshot_bestTradeoff_Ldom%4.0f_Lden%4.0f_hsi%3.0f_hcey%.3f_hau%.3f_n%.3f_alpha%.4f_%s.mph', ...
        Ldom_best, Lden_best, hsi_best, hcey_best, hau_best, n_base, alpha_best, timestamp));
    mphsave(model, snap_name);
    logline('Snapshot saved: %s\n', snap_name);
end

%% ==================================================================
%                     FULL: curva final (baseline n_base)
% ==================================================================
% FULL mantém a saída padrão “bonita” em um único n (n_base).
setParamScalar(model, PARAM_N, n_base);

[a_full, Rp_full, Rm_full, TM_full] = solveAndGetRplusRminus( ...
    model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
    0, alpha_full_step, 89, ...
    M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

runs_done_global = runs_done_global + 2;

if PLOT_LIVE
    updateLivePlot('FULL', Ldom_best, Lden_best, hsi_best, hcey_best, hau_best, ...
        a_full, TM_full, alpha_best, tmoke_best, Rp_full, Rm_full);
end

BestFull = table( ...
    repmat(Ldom_best,numel(a_full),1), repmat(Lden_best,numel(a_full),1), repmat(hsi_best,numel(a_full),1), ...
    repmat(hcey_best,numel(a_full),1), repmat(hau_best,numel(a_full),1), ...
    repmat(n_base,numel(a_full),1), ...
    a_full(:), Rp_full(:), Rm_full(:), TM_full(:), ...
    'VariableNames', {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm', ...
                      'n','alpha_deg','Rplus','Rminus','TMOKE'});

%% ==================================================================
%                               CSV / XLSX
% ==================================================================
if ~isempty(Tcoarse), writetable(Tcoarse, fullfile(homeDir,'tmoke_sens_5D_coarse.csv')); end
if ~isempty(Tfine),   writetable(Tfine,   fullfile(homeDir,'tmoke_sens_5D_fine.csv')); end
if ~isempty(Tsuper),  writetable(Tsuper,  fullfile(homeDir,'tmoke_sens_5D_super.csv')); end

writetable(BestDense, fullfile(homeDir,'tmoke_sens_bestTradeoff_dense_ALLn.csv'));
writetable(BestFull,  fullfile(homeDir,'tmoke_sens_bestTradeoff_full_baseline.csv'));

write_progress_xlsx(PROGRESS_XLSX, 'full', BestFull);

% Limpa checkpoint no final (execução “completa”)
if isfile(CKPT_FILE), delete(CKPT_FILE); end

%% ==================================================================
%                               PLOTS FINAIS
% ==================================================================
if MAKE_PLOTS
    % (1) TMOKE(alpha) para cada n (VALID)
    figure('Name','(1) TMOKE(alpha) for each n (VALID)','NumberTitle','off','Color','w');
    hold on; grid on;
    for i = 1:numel(n_valid)
        plot(alpha_all{i}, TM_all{i}, 'LineWidth', 1.2, 'DisplayName', sprintf('n=%.2f', n_valid(i)));
    end
    xlabel('\alpha [deg]'); ylabel('TMOKE(\alpha)');
    title(sprintf('Best trade-off | TMOKE(\\alpha) vs \\alpha (S_{dense}=%.4f deg/RIU)', S_dense));
    legend('Location','best');

    % (2) alpha_peak vs n + fit
    figure('Name','(2) alpha_peak vs n (VALID)','NumberTitle','off','Color','w');
    hold on; grid on;
    plot(n_valid, alpha_peaks, 'o', 'LineWidth', 1.5, 'DisplayName','dados');
    nfit = linspace(min(n_valid), max(n_valid), 100);
    plot(nfit, polyval(p,nfit), '-', 'LineWidth', 1.2, 'DisplayName', sprintf('fit (S=%.4f)', S_dense));
    xlabel('n'); ylabel('\alpha_{peak} [deg]');
    title('\alpha_{peak} em função de n');
    legend('Location','best');

    % (3) |TMOKE|max vs n
    figure('Name','(3) |TMOKE|max vs n (VALID)','NumberTitle','off','Color','w');
    grid on; hold on;
    plot(n_valid, TMmax_abs, 'o-', 'LineWidth', 1.5);
    xlabel('n'); ylabel('|TMOKE|_{max}');
    title('|TMOKE|_{max} vs n');

    % (4) Scatter (trade-off) no espaço (|TMOKE|, |S_est|)
    AllCand = Tcoarse;
    if ~isempty(Tfine),  AllCand = [AllCand; Tfine]; end %#ok<AGROW>
    if ~isempty(Tsuper), AllCand = [AllCand; Tsuper]; end %#ok<AGROW>

    figure('Name','(4) Trade-off map: |TMOKE| vs |S_est|','NumberTitle','off','Color','w');
    grid on; hold on;
    scatter(abs(AllCand.maxAbsTMOKE_base), abs(AllCand.S_est_deg_per_RIU), 25, 'filled');
    scatter(abs(best_tradeoff.maxAbsTMOKE_base), abs(best_tradeoff.S_est_deg_per_RIU), 120, 'filled');
    xlabel('|TMOKE|_{max} (baseline)'); ylabel('|S_{est}| [deg/RIU]');
    title('Candidatos (COARSE+FINE+SUPER) e best_tradeoff (destaque)');
end

%% ==================================================================
%                           SAVE ALL FIGURES
% ==================================================================
if SAVE_FIGS
    saveAllOpenFigures(FIG_DIR, "tmoke_sens_5d", FIG_FORMATS);
    fprintf('Figures saved in: %s\n', FIG_DIR);
end

%% ==================================================================
%                               SUMMARY
% ==================================================================
elapsed_total = toc(t0);
fprintf('\n===== SUMMARY =====\n');
fprintf('Best TRADE-OFF (SUPER):\n'); disp(best_tradeoff);
fprintf('Best TMOKE only (SUPER):\n'); disp(best_TMOKE);
fprintf('Best |S_est| only (SUPER):\n'); disp(best_S);
fprintf('VALID @ best_tradeoff: S_dense ≈ %.6f deg/RIU\n', S_dense);
fprintf('Runs done: %d | Elapsed: %s\n', runs_done_global, fmt_time_long(elapsed_total));

diary off;

%% =====================================================================
%                           LOCAL FUNCTIONS
% =====================================================================

function setParamNm(mdl, name, val_nm)
    mdl.param.set(name, sprintf('%.12g[nm]', val_nm));
end

function setParamScalar(mdl, name, val)
% Para parâmetros sem unidade (ex.: índice n)
    mdl.param.set(name, sprintf('%.12g', val));
end

function Tsel = selectTopK_tradeoff(T, K, colTM, colS)
% Seleciona TOP-K por score combinado: rank(|TM|) + rank(|S|)
% - rank 1 = melhor
    tm = abs(T.(colTM));
    ss = abs(T.(colS));

    rTM = tiedrank(-tm);
    rS  = tiedrank(-ss);

    score = rTM + rS;
    T.score_tradeoff = score; %#ok<AGROW>

    [~, ord] = sort(score, 'ascend');
    Tsel = T(ord(1:min(K,height(T))), :);
end

function Tsel = selectTopK_single(T, K, col)
% TOP-K por maior valor do col (sem abs)
    [~, ord] = sort(T.(col), 'descend');
    Tsel = T(ord(1:min(K,height(T))), :);
end

function Tsel = selectTopK_single_abs(T, K, col)
% TOP-K por maior abs(col)
    [~, ord] = sort(abs(T.(col)), 'descend');
    Tsel = T(ord(1:min(K,height(T))), :);
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
    setTwoParamSweep(mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mPlusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha1_deg, R1] = readAlphaAndRFromNamedTable(mdl, ttagPlus, Npts);

    % m = -1
    redirectAllNumericalsToTable(mdl, ttagMinus);
    clearTable(mdl, ttagMinus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mMinusStr);
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

function [alpha_deg, Rcol] = readAlphaAndRFromNamedTable(mdl, ttag, Npts)
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
        a = S.data(:,1); R = S.data(:,2);
    else
        hlow = lower(heads);
        aIdx = find(contains(hlow,'alpha'), 1, 'first');
        rIdx = find(contains(hlow,'reflectance') | contains(hlow,'total reflectance') | contains(hlow,' total r'), 1, 'first');
        if isempty(aIdx) || isempty(rIdx)
            if ~isfield(S,'data') || isempty(S.data) || size(S.data,2) < 2
                error('Alpha/Reflectance not found in %s and not enough columns.', ttag);
            end
            aIdx = 1; rIdx = 2;
        end
        a = S.data(:, aIdx); R = S.data(:, rIdx);
    end

    assert(numel(a) >= Npts && numel(R) >= Npts, ...
        'Table %s has %d rows; expected >= %d.', ttag, numel(a), Npts);

    a = a(end-Npts+1:end);
    R = R(end-Npts+1:end);

    if ~isempty(heads)
        hlow = lower(heads);
        if any(contains(hlow, '(rad)')) || any(contains(hlow, '[rad]')), a = a * 180/pi; end
    end

    alpha_deg = a;
    Rcol      = R;
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

function save_checkpoint(CKPT_FILE, stage, runs_done_global, done_points, payload)
    CKPT.stage            = char(stage);
    CKPT.runs_done_global = runs_done_global;
    CKPT.done_points      = done_points;
    CKPT.payload          = payload;
    tmp = [CKPT_FILE '.tmp'];
    save(tmp,'CKPT','-v7.3');
    movefile(tmp, CKPT_FILE, 'f');
end

function write_progress_xlsx(PROGRESS_XLSX, varargin)
    for i=1:2:numel(varargin)
        sheet = varargin{i}; T = varargin{i+1};
        if ~isempty(T)
            try
                writetable(T, PROGRESS_XLSX, 'Sheet', sheet);
            catch
                warning('write_progress_xlsx: failed sheet %s', sheet);
            end
        end
    end
end

function points_since_ckpt_new = maybe_checkpoint( ...
        stage, every_points, CKPT_FILE, PROGRESS_XLSX, ...
        runs_done_global, points_since_ckpt, done_points, payload, ...
        Tcoarse, Tfine, Tsuper, BestDense, BestFull)
    points_since_ckpt_new = points_since_ckpt;
    if points_since_ckpt >= every_points
        save_checkpoint(CKPT_FILE, stage, runs_done_global, done_points, payload);
        try
            write_progress_xlsx(PROGRESS_XLSX, ...
                'coarse', Tcoarse, 'fine', Tfine, 'super', Tsuper, ...
                'dense', BestDense, 'full', BestFull);
        catch
        end
        points_since_ckpt_new = 0;
    end
end

function updateLivePlot(stage, Ldom, Lden, hsi, hcey, hau, ...
                        alpha_deg, TMOKE_vec, alpha_peak, tmoke_peak, ...
                        Rplus, Rminus)
    persistent fig ax hTM hPeak hRp hRm inited
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
        hRp = plot(ax, nan, nan, '--', 'LineWidth', 1.0, 'DisplayName','R^+(\alpha)');
        hRm = plot(ax, nan, nan, ':',  'LineWidth', 1.0, 'DisplayName','R^-(\alpha)');
        ylabel(ax,'Reflectance');

        legend(ax,'Location','best');
        inited = true;
    end

    yyaxis(ax,'left');
    set(hTM,   'XData', alpha_deg, 'YData', TMOKE_vec);
    set(hPeak, 'XData', alpha_peak, 'YData', tmoke_peak);

    yyaxis(ax,'right');
    if nargin >= 12 && ~isempty(Rplus) && ~isempty(Rminus)
        set(hRp, 'XData', alpha_deg, 'YData', Rplus);
        set(hRm, 'XData', alpha_deg, 'YData', Rminus);
    else
        set(hRp, 'XData', nan, 'YData', nan);
        set(hRm, 'XData', nan, 'YData', nan);
    end

    title(ax, sprintf('%s | Ldom=%g | Lden=%g | hsi=%g | hcey=%.3f | hau=%.3f | |TM|=%.5f @ %.3f°', ...
        stage, Ldom, Lden, hsi, hcey, hau, abs(tmoke_peak), alpha_peak));

    drawnow('limitrate');
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
                warning('Falha ao salvar %s: %s', filePath, ME.message);
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
