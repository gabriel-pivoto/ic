% =====================================================================
% MATLAB + COMSOL — 5D Hierarchical Search with strict < 20k runs
% Params: (h_au, h_ceyig, L_domain, l_dente, h_si)
% Stages: COARSE -> FINE -> SUPERFINE -> dense validation -> final curve
% - Live plot every iteration: TMOKE (left) + R+ / R- (right)
% - Global ETA always visible: starts EST., becomes EXACT after seed clamp
% - Writes R+ in 'tblRplus' and R- in 'tblRminus' (cleared each run)
% =====================================================================

clear; clc; close all; format long; tic
import com.comsol.model.*; import com.comsol.model.util.*

% --------------------------- Paths/Project -------------------------
homeDir  = 'C:\Users\gabri\Documents\projetoIC';   % << adjust if needed
mph_file = fullfile(homeDir,'usandoMatlab.mph');
addpath(genpath(homeDir));

% ------------------------ Model Tags/Params ------------------------
STUDY_TAG   = 'std1';
PARAM_HAU   = 'h_au';
PARAM_HCEY  = 'h_ceyig';
PARAM_LDOM  = 'L_domain';
PARAM_LDEN  = 'l_dente';
PARAM_HSI   = 'h_si';
ALPHA_NAME  = 'alpha';
MSIGN_NAME  = 'm';
M_PLUS      = '1';
M_MINUS     = '-1';

% ---- RESULT TABLE TAGS (will be created if missing)
RPLUS_TABLE_TAG  = 'tblRplus';
RMINUS_TABLE_TAG = 'tblRminus';

% ------------------------ COARSE grids (exact) ---------------------
% Loop order: L_domain -> l_dente -> h_si -> h_ceyig -> h_au
Ldomain_grid_nm = 800:50:850;         % 3
ldente_grid_nm  = 500:50:600;         % 3
hsi_grid_nm     = [220, 240, 260];    % 3
hcey_grid_nm    = [100, 140];         % 2
hau_grid_nm     = 20:10:60;           % 5
% COARSE points = 3*3*3*2*5 = 270  -> runs = 270 * 2 (m=±1)

% ---------------------- Alpha ranges/steps -------------------------
alpha_coarse_rng = [0, 1.0, 89];  % full sweep for COARSE
alpha_fine_step  = 0.1;           % FINE window step
alpha_super_step = 0.01;          % SUPERFINE window step
alpha_dense_step = 0.01;          % final/validation
alpha_full_step  = 0.01;          % final curve

% ----------------------- TOP-K strategy ---------------------------
TOPK_COARSE = 1;   % refine around 1 seed
TOPK_FINE   = 1;

% --------------------- Refinement windows -------------------------
% FINE param windows
fine_hau_delta   = 2;    fine_hau_step   = 1;
fine_hcey_delta  = 5;    fine_hcey_step  = 1;
fine_Ldom_delta  = 10;   fine_Ldom_step  = 5;
fine_Lden_delta  = 10;   fine_Lden_step  = 5;
fine_hsi_delta   = 5;    fine_hsi_step   = 5;

% SUPERFINE param windows
super_hau_delta  = 1;    super_hau_step  = 0.5;
super_hcey_delta = 2;    super_hcey_step = 1;
super_Ldom_delta = 4;    super_Ldom_step = 2;
super_Lden_delta = 4;    super_Lden_step = 2;
super_hsi_delta  = 2;    super_hsi_step  = 1;

% ------------------ Alpha windows (centered at last peak) ----------
FINE_ALPHA_HALFSPAN   = 5;     % ±5° (step 0.1)
SUPER_ALPHA_HALFSPAN  = 2;     % ±2° (step 0.01)

% ----------------------- Outputs / Plots ---------------------------
SAVE_SNAPSHOT = true;
MAKE_PLOTS    = true;
PLOT_LIVE     = true;

% ------------------------- Load model ------------------------------
ModelUtil.clear;
model = mphload(mph_file);
ModelUtil.showProgress(true);

% ==================================================================
%           GLOBAL PLANNING (show estimate now, exact later)
% ==================================================================
t0 = tic;                            % global clock
runs_done_global = 0;
is_total_exact   = false;            % start with estimate
grand_total_target = NaN; %#ok<NASGU>

% ---- COARSE (exact) ----
nL = numel(Ldomain_grid_nm);
nD = numel(ldente_grid_nm);
nS = numel(hsi_grid_nm);
nG = numel(hcey_grid_nm);
nH = numel(hau_grid_nm);

coarse_total_pts  = nL*nD*nS*nG*nH;         % 270
coarse_total_runs = 2 * coarse_total_pts;   % m=±1

% ---- FINE (estimated before seed clamp) ----
fine_H = 1 + floor(2*fine_hau_delta   / fine_hau_step);
fine_G = 1 + floor(2*fine_hcey_delta  / fine_hcey_step);
fine_L = 1 + floor(2*fine_Ldom_delta  / fine_Ldom_step);
fine_D = 1 + floor(2*fine_Lden_delta  / fine_Lden_step);
fine_S = 1 + floor(2*fine_hsi_delta   / fine_hsi_step);
fine_pts_per_seed_est = fine_L*fine_D*fine_S*fine_G*fine_H;
fine_total_pts_est    = TOPK_COARSE * fine_pts_per_seed_est;
fine_total_runs_est   = 2 * fine_total_pts_est;

% ---- SUPERFINE (exact upfront) ----
super_H = 1 + floor(2*super_hau_delta   / super_hau_step);
super_G = 1 + floor(2*super_hcey_delta  / super_hcey_step);
super_L = 1 + floor(2*super_Ldom_delta  / super_Ldom_step);
super_D = 1 + floor(2*super_Lden_delta  / super_Lden_step);
super_S = 1 + floor(2*super_hsi_delta   / super_hsi_step);
super_pts_per_seed = super_L*super_D*super_S*super_G*super_H;
super_total_pts    = TOPK_FINE * super_pts_per_seed;
super_total_runs   = 2 * super_total_pts;

% ---- EXTRAS (exact) ----
extras_runs_fixed = 0;
extras_runs_fixed = extras_runs_fixed + 2;  % dense validation (m=±1)
if SAVE_SNAPSHOT
    extras_runs_fixed = extras_runs_fixed + 2; % snapshot (m=±1)
end
extras_runs_fixed = extras_runs_fixed + 2;  % final curve (m=±1)

% ---- Initial estimate for global total & header ----
grand_total_target = coarse_total_runs + fine_total_runs_est + super_total_runs + extras_runs_fixed;
is_total_exact = false;

fprintf('START\n');
fprintf('  COARSE (exact):    %d runs\n', coarse_total_runs);
fprintf('  FINE   (estimate): %d runs (TOP-%d)\n', fine_total_runs_est, TOPK_COARSE);
fprintf('  SUPER  (exact):    %d runs (TOP-%d)\n', super_total_runs, TOPK_FINE);
fprintf('  EXTRAS (exact):    %d runs\n', extras_runs_fixed);
fprintf('  TOTAL  (estimate): %d runs\n\n', grand_total_target);

% ==================================================================
%                              COARSE
% ==================================================================
fprintf('STAGE COARSE — EXACT: %d runs\n', coarse_total_runs);
stage_runs_start = runs_done_global;
stage_total_runs = coarse_total_runs;
stage_t0 = tic;

Rows = []; % [Ldom, Lden, h_si, h_cey, h_au, max|TM|, alpha*, TM*]

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

                    % Full alpha sweep (0:1:89) in COARSE
                    [a_deg, Rplus, Rminus, TM] = solveAndGetRplusRminus( ...
                        model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                        alpha_coarse_rng(1), alpha_coarse_rng(2), alpha_coarse_rng(3), ...
                        M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                    [pk, k] = max(abs(TM));
                    if k==1 || k==numel(TM)
                        logline('WARN COARSE peak at α-edge (α*=%.3f in [%.3f, %.3f])\n', a_deg(k), a_deg(1), a_deg(end));
                    end

                    if PLOT_LIVE
                        updateLivePlot('COARSE', ...
                            Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), ...
                            hcey_grid_nm(ig), hau_grid_nm(ih), ...
                            a_deg, TM, a_deg(k), TM(k), Rplus, Rminus);
                    end

                    Rows = [Rows; ...
                        Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), ...
                        hcey_grid_nm(ig), hau_grid_nm(ih), pk, a_deg(k), TM(k)]; %#ok<AGROW>

                    runs_done_global = runs_done_global + 2;

                    stage_done = runs_done_global - stage_runs_start;
                    [frac_stage, eta_stage]   = frac_eta(stage_done, stage_total_runs, toc(stage_t0));
                    [frac_global, eta_global] = frac_eta(runs_done_global, grand_total_target, toc(t0));
                    logline(['COARSE | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                             ' | max|TM|=%.5f @ α=%.2f°  [Stage %5.1f%% | ETA %s | Global %s%% | ETA %s]\n'], ...
                        Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), ...
                        hcey_grid_nm(ig), hau_grid_nm(ih), pk, a_deg(k), ...
                        100*frac_stage, fmt_time_long(eta_stage), fmt_pct(frac_global), fmt_eta_flag(eta_global, is_total_exact));
                end
            end
        end
    end
end

Tcoarse = array2table(Rows, 'VariableNames', ...
 {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm','maxAbsTMOKE','alpha_at_peak_deg','TMOKE_at_peak'});

% ---------------------- Select TOP-1 as seed ----------------------
[~, ordC] = sort(abs(Tcoarse.maxAbsTMOKE), 'descend');
kkeepC    = min(TOPK_COARSE, numel(ordC));
seeds     = Tcoarse(ordC(1:kkeepC), :);

% ==================================================================
%                        FINE PLANNING (EXACT NOW)
% ==================================================================
fine_total_pts = 0;
for s = 1:height(seeds)
    Hlist = max(min(hau_grid_nm),  seeds.h_au_nm(s)   - fine_hau_delta)  : fine_hau_step  : min(max(hau_grid_nm),  seeds.h_au_nm(s)   + fine_hau_delta);
    Glist = max(min(hcey_grid_nm), seeds.h_ceyig_nm(s)- fine_hcey_delta) : fine_hcey_step : min(max(hcey_grid_nm), seeds.h_ceyig_nm(s)+ fine_hcey_delta);
    Llist = max(min(Ldomain_grid_nm), seeds.L_domain_nm(s)- fine_Ldom_delta) : fine_Ldom_step : min(max(Ldomain_grid_nm), seeds.L_domain_nm(s)+ fine_Ldom_delta);
    Dlist = max(min(ldente_grid_nm),  seeds.l_dente_nm(s) - fine_Lden_delta) : fine_Lden_step : min(max(ldente_grid_nm),  seeds.l_dente_nm(s) + fine_Lden_delta);
    Slist = max(min(hsi_grid_nm),     seeds.h_si_nm(s)    - fine_hsi_delta)  : fine_hsi_step  : min(max(hsi_grid_nm),     seeds.h_si_nm(s)    + fine_hsi_delta);
    fine_total_pts = fine_total_pts + numel(Llist)*numel(Dlist)*numel(Slist)*numel(Glist)*numel(Hlist);
end
fine_total_runs = 2 * fine_total_pts;

% ---- GRAND TOTAL becomes EXACT from now on ----
grand_total_target = coarse_total_runs + fine_total_runs + super_total_runs + extras_runs_fixed;
is_total_exact = true;

fprintf('\nSTAGE FINE — EXACT: %d runs (TOP-%d)\n', fine_total_runs, height(seeds));
fprintf('GLOBAL TOTAL (EXACT): %d runs\n\n', grand_total_target);
if grand_total_target > 20000
    error('Planned total exceeds 20k runs (%d). Adjust deltas/grids.', grand_total_target);
end

% ==================================================================
%                                 FINE
% ==================================================================
FinePeaks = [];
stage_runs_start = runs_done_global;
stage_total_runs = fine_total_runs;
stage_t0 = tic;

for s = 1:height(seeds)
    % Center α window at the COARSE peak of this seed
    aCenter = seeds.alpha_at_peak_deg(s);

    Hlist = max(min(hau_grid_nm),  seeds.h_au_nm(s)   - fine_hau_delta)  : fine_hau_step  : min(max(hau_grid_nm),  seeds.h_au_nm(s)   + fine_hau_delta);
    Glist = max(min(hcey_grid_nm), seeds.h_ceyig_nm(s)- fine_hcey_delta) : fine_hcey_step : min(max(hcey_grid_nm), seeds.h_ceyig_nm(s)+ fine_hcey_delta);
    Llist = max(min(Ldomain_grid_nm), seeds.L_domain_nm(s)- fine_Ldom_delta) : fine_Ldom_step : min(max(Ldomain_grid_nm), seeds.L_domain_nm(s)+ fine_Ldom_delta);
    Dlist = max(min(ldente_grid_nm),  seeds.l_dente_nm(s) - fine_Lden_delta) : fine_Lden_step : min(max(ldente_grid_nm),  seeds.l_dente_nm(s) + fine_Lden_delta);
    Slist = max(min(hsi_grid_nm),     seeds.h_si_nm(s)    - fine_hsi_delta)  : fine_hsi_step  : min(max(hsi_grid_nm),     seeds.h_si_nm(s)    + fine_hsi_delta);

    % --- α window: ±5°, step 0.1 ---
    aStart = max(0,  aCenter - FINE_ALPHA_HALFSPAN);
    aStop  = min(89, aCenter + FINE_ALPHA_HALFSPAN);

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

                        [a_deg, Rplus, Rminus, TM] = solveAndGetRplusRminus( ...
                            model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                            aStart, alpha_fine_step, aStop, ...
                            M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                        [pk, k] = max(abs(TM));
                        if k==1 || k==numel(TM)
                            logline('WARN FINE peak at α-edge (α*=%.4f in [%.4f, %.4f])\n', a_deg(k), a_deg(1), a_deg(end));
                        end

                        FinePeaks = [FinePeaks; ...
                            Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), pk, a_deg(k), TM(k)]; %#ok<AGROW>

                        if PLOT_LIVE
                            updateLivePlot('FINE', ...
                                Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), ...
                                a_deg, TM, a_deg(k), TM(k), Rplus, Rminus);
                        end

                        runs_done_global = runs_done_global + 2;

                        stage_done = runs_done_global - stage_runs_start;
                        [frac_stage, eta_stage]   = frac_eta(stage_done, stage_total_runs, toc(stage_t0));
                        [frac_global, eta_global] = frac_eta(runs_done_global, grand_total_target, toc(t0));
                        logline(['FINE   | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                 ' | max|TM|=%.5f @ α=%.3f°  [Stage %5.1f%% | ETA %s | Global %5.1f%% | ETA %s]\n'], ...
                            Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), pk, a_deg(k), ...
                            100*frac_stage, fmt_time_long(eta_stage), 100*frac_global, fmt_time_long(eta_global));
                    end
                end
            end
        end
    end
end

Tfine = array2table(FinePeaks, 'VariableNames', ...
 {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm','maxAbsTMOKE','alpha_at_peak_deg','TMOKE_at_peak'});

% ------------------- TOP-1 for SUPERFINE --------------------------
[~, ordF]   = sort(abs(Tfine.maxAbsTMOKE), 'descend');
kkeepF      = min(TOPK_FINE, numel(ordF));
seedsSuper  = Tfine(ordF(1:kkeepF), :);

% ==================================================================
%                     SUPERFINE PLANNING (EXACT)
% ==================================================================
super_total_pts = 0;
for s = 1:height(seedsSuper)
    Hlist = (seedsSuper.h_au_nm(s)    - super_hau_delta)  : super_hau_step  : (seedsSuper.h_au_nm(s)    + super_hau_delta);
    Glist = (seedsSuper.h_ceyig_nm(s) - super_hcey_delta) : super_hcey_step : (seedsSuper.h_ceyig_nm(s) + super_hcey_delta);
    Llist = (seedsSuper.L_domain_nm(s)- super_Ldom_delta) : super_Ldom_step : (seedsSuper.L_domain_nm(s)+ super_Ldom_delta);
    Dlist = (seedsSuper.l_dente_nm(s) - super_Lden_delta) : super_Lden_step : (seedsSuper.l_dente_nm(s) + super_Lden_delta);
    Slist = (seedsSuper.h_si_nm(s)    - super_hsi_delta)  : super_hsi_step  : (seedsSuper.h_si_nm(s)    + super_hsi_delta);
    super_total_pts = super_total_pts + numel(Llist)*numel(Dlist)*numel(Slist)*numel(Glist)*numel(Hlist);
end
super_total_runs = 2 * super_total_pts;

fprintf('\nSTAGE SUPER — EXACT: %d runs (TOP-%d)\n', super_total_runs, kkeepF);
fprintf('GLOBAL TOTAL (EXACT): %d runs\n\n', grand_total_target);

% ==================================================================
%                             SUPERFINE
% ==================================================================
SuperPeaks = [];
stage_runs_start = runs_done_global;
stage_total_runs = super_total_runs;
stage_t0 = tic;

for s = 1:height(seedsSuper)
    % Center α window at the FINE peak
    aCenter = seedsSuper.alpha_at_peak_deg(s);

    Hlist = (seedsSuper.h_au_nm(s)    - super_hau_delta)  : super_hau_step  : (seedsSuper.h_au_nm(s)    + super_hau_delta);
    Glist = (seedsSuper.h_ceyig_nm(s) - super_hcey_delta) : super_hcey_step : (seedsSuper.h_ceyig_nm(s) + super_hcey_delta);
    Llist = (seedsSuper.L_domain_nm(s)- super_Ldom_delta) : super_Ldom_step : (seedsSuper.L_domain_nm(s)+ super_Ldom_delta);
    Dlist = (seedsSuper.l_dente_nm(s) - super_Lden_delta) : super_Lden_step : (seedsSuper.l_dente_nm(s) + super_Lden_delta);
    Slist = (seedsSuper.h_si_nm(s)    - super_hsi_delta)  : super_hsi_step  : (seedsSuper.h_si_nm(s)    + super_hsi_delta);

    % --- α window: ±2°, step 0.01 ---
    aStart = max(0,  aCenter - SUPER_ALPHA_HALFSPAN);
    aStop  = min(89, aCenter + SUPER_ALPHA_HALFSPAN);

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

                        [a_deg, Rplus, Rminus, TM] = solveAndGetRplusRminus( ...
                            model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, ...
                            aStart, alpha_super_step, aStop, ...
                            M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                        [pk, k] = max(abs(TM));
                        if k==1 || k==numel(TM)
                            logline('WARN SUPER peak at α-edge (α*=%.4f in [%.4f, %.4f])\n', a_deg(k), a_deg(1), a_deg(end));
                        end

                        SuperPeaks = [SuperPeaks; ...
                            Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), pk, a_deg(k), TM(k)]; %#ok<AGROW>

                        if PLOT_LIVE
                            updateLivePlot('SUPER', ...
                                Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), ...
                                a_deg, TM, a_deg(k), TM(k), Rplus, Rminus);
                        end

                        runs_done_global = runs_done_global + 2;

                        stage_done = runs_done_global - stage_runs_start;
                        [frac_stage, eta_stage]   = frac_eta(stage_done, stage_total_runs, toc(stage_t0));
                        [frac_global, eta_global] = frac_eta(runs_done_global, grand_total_target, toc(t0));
                        logline(['SUPER  | Ldom=%4.0f Lden=%4.0f h_si=%3.0f h_cey=%5.1f h_au=%5.1f' ...
                                 ' | max|TM|=%.5f @ α=%.4f°  [Stage %5.1f%% | ETA %s | Global %5.1f%% | ETA %s]\n'], ...
                            Llist(iL), Dlist(iD), Slist(iS), Glist(ig), Hlist(ih), pk, a_deg(k), ...
                            100*frac_stage, fmt_time_long(eta_stage), 100*frac_global, fmt_time_long(eta_global));
                    end
                end
            end
        end
    end
end

Tsuper = array2table(SuperPeaks, 'VariableNames', ...
 {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm','maxAbsTMOKE','alpha_at_peak_deg','TMOKE_at_peak'});

% ------------------ Best of SUPERFINE -----------------------------
[~, bestIdx] = max(abs(Tsuper.maxAbsTMOKE));
Ldom_best  = Tsuper.L_domain_nm(bestIdx);
Lden_best  = Tsuper.l_dente_nm(bestIdx);
hsi_best   = Tsuper.h_si_nm(bestIdx);
hcey_best  = Tsuper.h_ceyig_nm(bestIdx);
hau_best   = Tsuper.h_au_nm(bestIdx);
alpha_best = Tsuper.alpha_at_peak_deg(bestIdx);
tmoke_best = Tsuper.TMOKE_at_peak(bestIdx);

logline('SUPER  | BEST => Ldom=%4.0f nm | Lden=%4.0f nm | h_si=%3.0f nm | h_ceyig=%.3f nm | h_au=%.3f nm | |TM|=%.6f @ α≈%.4f°\n', ...
    Ldom_best, Lden_best, hsi_best, hcey_best, hau_best, abs(tmoke_best), alpha_best);

% ==================================================================
%                     Dense validation (α=0:0.01:89)
% ==================================================================
setParamNm(model, PARAM_LDOM, Ldom_best);
setParamNm(model, PARAM_LDEN, Lden_best);
setParamNm(model, PARAM_HSI,  hsi_best);
setParamNm(model, PARAM_HCEY, hcey_best);
setParamNm(model, PARAM_HAU,  hau_best);

[a_dense, Rp_dense, Rm_dense, TM_dense] = solveAndGetRplusRminus( ...
    model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, 0, alpha_dense_step, 89, M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

runs_done_global = runs_done_global + 2;

[~, kV]   = max(abs(TM_dense));
alpha_best  = a_dense(kV);
tmoke_best  = TM_dense(kV);

logline('VALID  | BEST(true) => Ldom=%4.0f nm | Lden=%4.0f nm | h_si=%3.0f nm | h_ceyig=%.3f nm | h_au=%.3f nm | |TM|=%.6f @ α=%.4f°\n', ...
    Ldom_best, Lden_best, hsi_best, hcey_best, hau_best, abs(tmoke_best), alpha_best);

if PLOT_LIVE
    updateLivePlot('VALID', Ldom_best, Lden_best, hsi_best, hcey_best, hau_best, ...
        a_dense, TM_dense, alpha_best, tmoke_best, Rp_dense, Rm_dense);
end

% ==================================================================
%                           Snapshot (optional)
% ==================================================================
if SAVE_SNAPSHOT
    setAlphaMSweep(model, STUDY_TAG, ALPHA_NAME, alpha_best, 0.01, alpha_best, MSIGN_NAME, sprintf('%s %s', M_PLUS, M_MINUS));
    model.study(STUDY_TAG).run;
    refreshDerivedValues(model);

    runs_done_global = runs_done_global + 2;

    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    snap_name = fullfile(homeDir, sprintf( ...
        'snapshot_best_Ldom%4.0f_Lden%4.0f_hsi%3.0f_hcey%.3f_hau%.3f_alpha%.4f_%s.mph', ...
        Ldom_best, Lden_best, hsi_best, hcey_best, hau_best, alpha_best, timestamp));
    mphsave(model, snap_name);
    logline('Snapshot saved: %s\n', snap_name);
end

% ==================================================================
%             Final "pretty" curve (α=0:0.01:89) at best 5D
% ==================================================================
[a_full, Rp_full, Rm_full, TM_full] = solveAndGetRplusRminus( ...
    model, STUDY_TAG, ALPHA_NAME, MSIGN_NAME, 0, alpha_full_step, 89, M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

runs_done_global = runs_done_global + 2;
logline('FULL   | Final/dense curves completed at best 5D\n');

if PLOT_LIVE
    updateLivePlot('FULL', Ldom_best, Lden_best, hsi_best, hcey_best, hau_best, ...
        a_full, TM_full, alpha_best, tmoke_best, Rp_full, Rm_full);
end

% ==================================================================
%                              CSVs
% ==================================================================
writetable(Tcoarse, fullfile(homeDir,'tmoke_grid_5D_coarse.csv'));
writetable(Tfine,   fullfile(homeDir,'tmoke_5D_fine_candidates.csv'));
writetable(Tsuper,  fullfile(homeDir,'tmoke_5D_superfine_candidates.csv'));

BestDense = table( ...
   repmat(Ldom_best,numel(a_dense),1), repmat(Lden_best,numel(a_dense),1), repmat(hsi_best,numel(a_dense),1), ...
   repmat(hcey_best,numel(a_dense),1), repmat(hau_best,numel(a_dense),1), ...
   a_dense(:), Rp_dense(:), Rm_dense(:), TM_dense(:), ...
   'VariableNames', {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm','alpha_deg','Rplus','Rminus','TMOKE'});

BestFull  = table( ...
   repmat(Ldom_best,numel(a_full),1), repmat(Lden_best,numel(a_full),1), repmat(hsi_best,numel(a_full),1), ...
   repmat(hcey_best,numel(a_full),1), repmat(hau_best,numel(a_full),1), ...
   a_full(:), Rp_full(:), Rm_full(:), TM_full(:), ...
   'VariableNames', {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm','alpha_deg','Rplus','Rminus','TMOKE'});

writetable(BestDense,fullfile(homeDir,'tmoke_best_5D_dense_alpha_0_89_step0p01.csv'));
writetable(BestFull, fullfile(homeDir,'tmoke_best_5D_full_alpha_0_89_step0p01.csv'));

% ==================================================================
%                           FINAL SUMMARY
% ==================================================================
elapsed_total = toc(t0);
[frac_global_end, ~] = frac_eta(runs_done_global, grand_total_target, elapsed_total);

fprintf('\nSummary:\n');
fprintf('  Best (validated): L_domain=%4.0f nm | l_dente=%4.0f nm | h_si=%3.0f nm | h_ceyig=%.3f nm | h_au=%.3f nm\n', ...
    Ldom_best, Lden_best, hsi_best, hcey_best, hau_best);
fprintf('  alpha* (deg)    : %.4f\n', alpha_best);
fprintf('  |TMOKE|*        : %.6f\n', abs(tmoke_best));
fprintf('  Runs            : %d / %d (Global %5.1f%%)\n', runs_done_global, grand_total_target, 100*frac_global_end);
fprintf('  Elapsed         : %s\n', fmt_time_long(elapsed_total));

% ==================================================================
%                           LOCAL FUNCTIONS
% ==================================================================
function setParamNm(mdl, name, val_nm)
    mdl.param.set(name, sprintf('%.12g[nm]', val_nm));
end

function [alpha_deg, Rplus, Rminus, TMOKE] = solveAndGetRplusRminus( ...
    mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mPlusStr, mMinusStr, ttagPlus, ttagMinus)
    % Run α sweep for m=+1 (-> ttagPlus) and m=-1 (-> ttagMinus).
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

    alpha_deg = a; Rcol = R;
end

function s = fmt_time_long(sec)
    if ~isfinite(sec) || sec < 0, s = '...'; return; end
    days = floor(sec/86400);
    rem  = sec - 86400*days;
    hrs  = floor(rem/3600);
    rem  = rem - 3600*hrs;
    mins = floor(rem/60);
    secs = floor(rem - 60*mins);
    if days > 0, s = sprintf('%dd %02d:%02d:%02d', days, hrs, mins, secs);
    else,        s = sprintf('%02d:%02d:%02d', hrs, mins, secs);
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

function updateLivePlot(stage, Ldom, Lden, hsi, hcey, hau, alpha_deg, TMOKE_vec, alpha_peak, tmoke_peak, Rplus, Rminus)
    % Live plot: TMOKE (left) and R+ / R- (right) vs alpha.
    persistent fig ax hTM hPeak hRp hRm inited
    if isempty(inited) || ~isvalidHandle(fig) || ~isvalidHandle(ax) || ...
       isempty(hTM) || ~isvalidHandle(hTM) || isempty(hPeak) || ~isvalidHandle(hPeak) || ...
       isempty(hRp) || ~isvalidHandle(hRp) || isempty(hRm) || ~isvalidHandle(hRm)

        fig = figure('Name','TMOKE vs \alpha (live)','NumberTitle','off','Tag','TMOKE_LIVE_FIG');
        ax  = axes('Parent',fig,'Tag','TMOKE_LIVE_AX'); grid(ax,'on'); hold(ax,'on');

        % Left axis: TMOKE + peak marker
        yyaxis(ax,'left');
        hTM   = plot(ax, nan, nan, '-', 'LineWidth', 1.2, 'DisplayName','TMOKE(\alpha)');
        hPeak = plot(ax, nan, nan, 'o', 'MarkerSize', 6, 'DisplayName','peak TMOKE');
        ylabel(ax,'TMOKE'); xlim(ax,[0 89]); xlabel(ax,'\alpha [deg]');

        % Right axis: R+ and R-
        yyaxis(ax,'right');
        hRp = plot(ax, nan, nan, '--', 'LineWidth', 1.0, 'DisplayName','R^+(\alpha)');
        hRm = plot(ax, nan, nan, ':',  'LineWidth', 1.0, 'DisplayName','R^-(\alpha)');
        ylabel(ax,'Reflectance');

        legend(ax,'Location','best');
        inited = true;
    end

    % Update data
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

    title(ax, sprintf('%s | L_{dom}=%g nm | l_{dente}=%g nm | h_{si}=%g nm | h_{ceyig}=%.3f nm | h_{au}=%.3f nm | |TM|_{peak}=%.5f @ %.3f°', ...
        stage, Ldom, Lden, hsi, hcey, hau, abs(tmoke_peak), alpha_peak));

    drawnow('limitrate');
end

function tf = isvalidHandle(h)
    tf = ~isempty(h) && isvalid(h);
end
