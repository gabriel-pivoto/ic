%% =====================================================================
% TMOKE 5D — Otimização de Sensibilidade vs n (deg/RIU)
%
% RESUMO DO QUE ESTE CÓDIGO FAZ
% ---------------------------------------------------------------------
% Este script automatiza uma busca (COARSE → FINE → VALID) pelos melhores
% parâmetros geométricos de uma estrutura (L_domain, l_dente, h_si, h_ceyig, h_au)
% visando maximizar a sensibilidade S (deg/RIU), onde:
%
%   1) Para cada geometria, o COMSOL roda um sweep em:
%      - alpha (ângulo de incidência, em graus)
%      - m (magnetização: +1 e -1)
%
%   2) Para cada n (índice externo), calcula-se TMOKE(alpha) = (R+ - R-) / (R+ + R-)
%
%   3) Para cada n, encontra-se alpha_peak = argmax(|TMOKE|)
%
%   4) Ajusta-se uma reta alpha_peak(n) e a sensibilidade é:
%      S = d(alpha_peak)/d(n)  [deg/RIU]
%
%   5) COARSE: varre uma grade grosseira de geometria (barata)
%   6) Seleciona TOP-K seeds por |S| para refinar
%   7) FINE: faz varredura local fina em torno do seed (mais cara)
%   8) VALID: re-roda o melhor ponto com alpha bem denso, e gera gráficos finais
%
% Além disso:
%   - Mantém CHECKPOINT em arquivo .mat para continuar execuções interrompidas
%   - Faz plots "ao vivo" por iteração:
%       (1) TMOKE(alpha) para n = [1.33 1.36 1.39]
%       (5) Heatmap TMOKE(alpha,n)
%   - No fim, gera TODOS os gráficos finais e salva em disco (PNG e PDF).
%
%
% FLUXOGRAMA (alto nível)
% ---------------------------------------------------------------------
%  [Start]
%     |
%     v
%  (Load model .mph)  ----> (Load checkpoint? resume stage?)
%     |
%     v
%  [COARSE stage]
%     |  loop em (Ldom, Lden, hsi, hcey, hau)
%     |    -> para cada geometria:
%     |        - para cada n:
%     |            - roda COMSOL (m=+1) e (m=-1)
%     |            - calcula TMOKE(alpha)
%     |            - encontra alpha_peak(n)
%     |        - fit linear alpha_peak(n) -> S
%     |        - plots live (TMOKE e heatmap)
%     |        - salva checkpoint periodicamente
%     v
%  (Select TOP-K seeds por |S|)
%     |
%     v
%  [FINE stage]
%     |  loop na vizinhança fina do seed
%     |    -> repete avaliação de S + plots + checkpoint
%     v
%  (Pick best por |S| global)
%     |
%     v
%  [VALID stage]
%     |  reavalia melhor ponto com alpha denso
%     |  gera gráficos finais + CSVs
%     v
%  [SAVE FIGURES]  -> salva TODAS as figuras abertas em PNG/PDF
%     |
%     v
%   [End]
% =====================================================================

clear; clc; close all; format long;
import com.comsol.model.*
import com.comsol.model.util.*

%% --------------------------- Paths/Model ------------------------------
% Pasta base do projeto (onde está o .mph e onde você quer salvar outputs)
homeDir  = 'C:\Users\gabri\Documents\projetoIC';   % ajuste se precisar
mph_file = fullfile(homeDir,'usandoMatlab.mph');
addpath(genpath(homeDir));

% Pasta e arquivo de checkpoint:
% - checkpoint salva progresso da busca e evita perder tempo se cair energia
CKPT_DIR  = fullfile(homeDir,'checkpoints');
if ~exist(CKPT_DIR,'dir'), mkdir(CKPT_DIR); end
CKPT_FILE = fullfile(CKPT_DIR,'tmoke_sens_5d_checkpoint.mat');

%% ---------------------- Saída de gráficos (SAVE) ----------------------
% Aqui definimos a pasta onde todos os gráficos serão salvos ao final.
% Dica: usar timestamp evita sobrescrever resultados anteriores.
SAVE_FIGS = true;
RUN_STAMP = datestr(now,'yyyy-mm-dd_HHMMSS');
FIG_DIR   = fullfile(homeDir,'plots', ['tmoke_sens_5d_' RUN_STAMP]);
if SAVE_FIGS && ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

% Formatos recomendados:
% - PNG: bom pra visualização rápida, relatórios, etc.
% - PDF (vetorial): ótimo pra artigos, zoom sem pixelar
FIG_FORMATS = {'png','pdf'};

%% ------------------------- Tags/params COMSOL -------------------------
% Importante: esses nomes precisam bater com os parâmetros e resultados
% definidos no seu arquivo COMSOL (.mph).
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

% Tabelas onde os Derived Values do COMSOL serão gravados por magnetização.
% Você precisa ter "Derived Values" (numerical) no COMSOL apontando para essas tabelas
% (ou este script redireciona automaticamente).
RPLUS_TABLE_TAG  = 'tblRplus';
RMINUS_TABLE_TAG = 'tblRminus';

%% ------------------------ Sensibilidade vs n --------------------------
% "na" define os índices externos usados pra medir sensibilidade.
% Quanto mais pontos, mais robusto o ajuste linear, porém mais caro.
na = [1.30 1.33 1.36];
% Primeiro roda n = 1.33; depois rastreia o mesmo pico em 1.30 e 1.36
% dentro de uma janela angular de +/-20 deg, respeitando alpha_start/alpha_stop.
tracking_reference_n     = 1.33;
tracking_half_window_deg = 20;

% Sweep de alpha (incidência)
alpha_start        = 0;
alpha_stop         = 89;

% ======= TESTE RÁPIDO (alpha mais grosso) =======
% Observação:
% - COARSE usa passo grosso pra reduzir custo
% - FINE usa passo intermediário pra refinar
% - VALID usa passo muito fino para consolidar (mais caro porém mais confiável)
alpha_step_coarse  = 5.0;
alpha_step_fine    = 2.0;
alpha_step_dense   = 0.5;   % (VALID)

%% ------------------------ Grade COARSE (TESTE RÁPIDO) -----------------
% Grade pequena para validar pipeline / estabilidade / tabelas COMSOL.
Ldomain_grid_nm = [800 850];     % 2
ldente_grid_nm  = [500 600];     % 2
hsi_grid_nm     = [220 260];     % 2
hcey_grid_nm    = [100 140];     % 2
hau_grid_nm     = [20 40];       % 2
% Total COARSE = 2*2*2*2*2 = 32 pontos

%% ------------------------ Janela FINE (TESTE RÁPIDO) ------------------
% Define uma vizinhança em torno do seed do COARSE, com passos pequenos.
% No "modo teste", deixamos poucos pontos para rodar rápido.
fine_hau_delta   = 2;    fine_hau_step   = 2;
fine_hcey_delta  = 5;    fine_hcey_step  = 5;
fine_Ldom_delta  = 10;   fine_Ldom_step  = 10;
fine_Lden_delta  = 10;   fine_Lden_step  = 10;
fine_hsi_delta   = 10;   fine_hsi_step   = 10;

TOPK_COARSE = 1;   % Seeds do COARSE que entram no FINE (no teste, 1)

%% ---------------------------- Checkpoint -----------------------------
% O checkpoint guarda:
% - stage atual (COARSE/FINE/VALID)
% - quantos pontos já foram processados naquele stage
% - tabelas Tcoarse/Tfine acumuladas
% - seeds usados no FINE
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
CKPT_EVERY_POINTS  = 5;  % salva a cada N pontos (teste)

%% --------------------------- Carrega modelo --------------------------
% Reinicia COMSOL e carrega o modelo.
% Dica: se você estiver iterando muito, manter o mph aberto e só variar params é bom.
ModelUtil.clear;
model = mphload(mph_file);
ModelUtil.showProgress(true);

t0 = tic;

%% =====================================================================
%                 ESTIMATIVA DE NÚMERO DE ITERAÇÕES
% =====================================================================
% Esse bloco só imprime uma noção do custo total (pontos e runs COMSOL).
nL = numel(Ldomain_grid_nm);
nD = numel(ldente_grid_nm);
nS = numel(hsi_grid_nm);
nG = numel(hcey_grid_nm);
nH = numel(hau_grid_nm);

coarse_total_pts  = nL*nD*nS*nG*nH;
coarse_total_runs = coarse_total_pts * numel(na) * 2; % 2 por m=+1/-1

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

%% =====================================================================
%                             COARSE
% =====================================================================
% Objetivo do COARSE:
% - Fazer uma varredura ampla e barata
% - Encontrar "regiões promissoras" (seeds) para refinamento
coarse_done_points = 0;

if RESUME && any(RES_STAGE == ["FINE","VALID"])
    fprintf('Pulando COARSE (já feito).\n');
else
    fprintf('\n===== STAGE COARSE =====\n');
    % Rows guarda uma linha por geometria:
    % [Ldom, Lden, h_si, h_cey, h_au, S, sign(S)]
    Rows = [];

    if RESUME && RES_STAGE == "COARSE"
        runs_done_global   = CKPT.runs_done_global;
        coarse_done_points = CKPT.done_points;
        if isfield(CKPT,'Rows'), Rows = CKPT.Rows; end
        if isfield(CKPT,'iter_global'), iter_global = CKPT.iter_global; end
        fprintf('>>> RESUME COARSE: já tinha %d pontos.\n', coarse_done_points);
    end

    % flat_idx permite retomar exatamente onde parou, contando combinações
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
                            % Já fizemos esse ponto antes (resume)
                            continue;
                        end

                        % ---- Avalia sensibilidade S e retorna TMOKE(α) ----
                        % Retornos:
                        %   Sval: slope do fit alpha_peak(n) [deg/RIU]
                        %   alpha_peaks: alpha_peak para cada n
                        %   alpha_grid: grid de alpha usado
                        %   TM_mat: TMOKE para cada n ao longo de alpha
                        [Sval, alpha_peaks, alpha_grid, TM_mat] = evaluate_tracked_sensitivity_with_tmoke( ...
                            model, STUDY_TAG, PARAM_N, na, tracking_reference_n, tracking_half_window_deg, ...
                            ALPHA_NAME, MSIGN_NAME, ...
                            alpha_start, alpha_step_coarse, alpha_stop, ...
                            M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

                        % Armazena resultado dessa geometria
                        Rows = [Rows; ...
                            Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), ...
                            hcey_grid_nm(ig), hau_grid_nm(ih), Sval, sign(Sval)]; %#ok<AGROW>

                        % Atualiza contadores globais (para debug e custo)
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

                        % ---- Plot "ao vivo" (útil para entender se TMOKE está ok) ----
                        geomLabel = sprintf('COARSE it=%d | Ldom=%.0f Lden=%.0f hsi=%.0f hcey=%.1f hau=%.1f', ...
                            iter_global, Ldomain_grid_nm(iL), ldente_grid_nm(iD), hsi_grid_nm(iS), hcey_grid_nm(ig), hau_grid_nm(ih));
                        livePlot_TMOKE_and_Heatmap(alpha_grid, na, TM_mat, geomLabel);

                        % ---- Plot online da evolução de S ----
                        updateSensitivityPlot(iter_global, Sval, 'COARSE');

                        % ---- Checkpoint periódico ----
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

    % Converte matriz Rows -> table (mais fácil de ordenar e salvar)
    Tcoarse = array2table(Rows, 'VariableNames', ...
        {'L_domain_nm','l_dente_nm','h_si_nm','h_ceyig_nm','h_au_nm','S_deg_per_RIU','S_sign'});

    % Salva checkpoint final do COARSE
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

%% =====================================================================
%                  Seleção de seeds (TOP-K) para FINE
% =====================================================================
% Critério: maior |S| (magnitude) → mais promissor (independe do sinal)
if isempty(Tcoarse)
    error('Tcoarse vazio (nada calculado em COARSE).');
end

[~, ordC] = sort(abs(Tcoarse.S_deg_per_RIU), 'descend');
kkeepC    = min(TOPK_COARSE, numel(ordC));
seeds     = Tcoarse(ordC(1:kkeepC), :);

fprintf('\nSeeds para FINE (TOP-%d por |S|):\n', kkeepC);
disp(seeds);

%% =====================================================================
%                             FINE
% =====================================================================
% Objetivo do FINE:
% - Refinar ao redor dos seeds, varrendo uma vizinhança local
% - Achar um máximo melhor sem varrer o espaço inteiro
fine_done_points = 0;

if RESUME && any(RES_STAGE == ["VALID"])
    fprintf('Pulando FINE (já feito).\n');
else
    fprintf('\n===== STAGE FINE =====\n');
    FineRows = [];

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
        % Monta listas locais ao redor do seed
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

                            [Sval, alpha_peaks, alpha_grid, TM_mat] = evaluate_tracked_sensitivity_with_tmoke( ...
                                model, STUDY_TAG, PARAM_N, na, tracking_reference_n, tracking_half_window_deg, ...
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

%% =====================================================================
%                       Escolhe melhor (global) por |S|
% =====================================================================
% Junta candidatos do COARSE e FINE e pega o melhor por |S|.
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

% Atualiza checkpoint para indicar que agora entraremos no VALID
CKPT.stage            = 'VALID';
CKPT.done_points      = 0;
CKPT.runs_done_global = runs_done_global;
CKPT.Tcoarse          = Tcoarse;
CKPT.Tfine            = Tfine;
CKPT.bestStruct       = bestStruct;
CKPT.iter_global      = iter_global;
save(CKPT_FILE,'CKPT','-v7.3');

%% =====================================================================
%                 VALID: curva densa + guarda dados p/ gráficos
% =====================================================================
% VALID é a confirmação final: alpha mais denso → alpha_peak mais confiável.
fprintf('\n===== STAGE VALID (reavalia melhor ponto) =====\n');

setParamNm(model, PARAM_LDOM, bestStruct.Ldom_best);
setParamNm(model, PARAM_LDEN, bestStruct.Lden_best);
setParamNm(model, PARAM_HSI,  bestStruct.hsi_best);
setParamNm(model, PARAM_HCEY, bestStruct.hcey_best);
setParamNm(model, PARAM_HAU,  bestStruct.hau_best);

[S_dense, alpha_peaks, alpha_grid_valid, TM_mat_valid, tracked_tmoke_abs, ...
    alpha_all, TM_all, Rplus_all, Rminus_all] = evaluate_tracked_sensitivity_with_tmoke( ...
        model, STUDY_TAG, PARAM_N, na, tracking_reference_n, tracking_half_window_deg, ...
        ALPHA_NAME, MSIGN_NAME, ...
        alpha_start, alpha_step_dense, alpha_stop, ...
        M_PLUS, M_MINUS, RPLUS_TABLE_TAG, RMINUS_TABLE_TAG);

runs_done_global = runs_done_global + numel(na)*2;
p = polyfit(na, alpha_peaks, 1);

fprintf('alpha_peak (deg) = ['); fprintf(' %.4f', alpha_peaks); fprintf(' ]\n');
fprintf('Sensibilidade (VALID) S ≈ %.6f deg/RIU\n', S_dense);

%% =====================================================================
%                 GRÁFICOS FINAIS (TODOS)
% =====================================================================

% (1) TMOKE(α) por n (melhor ponto)
figure('Name','(1) TMOKE vs alpha (VALID)','NumberTitle','off','Color','w');
hold on; grid on;
for i = 1:numel(na)
    plot(alpha_all{i}, TM_all{i}, 'LineWidth',1.2, ...
        'DisplayName', sprintf('n = %.2f', na(i)));
end
xlabel('\alpha [deg]');
ylabel('TMOKE(\alpha)');
title('TMOKE(\alpha) no melhor ponto (VALID)');
legend('Location','best');

% (2) alpha_peak vs n (com fit)
figure('Name','(2) alpha_{peak} vs n (VALID)','NumberTitle','off','Color','w');
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

% (3) |TMOKE| no pico rastreado vs n
figure('Name','(3) |TMOKE| no pico rastreado vs n (VALID)','NumberTitle','off','Color','w');
grid on; hold on;
plot(na, tracked_tmoke_abs, 'o-','LineWidth',1.5);
xlabel('Índice de refração n');
ylabel('|TMOKE| no pico rastreado');
title('Magnitude de |TMOKE| no pico rastreado vs n (VALID)');

% (5) Heatmap TMOKE(α,n) (VALID)
figure('Name','(5) Heatmap TMOKE(alpha,n) (VALID)','NumberTitle','off','Color','w');
imagesc(alpha_grid_valid, na, TM_mat_valid);
axis xy; colorbar;
xlabel('\alpha [deg]');
ylabel('n');
title('Heatmap TMOKE(\alpha,n) no melhor ponto (VALID)');

% (6) |S| histogram (COARSE e FINE)
figure('Name','(6) Histograma |S| (COARSE/FINE)','NumberTitle','off','Color','w');
hold on; grid on;
histogram(abs(Tcoarse.S_deg_per_RIU), 'DisplayName','COARSE |S|');
if ~isempty(Tfine)
    histogram(abs(Tfine.S_deg_per_RIU), 'DisplayName','FINE |S|');
end
xlabel('|S| [deg/RIU]');
ylabel('contagem');
title('Distribuição de |S|');
legend('Location','best');

% (7) Scatter |S| vs cada parâmetro (COARSE+FINE)
All = AllCand;
figure('Name','(7) |S| vs parâmetros (COARSE+FINE)','NumberTitle','off','Color','w');
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

%% --------------------------- CSVs ------------------------------------
% Exporta resultados numéricos para análise fora do MATLAB
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

%% =====================================================================
%                      SALVAR TODOS OS GRÁFICOS (FINAL)
% =====================================================================
% Este bloco salva todas as figuras abertas (inclui "live" e "S vs iteração")
% no final do script.
% Se você quiser salvar SOMENTE os finais, você pode:
%   - fechar os "live" antes, ou
%   - adaptar saveAllOpenFigures(...) para filtrar por Name.
if SAVE_FIGS
    saveAllOpenFigures(FIG_DIR, "tmoke_sens_5d", FIG_FORMATS);
end

%% =====================================================================
%                           FUNÇÕES LOCAIS
% =====================================================================

function setParamNm(mdl, name, val_nm)
    % Helper: garante unidade [nm] e formatação consistente.
    mdl.param.set(name, sprintf('%.12g[nm]', val_nm));
end

function [Sval, alpha_peaks, alpha_grid, TM_mat] = evaluate_sensitivity_with_tmoke( ...
        mdl, studyTag, PARAM_N, na, ...
        alphaName, mName, aStartDeg, aStepDeg, aStopDeg, ...
        mPlusStr, mMinusStr, ttagPlus, ttagMinus)
% Calcula sensibilidade S e também devolve TMOKE(alpha) para cada n.
%
% OUTPUTS:
%   alpha_grid: vetor de alpha (deg)
%   TM_mat:     matriz [numel(na) x numel(alpha_grid)] com TMOKE
%   alpha_peaks: alpha_peak para cada n
%   Sval: slope do ajuste linear alpha_peak(n) [deg/RIU]
%
% Nota prática:
% - Se TMOKE vier vazio, normalmente indica que:
%   (a) a tabela do COMSOL não foi preenchida (Derived Values não rodou),
%   (b) a coluna certa de reflectance não foi encontrada,
%   (c) ou o Study não executou como esperado.
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

        % Primeira iteração define o grid base (assumimos que é o mesmo para todos os n)
        if isempty(alpha_grid)
            alpha_grid = alpha_deg(:).';
            TM_mat = zeros(numel(na), numel(alpha_grid));
        end

        TM_mat(i,:) = TM(:).';

        % alpha_peak é o ângulo onde |TMOKE| atinge máximo
        [~, k] = max(abs(TM));
        alpha_peaks(i) = alpha_deg(k);
    end

    % Fit linear para obter S = d(alpha_peak)/d(n)
    p = polyfit(na, alpha_peaks, 1);
    Sval = p(1);   % slope [deg/RIU]
end

function [Sval, alpha_peaks, alpha_grid, TM_mat, tracked_tmoke_abs, ...
          alpha_all, TM_all, Rplus_all, Rminus_all] = evaluate_tracked_sensitivity_with_tmoke( ...
        mdl, studyTag, PARAM_N, na, trackingReferenceN, trackingHalfWindowDeg, ...
        alphaName, mName, aStartDeg, aStepDeg, aStopDeg, ...
        mPlusStr, mMinusStr, ttagPlus, ttagMinus)
% Avalia a sensibilidade rastreando sempre o mesmo pico de TMOKE.
% Ordem de execuÃ§Ã£o:
%   1) n = trackingReferenceN  -> escolhe o pico global de referÃªncia
%   2) n menores              -> procura o pico correspondente em +/- janela
%   3) n maiores              -> idem
    alpha_peaks = zeros(size(na));
    alpha_grid = [];
    TM_mat = [];
    tracked_tmoke_abs = zeros(size(na));
    alpha_all  = cell(size(na));
    TM_all     = cell(size(na));
    Rplus_all  = cell(size(na));
    Rminus_all = cell(size(na));

    eval_order = build_tracking_order(na, trackingReferenceN);
    reference_peak_alpha_deg = [];
    reference_peak_sign = 0;

    for order_pos = 1:numel(eval_order)
        i = eval_order(order_pos);
        n_val = na(i);
        mdl.param.set(PARAM_N, sprintf('%.12g', n_val));

        [alpha_deg, Rplus, Rminus, TM] = solveAndGetRplusRminus( ...
            mdl, studyTag, alphaName, mName, ...
            aStartDeg, aStepDeg, aStopDeg, ...
            mPlusStr, mMinusStr, ttagPlus, ttagMinus);

        if isempty(TM)
            error('TMOKE vazio na avaliaÃ§Ã£o de sensibilidade rastreada.');
        end

        if isempty(alpha_grid)
            alpha_grid = alpha_deg(:).';
            TM_mat = zeros(numel(na), numel(alpha_grid));
        else
            assert(numel(alpha_deg) == numel(alpha_grid) && max(abs(alpha_deg(:).' - alpha_grid)) < 1e-9, ...
                'Grids de alpha diferem entre os valores de n.');
        end

        TM_mat(i,:) = TM(:).';
        alpha_all{i}  = alpha_deg;
        TM_all{i}     = TM;
        Rplus_all{i}  = Rplus;
        Rminus_all{i} = Rminus;

        if isempty(reference_peak_alpha_deg)
            [~, k] = max(abs(TM));
            reference_peak_alpha_deg = alpha_deg(k);
            reference_peak_sign = sign(TM(k));
        else
            k = select_peak_in_tracking_window(alpha_deg, TM, reference_peak_alpha_deg, ...
                reference_peak_sign, trackingHalfWindowDeg, aStartDeg, aStopDeg);
        end

        alpha_peaks(i) = alpha_deg(k);
        tracked_tmoke_abs(i) = abs(TM(k));
    end

    p = polyfit(na, alpha_peaks, 1);
    Sval = p(1);
end

function eval_order = build_tracking_order(na, trackingReferenceN)
    tol = 1e-9;
    idx_ref = find(abs(na - trackingReferenceN) < tol, 1, 'first');
    if isempty(idx_ref)
        error('O valor de referÃªncia n = %.4f nÃ£o estÃ¡ presente em na.', trackingReferenceN);
    end

    idx_lower = find(na < trackingReferenceN - tol);
    idx_upper = find(na > trackingReferenceN + tol);

    [~, ord_lower] = sort(na(idx_lower), 'descend');
    [~, ord_upper] = sort(na(idx_upper), 'ascend');

    eval_order = [idx_ref, idx_lower(ord_lower), idx_upper(ord_upper)];
end

function peak_idx = select_peak_in_tracking_window(alpha_deg, TM, referencePeakAlphaDeg, ...
        referencePeakSign, trackingHalfWindowDeg, aStartDeg, aStopDeg)
    window_start_deg = max(aStartDeg, referencePeakAlphaDeg - trackingHalfWindowDeg);
    window_stop_deg  = min(aStopDeg,  referencePeakAlphaDeg + trackingHalfWindowDeg);

    in_window = alpha_deg >= window_start_deg - 1e-9 & alpha_deg <= window_stop_deg + 1e-9;
    if ~any(in_window)
        error('A janela de rastreio [%.3f, %.3f] nÃ£o contÃ©m pontos de alpha.', ...
            window_start_deg, window_stop_deg);
    end

    window_indices = find(in_window);
    candidate_indices = window_indices;
    if referencePeakSign ~= 0
        same_sign_mask = sign(TM(window_indices)) == referencePeakSign;
        if any(same_sign_mask)
            candidate_indices = window_indices(same_sign_mask);
        end
    end

    [~, local_idx] = max(abs(TM(candidate_indices)));
    peak_idx = candidate_indices(local_idx);
end

function [alpha_deg, Rplus, Rminus, TMOKE] = solveAndGetRplusRminus( ...
    mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, ...
    mPlusStr, mMinusStr, ttagPlus, ttagMinus)

    % Garante que as tabelas existem
    ensureTable(mdl, ttagPlus);
    ensureTable(mdl, ttagMinus);

    % Número esperado de pontos no sweep de alpha
    Npts = 1 + floor((aStopDeg - aStartDeg)/aStepDeg + 1e-9);

    % -------------------- m = +1 --------------------
    % Redireciona todo Derived Values para ttagPlus, limpa tabela e roda
    redirectAllNumericalsToTable(mdl, ttagPlus);
    clearTable(mdl, ttagPlus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, ...
                     aStartDeg, aStepDeg, aStopDeg, mPlusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha1_deg, R1] = readAlphaAndRFromNamedTable(mdl, ttagPlus, Npts);

    % -------------------- m = -1 --------------------
    redirectAllNumericalsToTable(mdl, ttagMinus);
    clearTable(mdl, ttagMinus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, ...
                     aStartDeg, aStepDeg, aStopDeg, mMinusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha2_deg, R2] = readAlphaAndRFromNamedTable(mdl, ttagMinus, Npts);

    % Segurança: garante que grids batem
    assert(numel(alpha1_deg)==Npts && numel(alpha2_deg)==Npts, 'Unexpected α sweep length.');
    assert(max(abs(alpha1_deg - alpha2_deg)) < 1e-9, 'α grids differ between m=+1 and m=-1.');

    alpha_deg = alpha1_deg;
    Rplus  = R1;
    Rminus = R2;

    % TMOKE = (R+ - R-) / (R+ + R-)
    % Proteção contra divisão por ~0
    denom = Rplus + Rminus;
    denom(abs(denom) < 1e-9) = 1e-9;
    TMOKE = (Rplus - Rminus) ./ denom;
end

function setTwoParamSweep(mdl, studyTag, alphaName, mName, ...
                          aStartDeg, aStepDeg, aStopDeg, mStr)
    % Ajusta o sweep paramétrico do Study para (alpha, m)
    % Importante: o COMSOL espera unidades nos strings do range.
    ptag = getParametricTag(mdl, studyTag);
    mdl.study(studyTag).feature(ptag).set('pname', {alphaName, mName});
    mdl.study(studyTag).feature(ptag).set('punit', {'deg','1'});
    alist = sprintf('range(%.12g[deg], %.12g[deg], %.12g[deg])', ...
        aStartDeg, aStepDeg, aStopDeg);
    mdl.study(studyTag).feature(ptag).set('plistarr', {alist, mStr});
end

function ptag = getParametricTag(mdl, studyTag)
    % Descobre o tag do "Parametric Sweep" dentro do Study.
    % Se não achar, assume 'param' (fallback comum).
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
    % Força atualização dos "Derived Values" (numerical nodes) para preencher tabelas.
    try
        ntags = cell(mdl.result().numerical().tags());
        for k = 1:numel(ntags)
            mdl.result().numerical(ntags{k}).setResult;
        end
    catch
    end
end

function ensureTable(mdl, ttag)
    % Garante que a tabela existe no Result node do COMSOL.
    try
        mdl.result().table(ttag);
    catch
        try mdl.result().table().create(ttag, 'Table'); catch, end
    end
end

function clearTable(mdl, ttag)
    % Limpa dados antigos (evita acumular linhas e bagunçar leitura).
    try mdl.result().table(ttag).clearTableData; catch, end
end

function redirectAllNumericalsToTable(mdl, ttag)
    % Direciona todos os "numerical" nodes para gravarem na tabela ttag.
    % Isso é útil quando o COMSOL tem vários Derived Values e você quer controlar output.
    try
        ntags = cell(mdl.result().numerical().tags());
        for k = 1:numel(ntags)
            try mdl.result().numerical(ntags{k}).set('table', ttag); catch, end
        end
    catch
    end
end

function [alpha_deg, Rcol] = readAlphaAndRFromNamedTable(mdl, ttag, Npts)
    % Lê da tabela do COMSOL (mphtable) e tenta localizar colunas:
    % - alpha (ou algo contendo 'alpha')
    % - reflectance (ou variações no header)
    % Se não achar header, pega colunas 1 e 2 como fallback.

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

    % Pega apenas os últimos Npts (caso a tabela tenha acumulado runs)
    assert(numel(a) >= Npts && numel(R) >= Npts, ...
        'Table %s has %d rows; expected >= %d.', ttag, numel(a), Npts);
    a = a(end-Npts+1:end);
    R = R(end-Npts+1:end);

    % Conversão rad->deg se necessário (se tabela veio em rad)
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
        fig = figure('Name','Sensibilidade vs Iteração','NumberTitle','off','Color','w');
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
        fig = figure('Name','LIVE: (1) TMOKE(alpha) + (5) Heatmap','NumberTitle','off','Color','w');
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

%% ======================= FUNÇÕES DE SALVAMENTO =======================
function saveAllOpenFigures(outDir, prefix, formats)
% Salva TODAS as figuras abertas no MATLAB.
% - outDir: pasta de saída
% - prefix: prefixo no nome do arquivo
% - formats: cell array, ex: {'png','pdf'}
%
% Estratégia:
% - usa o Name da figura (quando existe) para gerar nome amigável
% - sanitiza caracteres inválidos em arquivo
% - exportgraphics (preferível) com fallback para saveas

    if ~exist(outDir,'dir'), mkdir(outDir); end

    figs = findall(0,'Type','figure');

    % Opcional: ordenar por número para ficar previsível
    [~, idx] = sort([figs.Number]);
    figs = figs(idx);

    for i = 1:numel(figs)
        fig = figs(i);

        figName = string(get(fig,'Name'));
        if strlength(figName)==0
            figName = "Figure_" + string(fig.Number);
        end

        base = string(prefix) + "__" + string(sanitizeFilename(figName));

        % Algumas versões do MATLAB podem exportar melhor com renderer padrão
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
                    case "eps"
                        exportgraphics(fig, filePath, 'ContentType','vector');
                    otherwise
                        saveas(fig, filePath);
                end
            catch ME
                warning('Falha ao salvar %s: %s', filePath, ME.message);
            end
        end
    end

    fprintf('>> %d figuras salvas em: %s\n', numel(figs), outDir);
end

function s = sanitizeFilename(str)
% Remove caracteres problemáticos para nome de arquivo (Windows-friendly)
    s = regexprep(char(str), '[^\w\d\-]+', '_'); % troca não-alfanum por "_"
    s = regexprep(s, '_+', '_');                % comprime "__" -> "_"
    s = regexprep(s, '^_|_$', '');              % remove "_" no começo/fim
end
