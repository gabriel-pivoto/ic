%% =====================================================================
% Report the TOP-3 geometries from a (possibly unfinished) run
% ---------------------------------------------------------------------
% Reads the periodic checkpoint written by senseAndTmoke_semCeyig.m,
% rebuilds the results table for the most advanced stage available, and
% prints THREE rankings to the terminal:
%
%   1) TOP-3 by |TMOKE|          (best TMOKE)
%   2) TOP-3 by |S|              (best angular sensitivity)
%   3) TOP-3 by trade-off        (rank|TMOKE| + rank|S|, same as main run)
%
% It does NOT call COMSOL. It only reads the checkpoint .mat and applies
% the SAME selection logic as the main script.
% =====================================================================
clear; clc; format long;

%% --------------------------- Config --------------------------------
projectRootDir     = 'D:\Gabriel Pivoto\projetoIC';
checkpointFilePath = fullfile(projectRootDir, 'checkpoints', ...
                              'tmoke_sens_4d_sem_ceyig_checkpoint.mat');

topN = 3;   % how many to show per ranking

columnNames = {'L_domain_nm','l_dente_nm','h_si_nm','h_au_nm', ...
               'maxAbsTMOKE_base','alpha_peak_base_deg','TMOKE_at_peak_base', ...
               'alpha_peak_high_deg','S_est_deg_per_RIU'};

%% ----------------------- Load checkpoint ---------------------------
if ~isfile(checkpointFilePath)
    error('Checkpoint not found: %s', checkpointFilePath);
end
S = load(checkpointFilePath, 'checkpointData');
cp = S.checkpointData;
payload = cp.payload;

fprintf('Checkpoint loaded.\n');
fprintf('  stage saved   : %s\n', string(cp.stage));
fprintf('  done_points   : %d\n', cp.done_points);
fprintf('  runs (global) : %d\n\n', cp.runsCompletedGlobal);

%% --------- Pick the most advanced results table available ----------
% Priority: SUPER (partial superRows or finished superResultsTable) ->
%           FINE  -> COARSE. We want the deepest refinement reached.
resultsTable = [];
stageUsed = '';

if isfield(payload,'superResultsTable') && ~isempty(payload.superResultsTable)
    resultsTable = payload.superResultsTable;       stageUsed = 'SUPER (finished)';
elseif isfield(payload,'superRows') && ~isempty(payload.superRows)
    resultsTable = array2table(payload.superRows, 'VariableNames', columnNames);
    stageUsed = 'SUPER (partial)';
elseif isfield(payload,'fineResultsTable') && ~isempty(payload.fineResultsTable)
    resultsTable = payload.fineResultsTable;        stageUsed = 'FINE (finished)';
elseif isfield(payload,'fineRows') && ~isempty(payload.fineRows)
    resultsTable = array2table(payload.fineRows, 'VariableNames', columnNames);
    stageUsed = 'FINE (partial)';
elseif isfield(payload,'coarseResultsTable') && ~isempty(payload.coarseResultsTable)
    resultsTable = payload.coarseResultsTable;      stageUsed = 'COARSE (finished)';
elseif isfield(payload,'coarseRows') && ~isempty(payload.coarseRows)
    resultsTable = array2table(payload.coarseRows, 'VariableNames', columnNames);
    stageUsed = 'COARSE (partial)';
else
    error('No usable results found in checkpoint payload.');
end

if ~istable(resultsTable)
    resultsTable = array2table(resultsTable, 'VariableNames', columnNames);
end

fprintf('Using stage: %s | %d points evaluated so far.\n', stageUsed, height(resultsTable));

%% ------------------------- Rankings --------------------------------
k = min(topN, height(resultsTable));

topTmoke    = selectTopK_single    (resultsTable, k, 'maxAbsTMOKE_base');
topSens     = selectTopK_single_abs(resultsTable, k, 'S_est_deg_per_RIU');
topTradeoff = selectTopK_tradeoff  (resultsTable, k, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');

%% ------------------------- Report ----------------------------------
fprintf('\n========================================================\n');
fprintf('  TOP-%d  |  stage: %s\n', k, stageUsed);
fprintf('========================================================\n');

fprintf('\n----- TOP-%d por |TMOKE| (melhores TMOKE) -----\n', k);
printRanking(topTmoke, false);

fprintf('\n----- TOP-%d por |S| (melhores sensibilidades) -----\n', k);
printRanking(topSens, false);

fprintf('\n----- TOP-%d por TRADE-OFF (rank|TMOKE| + rank|S|) -----\n', k);
printRanking(topTradeoff, true);

fprintf('\n');

%% ===================== local functions =============================
function printRanking(T, showScore)
    % Compact, terminal-friendly listing of a ranked table.
    for i = 1:height(T)
        fprintf(['#%d | Ldom=%4.0f  Lden=%4.0f  h_si=%3.0f  h_au=%6.2f nm' ...
                 ' | |TMOKE|=%.5f @ alpha=%.3f deg | S=%+.4f deg/RIU'], ...
            i, T.L_domain_nm(i), T.l_dente_nm(i), T.h_si_nm(i), T.h_au_nm(i), ...
            T.maxAbsTMOKE_base(i), T.alpha_peak_base_deg(i), T.S_est_deg_per_RIU(i));
        if showScore && ismember('score_tradeoff', T.Properties.VariableNames)
            fprintf(' | score=%g', T.score_tradeoff(i));
        end
        fprintf('\n');
    end
end

% Selection helpers -- copied verbatim from senseAndTmoke_semCeyig.m to
% keep the exact same ranking behaviour.
function Tsel = selectTopK_tradeoff(T, K, colTM, colS)
    tm = abs(T.(colTM));
    ss = abs(T.(colS));
    rTM = tiedrank(-tm);
    rS  = tiedrank(-ss);
    score = rTM + rS;
    T.score_tradeoff = score;
    [~, ord] = sort(score, 'ascend');
    Tsel = T(ord(1:min(K,height(T))), :);
end

function Tsel = selectTopK_single(T, K, col)
    [~, ord] = sort(T.(col), 'descend');
    Tsel = T(ord(1:min(K,height(T))), :);
end

function Tsel = selectTopK_single_abs(T, K, col)
    [~, ord] = sort(abs(T.(col)), 'descend');
    Tsel = T(ord(1:min(K,height(T))), :);
end
