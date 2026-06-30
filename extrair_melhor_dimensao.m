%% =====================================================================
% Extract the BEST geometry found SO FAR from a (possibly unfinished) run
% ---------------------------------------------------------------------
% Use case:
%   A search (senseAndTmoke_semCeyig.m) is running and is already in the
%   SUPER stage, but you need to stop it before it finishes. This script
%   reads the periodic checkpoint, rebuilds the results table for the most
%   advanced stage available, and reports the best geometry so far.
%
% It does NOT call COMSOL. It only reads the checkpoint .mat and applies
% the SAME selection logic as the main script (rank(|TMOKE|)+rank(|S|)).
%
% Output:
%   - prints the best trade-off / best |TMOKE| / best |S| geometries
%   - saves melhor_dimensao.mat (struct `best`) + melhor_dimensao.txt
%     so you can later run a separate plotting script using those numbers.
% =====================================================================
clear; clc; format long;

%% --------------------------- Config --------------------------------
projectRootDir     = 'C:\Users\gabri\Documents\projetoIC';
checkpointFilePath = fullfile(projectRootDir, 'checkpoints', ...
                              'tmoke_sens_4d_sem_ceyig_checkpoint.mat');

% Where to write the extracted geometry (same folder as this script).
outputMatPath = fullfile(fileparts(mfilename('fullpath')), 'melhor_dimensao.mat');
outputTxtPath = fullfile(fileparts(mfilename('fullpath')), 'melhor_dimensao.txt');

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

% array2table from a numeric matrix keeps it a table already; ensure table type.
if ~istable(resultsTable)
    resultsTable = array2table(resultsTable, 'VariableNames', columnNames);
end

fprintf('Using stage: %s | %d points evaluated so far.\n\n', stageUsed, height(resultsTable));

%% ------------------------- Selection -------------------------------
bestTradeoff = selectTopK_tradeoff(resultsTable, 1, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');
bestTmoke    = selectTopK_single (resultsTable, 1, 'maxAbsTMOKE_base');
bestSens     = selectTopK_single_abs(resultsTable, 1, 'S_est_deg_per_RIU');

fprintf('===== BEST SO FAR (%s) =====\n', stageUsed);
fprintf('-- Trade-off (rank|TMOKE| + rank|S|) --\n'); disp(bestTradeoff);
fprintf('-- Best |TMOKE| only --\n');                 disp(bestTmoke);
fprintf('-- Best |S| only --\n');                     disp(bestSens);

%% ------------------- Save the best geometry ------------------------
% Primary pick = trade-off candidate (same criterion the main script uses).
best = struct();
best.stageUsed   = stageUsed;
best.L_domain_nm = bestTradeoff.L_domain_nm(1);
best.l_dente_nm  = bestTradeoff.l_dente_nm(1);
best.h_si_nm     = bestTradeoff.h_si_nm(1);
best.h_au_nm     = bestTradeoff.h_au_nm(1);
best.alpha_peak_base_deg = bestTradeoff.alpha_peak_base_deg(1);
best.maxAbsTMOKE_base    = bestTradeoff.maxAbsTMOKE_base(1);
best.S_est_deg_per_RIU   = bestTradeoff.S_est_deg_per_RIU(1);

best.bestTradeoffCandidate = bestTradeoff;
best.bestTmokeCandidate    = bestTmoke;
best.bestSensitivityCandidate = bestSens;

save(outputMatPath, 'best');

fid = fopen(outputTxtPath, 'w');
fprintf(fid, 'BEST GEOMETRY SO FAR (%s)\n', stageUsed);
fprintf(fid, 'L_domain_nm = %g\n', best.L_domain_nm);
fprintf(fid, 'l_dente_nm  = %g\n', best.l_dente_nm);
fprintf(fid, 'h_si_nm     = %g\n', best.h_si_nm);
fprintf(fid, 'h_au_nm     = %g\n', best.h_au_nm);
fprintf(fid, 'alpha_peak_base_deg = %g\n', best.alpha_peak_base_deg);
fprintf(fid, 'maxAbsTMOKE_base    = %g\n', best.maxAbsTMOKE_base);
fprintf(fid, 'S_est_deg_per_RIU   = %g\n', best.S_est_deg_per_RIU);
fclose(fid);

fprintf('\nSaved best geometry to:\n  %s\n  %s\n', outputMatPath, outputTxtPath);
fprintf('Use these (L_domain, l_dente, h_si, h_au) in your plotting script.\n');

%% ===================== local functions =============================
% Copied verbatim from senseAndTmoke_semCeyig.m to keep the exact same
% selection behaviour.
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
