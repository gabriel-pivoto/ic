%% =====================================================================
% Draw the FINAL geometry as a COMSOL-like cross-section (schematic)
% ---------------------------------------------------------------------
% Reads the checkpoint written by senseAndTmoke_semCeyig.m, picks the
% final geometry (best trade-off by default), and draws a to-scale
% cross-section of the simplified structure:
%
%       Air  (superstrate)
%       Au   tooth  (l_dente wide, h_au tall)  -- centered per period
%       SiO2 block  (h_si tall)                -- substrate
%
% Period = L_domain. The figure mimics a COMSOL Graphics view (nm axes,
% material colors, legend, dimension callouts, incidence angle) and is
% saved to <projectRootDir>\graficos.
%
% It does NOT call COMSOL. It only reads the checkpoint .mat.
% =====================================================================
clear; clc; format long;

%% --------------------------- Config --------------------------------
projectRootDir     = 'D:\Gabriel Pivoto\projetoIC';
checkpointFilePath = fullfile(projectRootDir, 'checkpoints', ...
                              'tmoke_sens_4d_sem_ceyig_checkpoint.mat');
graphicsOutputDir  = fullfile(projectRootDir, 'graficos');

whichBest = 'tradeoff';   % 'tradeoff' | 'tmoke' | 'sens'
nCells    = 1;            % number of periods (teeth) to draw

% COMSOL-like material colors
colAir = [0.945 0.945 0.960];
colSiO2 = [0.62 0.80 0.92];
colAu  = [0.85 0.65 0.13];

columnNames = {'L_domain_nm','l_dente_nm','h_si_nm','h_au_nm', ...
               'maxAbsTMOKE_base','alpha_peak_base_deg','TMOKE_at_peak_base', ...
               'alpha_peak_high_deg','S_est_deg_per_RIU'};

%% ----------------------- Load checkpoint ---------------------------
if ~isfile(checkpointFilePath)
    error('Checkpoint not found: %s', checkpointFilePath);
end
S = load(checkpointFilePath, 'checkpointData');
cp = S.checkpointData; payload = cp.payload;

resultsTable = pickMostAdvancedTable(payload, columnNames);
stageUsed = resultsTable.Properties.UserData;   % stage label stashed below

%% ------------------------- Selection -------------------------------
switch lower(whichBest)
    case 'tradeoff'
        best = selectTopK_tradeoff(resultsTable, 1, 'maxAbsTMOKE_base', 'S_est_deg_per_RIU');
        bestLabel = 'melhor trade-off';
    case 'tmoke'
        best = selectTopK_single(resultsTable, 1, 'maxAbsTMOKE_base');
        bestLabel = 'melhor |TMOKE|';
    case 'sens'
        best = selectTopK_single_abs(resultsTable, 1, 'S_est_deg_per_RIU');
        bestLabel = 'melhor sensibilidade';
    otherwise
        error('whichBest must be tradeoff | tmoke | sens');
end

Ldom  = best.L_domain_nm(1);        % period [nm]
Lden  = best.l_dente_nm(1);         % gold tooth width [nm]
hSi   = best.h_si_nm(1);            % SiO2 block height [nm]
hAu   = best.h_au_nm(1);            % gold tooth height [nm]
alpha = best.alpha_peak_base_deg(1);
tmk   = best.maxAbsTMOKE_base(1);
Sens  = best.S_est_deg_per_RIU(1);

fprintf('Geometria final (%s) | stage: %s\n', bestLabel, stageUsed);
fprintf('  L_domain = %.2f nm | l_dente = %.2f nm | h_si = %.2f nm | h_au = %.2f nm\n', ...
        Ldom, Lden, hSi, hAu);
fprintf('  |TMOKE| = %.5f @ alpha = %.3f deg | S = %+.4f deg/RIU\n', tmk, alpha, Sens);

%% ------------------------- Draw figure -----------------------------
fig = figure('Color','w','Position',[100 100 760 720]);
ax  = axes(fig); hold(ax,'on'); box(ax,'on');

xHalf   = nCells*Ldom/2;
xLeft   = -xHalf;  xRight = xHalf;
airTop  = hSi + hAu + max([0.55*hSi, 4*hAu, 250]);

% --- material patches (order: air background -> substrate -> teeth) ---
rectangle(ax,'Position',[xLeft 0 (xRight-xLeft) airTop], ...
    'FaceColor',colAir,'EdgeColor','none');                          % Air
rectangle(ax,'Position',[xLeft 0 (xRight-xLeft) hSi], ...
    'FaceColor',colSiO2,'EdgeColor',[0.35 0.35 0.45],'LineWidth',0.75); % SiO2

for k = 0:(nCells-1)
    xc = xLeft + (k+0.5)*Ldom;                 % tooth center in this cell
    rectangle(ax,'Position',[xc-Lden/2, hSi, Lden, hAu], ...
        'FaceColor',colAu,'EdgeColor',[0.4 0.3 0.05],'LineWidth',0.75); % Au tooth
end

% --- dimension callouts (all placed in the AIR, clear of the shapes) --
cx = xLeft + floor(nCells/2)*Ldom;             % left edge of central cell
drawHDim(ax, cx, cx+Ldom, -0.14*hSi,            sprintf('L_{domain} = %.0f nm', Ldom), false); % below block
drawHDim(ax, -Lden/2, Lden/2, hSi+hAu+0.18*hSi, sprintf('l_{dente} = %.0f nm', Lden), true);   % above tooth
drawVDim(ax, xLeft-0.06*Ldom, 0, hSi,           sprintf('h_{si} = %.0f nm', hSi),  xLeft-0.13*Ldom);   % left of block
drawVDim(ax, Lden/2+0.10*Ldom, hSi, hSi+hAu,    sprintf('h_{au} = %.1f nm', hAu),  Lden/2+0.17*Ldom);  % right of tooth

% --- material legend (placed OUTSIDE the axes to avoid overlap) -------
hAir  = patch(ax,NaN,NaN,colAir ,'EdgeColor',[0.5 0.5 0.5],'DisplayName','Air');
hSiO2 = patch(ax,NaN,NaN,colSiO2,'EdgeColor',[0.35 0.35 0.45],'DisplayName','SiO_2');
hGold = patch(ax,NaN,NaN,colAu  ,'EdgeColor',[0.4 0.3 0.05],'DisplayName','Au');
legend(ax,[hGold hSiO2 hAir],'Location','northeastoutside','FontSize',10);

% --- axes cosmetics ---------------------------------------------------
axis(ax,'equal');
xlim(ax,[xLeft-0.30*Ldom, xRight+0.30*Ldom]);
ylim(ax,[-0.30*hSi, airTop]);
xlabel(ax,'x [nm]'); ylabel(ax,'y [nm]');
title(ax, {sprintf('Geometria final (%s) \\alpha_{peak} = %.3f\\circ', bestLabel, alpha), ...
    sprintf('|TMOKE| = %.5f   |   S = %+.4f deg/RIU', tmk, Sens)}, ...
    'FontWeight','bold');
set(ax,'Layer','top','FontSize',10);

%% ------------------------- Save PNG --------------------------------
if ~exist(graphicsOutputDir,'dir'), mkdir(graphicsOutputDir); end
stamp = datestr(now,'yyyymmdd_HHMMSS');
fileBase = sprintf('geometria_final_%s_Ldom%.0f_Lden%.0f_hsi%.0f_hau%.1f_%s', ...
    lower(whichBest), Ldom, Lden, hSi, hAu, stamp);
pngPath = fullfile(graphicsOutputDir, [fileBase '.png']);
figPath = fullfile(graphicsOutputDir, [fileBase '.fig']);

try
    exportgraphics(fig, pngPath, 'Resolution', 300);
catch
    print(fig, pngPath, '-dpng', '-r300');
end
try, savefig(fig, figPath); catch, end

fprintf('\nFigura salva em:\n  %s\n', pngPath);

%% ===================== local functions =============================
function T = pickMostAdvancedTable(payload, columnNames)
    % Priority: SUPER -> FINE -> COARSE (finished or partial rows).
    T = []; stageUsed = '';
    if isfield(payload,'superResultsTable') && ~isempty(payload.superResultsTable)
        T = payload.superResultsTable;                              stageUsed = 'SUPER (finished)';
    elseif isfield(payload,'superRows') && ~isempty(payload.superRows)
        T = array2table(payload.superRows,'VariableNames',columnNames); stageUsed = 'SUPER (partial)';
    elseif isfield(payload,'fineResultsTable') && ~isempty(payload.fineResultsTable)
        T = payload.fineResultsTable;                               stageUsed = 'FINE (finished)';
    elseif isfield(payload,'fineRows') && ~isempty(payload.fineRows)
        T = array2table(payload.fineRows,'VariableNames',columnNames);   stageUsed = 'FINE (partial)';
    elseif isfield(payload,'coarseResultsTable') && ~isempty(payload.coarseResultsTable)
        T = payload.coarseResultsTable;                             stageUsed = 'COARSE (finished)';
    elseif isfield(payload,'coarseRows') && ~isempty(payload.coarseRows)
        T = array2table(payload.coarseRows,'VariableNames',columnNames); stageUsed = 'COARSE (partial)';
    else
        error('No usable results found in checkpoint payload.');
    end
    if ~istable(T), T = array2table(T,'VariableNames',columnNames); end
    T.Properties.UserData = stageUsed;
end

function drawHDim(ax, x1, x2, y, label, above)
    % Horizontal dimension line with end ticks; label sits above/below it.
    tick = 0.03*abs(x2-x1) + eps;
    plot(ax,[x1 x2],[y y],'k-','LineWidth',1);
    plot(ax,[x1 x1],[y-tick y+tick],'k-','LineWidth',1);
    plot(ax,[x2 x2],[y-tick y+tick],'k-','LineWidth',1);
    if above, va = 'bottom'; else, va = 'top'; end
    text(ax,(x1+x2)/2, y, label, 'HorizontalAlignment','center', ...
        'VerticalAlignment',va,'FontSize',9,'BackgroundColor','none','Margin',1);
end

function drawVDim(ax, x, y1, y2, label, textX)
    % Vertical dimension line at x with end ticks; rotated label at textX.
    tick = 0.10*abs(x-textX) + eps;
    plot(ax,[x x],[y1 y2],'k-','LineWidth',1);
    plot(ax,[x-tick x+tick],[y1 y1],'k-','LineWidth',1);
    plot(ax,[x-tick x+tick],[y2 y2],'k-','LineWidth',1);
    text(ax, textX, (y1+y2)/2, label, 'HorizontalAlignment','center', ...
        'VerticalAlignment','middle','Rotation',90,'FontSize',9, ...
        'BackgroundColor','none','Margin',1);
end

% Selection helpers -- identical to senseAndTmoke_semCeyig.m
function Tsel = selectTopK_tradeoff(T, K, colTM, colS)
    rTM = tiedrank(-abs(T.(colTM)));
    rS  = tiedrank(-abs(T.(colS)));
    T.score_tradeoff = rTM + rS;
    [~, ord] = sort(T.score_tradeoff, 'ascend');
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
