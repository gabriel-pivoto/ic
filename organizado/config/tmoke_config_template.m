function cfg = tmoke_config_template()
% Template configuration for TMOKE + sensitivity search.
% Copy/override via tmoke_config_local.m for machine-specific paths.

    templateDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(templateDir); % root = folder "pasta para essa organizaçaõ"

    cfg = struct();

    % Paths
    cfg.comsolProjectRootDir = projectRoot;                               % Where supporting scripts/models live
    cfg.comsolModelFile      = fullfile(projectRoot, 'usandoMatlab.mph'); % Path to COMSOL .mph model (edit this!)
    cfg.outputDir            = fullfile(projectRoot, 'output');           % Base output dir for CSV/plots/checkpoints
    cfg.checkpointDir        = fullfile(cfg.outputDir, 'checkpoints');
    cfg.plotsRootDir         = fullfile(cfg.outputDir, 'plots');
    cfg.progressWorkbookPath = fullfile(cfg.checkpointDir, 'tmoke_sens_progress.xlsx');

    % Run budget / checkpoint
    cfg.MAX_RUNS = 20000;
    cfg.checkpointEveryPoints = 10;

    % Refractive indices
    cfg.fastRefractiveIndexSamples = [1.33, 1.39];
    cfg.baselineRefractiveIndex    = cfg.fastRefractiveIndexSamples(1);
    cfg.validationRefractiveIndexList = [1.33, 1.36, 1.39];

    % Coarse grids (nm)
    cfg.domainPeriodGridNm = 800:50:850;
    cfg.toothWidthGridNm   = 500:50:600;
    cfg.siliconHeightGridNm = [220, 240, 260];
    cfg.ceyigHeightGridNm   = [100, 140];
    cfg.goldHeightGridNm    = 20:10:60;

    % Alpha steps
    cfg.alphaCoarseRange = [0, 1.0, 89];
    cfg.alphaFineStep  = 0.1;
    cfg.alphaSuperStep = 0.01;
    cfg.alphaDenseStep = 0.01;
    cfg.alphaFullStep  = 0.01;

    % Selection strategy
    cfg.topKCoarse = 1;
    cfg.topKFine   = 1;

    % Refinement windows
    cfg.fineGoldHeightDelta   = 2;    cfg.fineGoldHeightStep   = 1;
    cfg.fineCeyigHeightDelta  = 5;    cfg.fineCeyigHeightStep  = 1;
    cfg.fineDomainPeriodDelta = 10;   cfg.fineDomainPeriodStep = 5;
    cfg.fineToothWidthDelta   = 10;   cfg.fineToothWidthStep   = 5;
    cfg.fineSiliconHeightDelta= 5;    cfg.fineSiliconHeightStep= 5;

    cfg.superGoldHeightDelta  = 1;    cfg.superGoldHeightStep  = 0.5;
    cfg.superCeyigHeightDelta = 2;    cfg.superCeyigHeightStep = 1;
    cfg.superDomainPeriodDelta = 4;   cfg.superDomainPeriodStep = 2;
    cfg.superToothWidthDelta   = 4;   cfg.superToothWidthStep   = 2;
    cfg.superSiliconHeightDelta= 2;   cfg.superSiliconHeightStep= 1;

    % Alpha windows
    cfg.fineAlphaHalfSpan   = 5;
    cfg.superAlphaHalfSpan  = 2;
    cfg.fineAlphaHalfSpanSensitivity  = cfg.fineAlphaHalfSpan  + 1; % +/-6 deg
    cfg.superAlphaHalfSpanSensitivity = cfg.superAlphaHalfSpan + 2; % +/-4 deg

    % Output flags
    cfg.SAVE_SNAPSHOT = true;
    cfg.MAKE_PLOTS    = true;
    cfg.PLOT_LIVE     = true;
    cfg.SAVE_FIGS     = true;
    cfg.FIG_FORMATS   = {'png','pdf'};
end
