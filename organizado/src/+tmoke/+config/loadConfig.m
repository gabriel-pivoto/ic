function cfg = loadConfig()
% Load TMOKE config: template + optional local overrides.
    templateDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'config');

    addpath(templateDir);
    cfg = tmoke_config_template();

    localCfgPath = fullfile(templateDir, 'tmoke_config_local.m');
    if exist(localCfgPath, 'file') == 2
        try
            localCfg = tmoke_config_local();
            cfg = mergeStructs(cfg, localCfg);
        catch ME
            warning('Failed to load tmoke_config_local.m: %s', ME.message);
        end
    end

    % Ensure derived paths exist in config
    if ~isfield(cfg,'checkpointDir') || isempty(cfg.checkpointDir)
        cfg.checkpointDir = fullfile(cfg.outputDir, 'checkpoints');
    end
    if ~isfield(cfg,'plotsRootDir') || isempty(cfg.plotsRootDir)
        cfg.plotsRootDir = fullfile(cfg.outputDir, 'plots');
    end
    if ~isfield(cfg,'progressWorkbookPath') || isempty(cfg.progressWorkbookPath)
        cfg.progressWorkbookPath = fullfile(cfg.checkpointDir, 'tmoke_sens_progress.xlsx');
    end

    % Derived refractive index baseline if missing
    if ~isfield(cfg,'baselineRefractiveIndex') || isempty(cfg.baselineRefractiveIndex)
        cfg.baselineRefractiveIndex = cfg.fastRefractiveIndexSamples(1);
    end
end

function merged = mergeStructs(base, override)
% Shallow merge: override fields replace base fields when present.
    merged = base;
    if isempty(override)
        return;
    end
    fields = fieldnames(override);
    for i = 1:numel(fields)
        merged.(fields{i}) = override.(fields{i});
    end
end
