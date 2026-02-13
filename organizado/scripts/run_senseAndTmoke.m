function run_senseAndTmoke()
% Entry point to run the TMOKE + sensitivity search with packaged helpers.
    runnerDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(runnerDir);

    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'config'));

    cfg = tmoke.config.loadConfig();
    tmoke.runSenseAndTmoke(cfg);
end
