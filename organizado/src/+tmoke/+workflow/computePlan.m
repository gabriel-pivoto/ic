function plan = computePlan(cfg)
% Compute stage run counts and ETA planning constants from config.
plan = struct();

plan.numDomainPeriodPoints  = numel(cfg.domainPeriodGridNm);
plan.numToothWidthPoints    = numel(cfg.toothWidthGridNm);
plan.numSiliconHeightPoints = numel(cfg.siliconHeightGridNm);
plan.numCeyigHeightPoints   = numel(cfg.ceyigHeightGridNm);
plan.numGoldHeightPoints    = numel(cfg.goldHeightGridNm);

% Two runs (m = +/-1) for each refractive index sampled during the search.
plan.runsPerSearchPoint = 2 * numel(cfg.fastRefractiveIndexSamples);

plan.coarseTotalPoints = plan.numDomainPeriodPoints * plan.numToothWidthPoints * ...
    plan.numSiliconHeightPoints * plan.numCeyigHeightPoints * plan.numGoldHeightPoints;
plan.coarseTotalRuns = plan.runsPerSearchPoint * plan.coarseTotalPoints;

% Fine estimates before clamping to limits.
fineGoldHeightCount     = 1 + floor(2 * cfg.fineGoldHeightDelta     / cfg.fineGoldHeightStep);
fineCeyigHeightCount    = 1 + floor(2 * cfg.fineCeyigHeightDelta    / cfg.fineCeyigHeightStep);
fineDomainPeriodCount   = 1 + floor(2 * cfg.fineDomainPeriodDelta   / cfg.fineDomainPeriodStep);
fineToothWidthCount     = 1 + floor(2 * cfg.fineToothWidthDelta     / cfg.fineToothWidthStep);
fineSiliconHeightCount  = 1 + floor(2 * cfg.fineSiliconHeightDelta  / cfg.fineSiliconHeightStep);
plan.finePointsPerSeedEstimate = fineDomainPeriodCount * fineToothWidthCount * fineSiliconHeightCount * ...
    fineCeyigHeightCount * fineGoldHeightCount;
plan.fineTotalPointsEstimate = cfg.topKCoarse * plan.finePointsPerSeedEstimate;
plan.fineTotalRunsEstimate = plan.runsPerSearchPoint * plan.fineTotalPointsEstimate;

% Super counts are exact from the start.
superGoldHeightCount     = 1 + floor(2 * cfg.superGoldHeightDelta     / cfg.superGoldHeightStep);
superCeyigHeightCount    = 1 + floor(2 * cfg.superCeyigHeightDelta    / cfg.superCeyigHeightStep);
superDomainPeriodCount   = 1 + floor(2 * cfg.superDomainPeriodDelta   / cfg.superDomainPeriodStep);
superToothWidthCount     = 1 + floor(2 * cfg.superToothWidthDelta     / cfg.superToothWidthStep);
superSiliconHeightCount  = 1 + floor(2 * cfg.superSiliconHeightDelta  / cfg.superSiliconHeightStep);
plan.superPointsPerSeed = superDomainPeriodCount * superToothWidthCount * superSiliconHeightCount * ...
    superCeyigHeightCount * superGoldHeightCount;
plan.superTotalPoints = cfg.topKFine * plan.superPointsPerSeed;
plan.superTotalRuns = plan.runsPerSearchPoint * plan.superTotalPoints;

% Fixed extra runs (validation, optional snapshot, baseline full).
plan.extraRunsFixed = 2 * numel(cfg.validationRefractiveIndexList) + (cfg.SAVE_SNAPSHOT * 2) + 2;

plan.globalRunTargetEstimate = plan.coarseTotalRuns + plan.fineTotalRunsEstimate + ...
    plan.superTotalRuns + plan.extraRunsFixed;
plan.isTotalRunEstimateExact = false;

end
