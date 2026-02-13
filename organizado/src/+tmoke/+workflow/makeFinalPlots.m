function makeFinalPlots(cfg, paths, validation, coarseResultsTable, fineResultsTable, superResultsTable, bestTradeoffCandidate)
% Generate validation/trade-off figures and optionally save them.
import tmoke.plot.*;

if cfg.MAKE_PLOTS
    validationRefractiveIndexList = cfg.validationRefractiveIndexList;
    sensitivityDense = validation.sensitivityDense;
    alphaGridsByN = validation.alphaGridsByN;
    tmokeCurvesByN = validation.tmokeCurvesByN;
    alphaPeakDegreesByN = validation.alphaPeakDegreesByN;
    alphaVsNLinearFit = validation.alphaVsNLinearFit;
    tmokeMaxAbsByN = validation.tmokeMaxAbsByN;

    figure('Name','(1) TMOKE(alpha) for each n (VALID)','NumberTitle','off','Color','w');
    hold on; grid on;
    for i = 1:numel(validationRefractiveIndexList)
        plot(alphaGridsByN{i}, tmokeCurvesByN{i}, 'LineWidth', 1.2, 'DisplayName', sprintf('n=%.2f', validationRefractiveIndexList(i)));
    end
    xlabel('\alpha [deg]'); ylabel('TMOKE(\alpha)');
    title(sprintf('Best trade-off | TMOKE(\\alpha) vs \\alpha (S_{dense}=%.4f deg/RIU)', sensitivityDense));
    legend('Location','best');

    figure('Name','(2) alpha_peak vs n (VALID)','NumberTitle','off','Color','w');
    hold on; grid on;
    plot(validationRefractiveIndexList, alphaPeakDegreesByN, 'o', 'LineWidth', 1.5, 'DisplayName','data');
    nfit = linspace(min(validationRefractiveIndexList), max(validationRefractiveIndexList), 100);
    plot(nfit, polyval(alphaVsNLinearFit,nfit), '-', 'LineWidth', 1.2, 'DisplayName', sprintf('fit (S=%.4f)', sensitivityDense));
    xlabel('n'); ylabel('\alpha_{peak} [deg]');
    title('alpha_{peak} as a function of n');
    legend('Location','best');

    figure('Name','(3) |TMOKE|max vs n (VALID)','NumberTitle','off','Color','w');
    grid on; hold on;
    plot(validationRefractiveIndexList, tmokeMaxAbsByN, 'o-', 'LineWidth', 1.5);
    xlabel('n'); ylabel('|TMOKE|_{max}');
    title('|TMOKE|_{max} vs n');

    allCandidateResults = coarseResultsTable;
    if ~isempty(fineResultsTable),  allCandidateResults = [allCandidateResults; fineResultsTable]; end %#ok<AGROW>
    if ~isempty(superResultsTable), allCandidateResults = [allCandidateResults; superResultsTable]; end %#ok<AGROW>

    figure('Name','(4) Trade-off map: |TMOKE| vs |sensitivityEstimateFast|','NumberTitle','off','Color','w');
    grid on; hold on;
    scatter(abs(allCandidateResults.maxAbsTMOKE_base), abs(allCandidateResults.S_est_deg_per_RIU), 25, 'filled');
    scatter(abs(bestTradeoffCandidate.maxAbsTMOKE_base), abs(bestTradeoffCandidate.S_est_deg_per_RIU), 120, 'filled');
    xlabel('|TMOKE|_{max} (baseline)'); ylabel('|S_{est}| [deg/RIU]');
    title('Candidates (COARSE+FINE+SUPER) with bestTradeoffCandidate highlighted');
end

if cfg.SAVE_FIGS
    saveAllOpenFigures(paths.figureOutputDirectory, "tmoke_sens_5d", cfg.FIG_FORMATS);
    fprintf('Figures saved in: %s\n', paths.figureOutputDirectory);
end
end
