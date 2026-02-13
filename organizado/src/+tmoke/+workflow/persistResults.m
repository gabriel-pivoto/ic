function persistResults(paths, coarseResultsTable, fineResultsTable, superResultsTable, bestDenseTable, bestFullTable)
% Save CSV/XLSX outputs and clear checkpoint when done.
import tmoke.checkpoint.*;
if ~isempty(coarseResultsTable), writetable(coarseResultsTable, fullfile(paths.projectRootDir,'tmoke_sens_5D_coarse.csv')); end
if ~isempty(fineResultsTable),   writetable(fineResultsTable,   fullfile(paths.projectRootDir,'tmoke_sens_5D_fine.csv')); end
if ~isempty(superResultsTable),  writetable(superResultsTable,  fullfile(paths.projectRootDir,'tmoke_sens_5D_super.csv')); end

writetable(bestDenseTable, fullfile(paths.projectRootDir,'tmoke_sens_bestTradeoff_dense_ALLn.csv'));
writetable(bestFullTable,  fullfile(paths.projectRootDir,'tmoke_sens_bestTradeoff_full_baseline.csv'));

write_progress_xlsx(paths.progressWorkbookPath, 'full', bestFullTable);

if isfile(paths.checkpointFilePath), delete(paths.checkpointFilePath); end
end
