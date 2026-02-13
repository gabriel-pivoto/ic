function printSummary(bestTradeoffCandidate, bestTmokeCandidate, bestSensitivityCandidate, sensitivityDense, runsCompletedGlobal, elapsedSeconds)
% Final console summary for the TMOKE run.
import tmoke.util.*;

fprintf('\n===== SUMMARY =====\n');
fprintf('Best TRADE-OFF (SUPER):\n'); disp(bestTradeoffCandidate);
fprintf('Best TMOKE only (SUPER):\n'); disp(bestTmokeCandidate);
fprintf('Best |sensitivityEstimateFast| only (SUPER):\n'); disp(bestSensitivityCandidate);
fprintf('VALID @ bestTradeoffCandidate: sensitivityDense ~= %.6f deg/RIU\n', sensitivityDense);
fprintf('Runs done: %d | Elapsed: %s\n', runsCompletedGlobal, fmt_time_long(elapsedSeconds));
end
