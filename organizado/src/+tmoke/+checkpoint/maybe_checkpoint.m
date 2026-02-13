function points_since_ckpt_new = maybe_checkpoint( ...
        stage, every_points, checkpointFilePath, progressWorkbookPath, ...
        runsCompletedGlobal, pointsSinceCheckpoint, done_points, payload, ...
        coarseResultsTable, fineResultsTable, superResultsTable, bestDenseTable, bestFullTable)
% Trigger periodic checkpoint saves and progress workbook updates.
    import tmoke.checkpoint.save_checkpoint
    import tmoke.checkpoint.write_progress_xlsx

    points_since_ckpt_new = pointsSinceCheckpoint;
    if pointsSinceCheckpoint >= every_points
        save_checkpoint(checkpointFilePath, stage, runsCompletedGlobal, done_points, payload);
        try
            write_progress_xlsx(progressWorkbookPath, ...
                'coarse', coarseResultsTable, 'fine', fineResultsTable, 'super', superResultsTable, ...
                'dense', bestDenseTable, 'full', bestFullTable);
        catch
        end
        points_since_ckpt_new = 0;
    end
end
