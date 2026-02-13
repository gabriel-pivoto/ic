function save_checkpoint(checkpointFilePath, stage, runsCompletedGlobal, done_points, payload)
% Persist a checkpoint struct to disk (atomic save via .tmp).
    checkpointData.stage            = char(stage);
    checkpointData.runsCompletedGlobal = runsCompletedGlobal;
    checkpointData.done_points      = done_points;
    checkpointData.payload          = payload;
    tmp = [checkpointFilePath '.tmp'];
    save(tmp,'checkpointData','-v7.3');
    movefile(tmp, checkpointFilePath, 'f');
end
