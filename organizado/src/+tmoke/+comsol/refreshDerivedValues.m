function refreshDerivedValues(mdl)
% Refresh COMSOL numerical results to ensure tables are updated after a run.
    try
        ntags = cell(mdl.result().numerical().tags());
        for k = 1:numel(ntags)
            mdl.result().numerical(ntags{k}).setResult;
        end
    catch
    end
end
