function redirectAllNumericalsToTable(mdl, ttag)
% Redirect all numerical results to the specified table.
    try
        ntags = cell(mdl.result().numerical().tags());
        for k = 1:numel(ntags)
            try
                mdl.result().numerical(ntags{k}).set('table', ttag);
            catch
            end
        end
    catch
    end
end
