function clearTable(mdl, ttag)
% Clear data from a COMSOL result table.
    try
        mdl.result().table(ttag).clearTableData;
    catch
    end
end
