function ensureTable(mdl, ttag)
% Ensure a COMSOL result table exists.
    try
        mdl.result().table(ttag);
    catch
        try
            mdl.result().table().create(ttag, 'Table');
        catch
        end
    end
end
