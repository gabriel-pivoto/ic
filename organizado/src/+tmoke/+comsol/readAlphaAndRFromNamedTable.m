function [alpha_deg, Rcol] = readAlphaAndRFromNamedTable(mdl, ttag, Npts)
% Read alpha and reflectance columns from a COMSOL table.
    S = mphtable(mdl, ttag);

    heads = [];
    if isfield(S,'colhead') && ~isempty(S.colhead), heads = string(S.colhead); end
    if isempty(heads) && isfield(S,'head') && ~isempty(S.head), heads = string(S.head); end
    if isempty(heads) && isfield(S,'header') && ~isempty(S.header), heads = string(S.header); end
    if isempty(heads) && isfield(S,'colnames') && ~isempty(S.colnames), heads = string(S.colnames); end

    if isempty(heads)
        if ~isfield(S,'data') || isempty(S.data) || size(S.data,2) < 2
            error('Table %s empty after run. Check Derived Values.', ttag);
        end
        a = S.data(:,1); R = S.data(:,2);
    else
        hlow = lower(heads);
        aIdx = find(contains(hlow,'alpha'), 1, 'first');
        rIdx = find(contains(hlow,'reflectance') | contains(hlow,'total reflectance') | contains(hlow,' total r'), 1, 'first');
        if isempty(aIdx) || isempty(rIdx)
            if ~isfield(S,'data') || isempty(S.data) || size(S.data,2) < 2
                error('Alpha/Reflectance not found in %s and not enough columns.', ttag);
            end
            aIdx = 1; rIdx = 2;
        end
        a = S.data(:, aIdx); R = S.data(:, rIdx);
    end

    assert(numel(a) >= Npts && numel(R) >= Npts, ...
        'Table %s has %d rows; expected >= %d.', ttag, numel(a), Npts);

    a = a(end-Npts+1:end);
    R = R(end-Npts+1:end);

    if ~isempty(heads)
        hlow = lower(heads);
        if any(contains(hlow, '(rad)')) || any(contains(hlow, '[rad]'))
            a = a * 180/pi;
        end
    end

    alpha_deg = a;
    Rcol      = R;
end
