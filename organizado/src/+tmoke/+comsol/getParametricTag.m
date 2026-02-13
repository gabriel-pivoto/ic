function ptag = getParametricTag(mdl, studyTag)
% Find the parametric sweep feature tag for the given study.
    ptag = 'param';
    try
        fts = cell(mdl.study(studyTag).feature().tags());
        for i = 1:numel(fts)
            typ = char(mdl.study(studyTag).feature(fts{i}).featureType());
            if contains(lower(typ), 'param')
                ptag = fts{i};
                break;
            end
        end
    catch
    end
end
