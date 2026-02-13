function setAlphaMSweep(mdl, studyTag, alphaName, aStartDeg, aStepDeg, aStopDeg, mName, mListStr)
% Configure a sweep over alpha with a list of magnetization values.
    import tmoke.comsol.getParametricTag
    ptag = getParametricTag(mdl, studyTag);
    mdl.study(studyTag).feature(ptag).set('pname', {alphaName, mName});
    mdl.study(studyTag).feature(ptag).set('punit', {'deg','1'});
    alist = sprintf('range(%.12g[deg], %.12g[deg], %.12g[deg])', aStartDeg, aStepDeg, aStopDeg);
    mdl.study(studyTag).feature(ptag).set('plistarr', {alist, mListStr});
end
