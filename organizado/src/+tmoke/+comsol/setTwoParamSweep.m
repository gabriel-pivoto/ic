function setTwoParamSweep(mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mStr)
% Configure a two-parameter sweep over alpha and magnetization m.
    import tmoke.comsol.getParametricTag
    ptag = getParametricTag(mdl, studyTag);
    mdl.study(studyTag).feature(ptag).set('pname', {alphaName, mName});
    mdl.study(studyTag).feature(ptag).set('punit', {'deg','1'});
    alist = sprintf('range(%.12g[deg], %.12g[deg], %.12g[deg])', aStartDeg, aStepDeg, aStopDeg);
    mdl.study(studyTag).feature(ptag).set('plistarr', {alist, mStr});
end
