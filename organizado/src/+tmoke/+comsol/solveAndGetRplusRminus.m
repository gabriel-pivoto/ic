function [alpha_deg, Rplus, Rminus, TMOKE] = solveAndGetRplusRminus( ...
    mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, ...
    mPlusStr, mMinusStr, ttagPlus, ttagMinus)
% Run COMSOL sweeps for m = +1/-1, read reflectance tables, and compute TMOKE.

    import tmoke.comsol.*;

    ensureTable(mdl, ttagPlus);
    ensureTable(mdl, ttagMinus);

    Npts = 1 + floor((aStopDeg - aStartDeg)/aStepDeg + 1e-9);

    % m = +1
    redirectAllNumericalsToTable(mdl, ttagPlus);
    clearTable(mdl, ttagPlus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mPlusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha1_deg, R1] = readAlphaAndRFromNamedTable(mdl, ttagPlus, Npts);

    % m = -1
    redirectAllNumericalsToTable(mdl, ttagMinus);
    clearTable(mdl, ttagMinus);
    setTwoParamSweep(mdl, studyTag, alphaName, mName, aStartDeg, aStepDeg, aStopDeg, mMinusStr);
    mdl.study(studyTag).run;
    refreshDerivedValues(mdl);
    [alpha2_deg, R2] = readAlphaAndRFromNamedTable(mdl, ttagMinus, Npts);

    assert(numel(alpha1_deg)==Npts && numel(alpha2_deg)==Npts, 'Unexpected alpha sweep length.');
    assert(max(abs(alpha1_deg - alpha2_deg)) < 1e-9, 'alpha grids differ between m=+1 and m=-1.');

    alpha_deg = alpha1_deg;
    Rplus  = R1;
    Rminus = R2;

    denom = Rplus + Rminus;
    denom(abs(denom) < 1e-9) = 1e-9;
    TMOKE = (Rplus - Rminus) ./ denom;
end
