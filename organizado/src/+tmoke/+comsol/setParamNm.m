function setParamNm(mdl, name, val_nm)
% Set a COMSOL parameter carrying nanometer units.
    mdl.param.set(name, sprintf('%.12g[nm]', val_nm));
end
