function setParamScalar(mdl, name, val)
% Set a unitless COMSOL parameter (e.g., refractive index).
    mdl.param.set(name, sprintf('%.12g', val));
end
