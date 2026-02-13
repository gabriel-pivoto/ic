function tf = isvalidHandle(h)
% Small guard to test MATLAB graphics handles.
    tf = ~isempty(h) && isvalid(h);
end
