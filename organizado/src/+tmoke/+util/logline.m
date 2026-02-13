function logline(varargin)
% Print formatted log line and flush graphics.
    fprintf(varargin{:});
    drawnow('limitrate');
end
