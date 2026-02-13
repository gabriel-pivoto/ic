function write_progress_xlsx(progressWorkbookPath, varargin)
% Write one or more tables to the Excel progress workbook (safe if Excel unavailable).
    for i = 1:2:numel(varargin)
        sheet = varargin{i};
        T = varargin{i+1};
        if ~isempty(T)
            try
                writetable(T, progressWorkbookPath, 'Sheet', sheet);
            catch
                warning('write_progress_xlsx: failed sheet %s', sheet);
            end
        end
    end
end
