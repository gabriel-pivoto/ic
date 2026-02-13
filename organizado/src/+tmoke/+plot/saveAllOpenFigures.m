function saveAllOpenFigures(outDir, prefix, formats)
% Save all open figures to the specified directory with given formats.
    import tmoke.util.sanitizeFilename

    if ~exist(outDir,'dir'), mkdir(outDir); end
    figs = findall(0,'Type','figure');
    if isempty(figs)
        fprintf('>> No figures to save.\n');
        return;
    end

    [~, idx] = sort([figs.Number]);
    figs = figs(idx);

    for i = 1:numel(figs)
        fig = figs(i);
        figName = string(get(fig,'Name'));
        figTag  = string(get(fig,'Tag'));
        if strlength(figName)==0
            if strlength(figTag)>0
                figName = figTag;
            else
                figName = "Figure_" + string(fig.Number);
            end
        end

        base = string(prefix) + "__" + string(sanitizeFilename(figName));
        try
            set(fig,'Color','w');
        catch
        end

        for k = 1:numel(formats)
            ext = lower(string(formats{k}));
            filePath = fullfile(outDir, base + "." + ext);
            try
                switch ext
                    case "png"
                        exportgraphics(fig, filePath, 'Resolution', 300);
                    case "pdf"
                        exportgraphics(fig, filePath, 'ContentType','vector');
                    otherwise
                        saveas(fig, filePath);
                end
            catch ME
                warning('Failed to save %s: %s', filePath, ME.message);
            end
        end
    end

    fprintf('>> %d figure(s) saved in: %s\n', numel(figs), outDir);
end
