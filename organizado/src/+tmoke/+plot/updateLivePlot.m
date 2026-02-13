function updateLivePlot(stage, Ldom, Lden, hsi, hcey, hau, ...
                        alpha_deg, TMOKE_vec, alpha_peak, tmoke_peak, ...
                        Rplus, Rminus)
% Live plot updater for TMOKE and reflectance during sweeps.
    import tmoke.plot.isvalidHandle

    persistent fig ax hTM hPeak hRp hRm inited
    if isempty(inited) || ~isvalidHandle(fig) || ~isvalidHandle(ax) || ...
       isempty(hTM)   || ~isvalidHandle(hTM)   || ...
       isempty(hPeak) || ~isvalidHandle(hPeak) || ...
       isempty(hRp)   || ~isvalidHandle(hRp)   || ...
       isempty(hRm)   || ~isvalidHandle(hRm)

        fig = figure('Name','TMOKE vs \alpha (live)', ...
                     'NumberTitle','off','Tag','TMOKE_LIVE_FIG', 'Color','w');
        ax  = axes('Parent',fig,'Tag','TMOKE_LIVE_AX');
        grid(ax,'on'); hold(ax,'on');

        yyaxis(ax,'left');
        hTM   = plot(ax, nan, nan, '-', 'LineWidth', 1.2, 'DisplayName','TMOKE(\alpha)');
        hPeak = plot(ax, nan, nan, 'o', 'MarkerSize', 6, 'DisplayName','peak TMOKE');
        ylabel(ax,'TMOKE'); xlim(ax,[0 89]); xlabel(ax,'\alpha [deg]');

        yyaxis(ax,'right');
        hRp = plot(ax, nan, nan, '--', 'LineWidth', 1.0, 'DisplayName','R^+(\alpha)');
        hRm = plot(ax, nan, nan, ':',  'LineWidth', 1.0, 'DisplayName','R^-(\alpha)');
        ylabel(ax,'Reflectance');

        legend(ax,'Location','best');
        inited = true;
    end

    yyaxis(ax,'left');
    set(hTM,   'XData', alpha_deg, 'YData', TMOKE_vec);
    set(hPeak, 'XData', alpha_peak, 'YData', tmoke_peak);

    yyaxis(ax,'right');
    if nargin >= 12 && ~isempty(Rplus) && ~isempty(Rminus)
        set(hRp, 'XData', alpha_deg, 'YData', Rplus);
        set(hRm, 'XData', alpha_deg, 'YData', Rminus);
    else
        set(hRp, 'XData', nan, 'YData', nan);
        set(hRm, 'XData', nan, 'YData', nan);
    end

    title(ax, sprintf('%s | Ldom=%g | Lden=%g | hsi=%g | hcey=%.3f | hau=%.3f | |TM|=%.5f @ %.3f deg', ...
        stage, Ldom, Lden, hsi, hcey, hau, abs(tmoke_peak), alpha_peak));

    drawnow('limitrate');
end
