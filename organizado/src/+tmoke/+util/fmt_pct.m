function s = fmt_pct(frac)
% Format a fraction as percentage string with one decimal.
    if ~isfinite(frac)
        s = 'N/A';
    else
        s = sprintf('%5.1f', 100*frac);
    end
end
