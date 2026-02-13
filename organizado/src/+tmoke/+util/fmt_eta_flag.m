function s = fmt_eta_flag(eta, is_exact)
% Format ETA with optional "(est.)" tag.
    import tmoke.util.fmt_time_long
    if ~isfinite(eta)
        s = 'N/A';
    else
        s = fmt_time_long(eta);
    end
    if ~is_exact
        s = [s ' (est.)'];
    end
end
