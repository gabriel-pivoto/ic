function s = fmt_time_long(sec)
% Format seconds into d hh:mm:ss.
    if ~isfinite(sec) || sec < 0, s = '...'; return; end
    days = floor(sec/86400);
    rem  = sec - 86400*days;
    hrs  = floor(rem/3600);
    rem  = rem - 3600*hrs;
    mins = floor(rem/60);
    secs = floor(rem - 60*mins);
    if days > 0
        s = sprintf('%dd %02d:%02d:%02d', days, hrs, mins, secs);
    else
        s = sprintf('%02d:%02d:%02d', hrs, mins, secs);
    end
end
