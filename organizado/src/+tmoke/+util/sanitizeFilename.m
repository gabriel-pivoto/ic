function s = sanitizeFilename(str)
% Replace unsafe filename characters with underscores.
    s = regexprep(char(str), '[^\\w\\d\\-]+', '_');
    s = regexprep(s, '_+', '_');
    s = regexprep(s, '^_|_$', '');
end
