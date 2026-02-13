function Tsel = selectTopK_single_abs(T, K, col)
% Select TOP-K rows by the largest absolute value of the specified column.
    [~, ord] = sort(abs(T.(col)), 'descend');
    Tsel = T(ord(1:min(K,height(T))), :);
end
