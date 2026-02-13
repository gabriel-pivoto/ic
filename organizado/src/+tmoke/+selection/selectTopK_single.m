function Tsel = selectTopK_single(T, K, col)
% Select TOP-K rows by the largest value of the specified column.
    [~, ord] = sort(T.(col), 'descend');
    Tsel = T(ord(1:min(K,height(T))), :);
end
