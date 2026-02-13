function Tsel = selectTopK_tradeoff(T, K, colTM, colS)
% Select TOP-K by combined score: rank(|TM|) + rank(|S|)
% - rank 1 = best
    tm = abs(T.(colTM));
    ss = abs(T.(colS));

    rTM = tiedrank(-tm);
    rS  = tiedrank(-ss);

    score = rTM + rS;
    T.score_tradeoff = score; %#ok<AGROW>

    [~, ord] = sort(score, 'ascend');
    Tsel = T(ord(1:min(K,height(T))), :);
end
