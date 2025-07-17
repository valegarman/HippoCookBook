
function r = assortativity_by_cluster(W, Ci)
    n = length(Ci);
    same = Ci(:) == Ci(:)';
    upper = triu(true(n), 1);  % upper triangle to avoid double-counting
    weights = W(upper);
    same_clust = same(upper);
    
    r = corr(double(same_clust), weights, 'Type', 'Spearman');  % or Pearson
end