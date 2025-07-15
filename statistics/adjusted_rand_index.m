
function ari = adjusted_rand_index(labels_true, labels_pred)
    % Based on Hubert & Arabie 1985

    n = length(labels_true);
    contingency = crosstab(labels_true, labels_pred);
    sum_comb_c = sum(nchoosek_vector(sum(contingency,2)));
    sum_comb_k = sum(nchoosek_vector(sum(contingency,1)));
    sum_comb = sum(nchoosek_vector(contingency(:)));
    total_pairs = nchoosek(n,2);

    expected_index = sum_comb_c * sum_comb_k / total_pairs;
    max_index = (sum_comb_c + sum_comb_k) / 2;
    ari = (sum_comb - expected_index) / (max_index - expected_index);
end

function result = nchoosek_vector(v)
    result = arrayfun(@(x) nchoosek_safe(x,2), v);
end

function c = nchoosek_safe(n, k)
    if n >= k && mod(n,1) == 0
        c = nchoosek(n, k);
    else
        c = 0;
    end
end
