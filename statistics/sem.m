function stderror = sem(data)
    % columnwise standard error
    stderror = nanstd(data)./sqrt(sum(~isnan(data)));
end