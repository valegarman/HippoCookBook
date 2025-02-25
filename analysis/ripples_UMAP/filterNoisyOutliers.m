function noiseIdx = filterNoisyOutliers(data)
    % Compute pairwise distances
    D = pdist2(data, data);
    
    % Set diagonal to NaN
    D(1:size(D,1)+1:end) = NaN;
    
    % Calculate nearest neighbor distances
    nnDist = sum(D < prctile(D(:), 5), 2);
    
    % Identify noisy points
    noiseIdx = nnDist < prctile(nnDist, 20);
end