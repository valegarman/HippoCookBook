function ID = computeIntrinsicDimension(rippleData, k)
    % Extract manifold data
    manifold = rippleData;
    
    % Compute intrinsic dimension
    ID = round(mean(abids(manifold, k)));
    fprintf('Intrinsic dimension = %d\n', ID);
end
