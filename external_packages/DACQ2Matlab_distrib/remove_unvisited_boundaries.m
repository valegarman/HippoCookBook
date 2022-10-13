function map = remove_unvisited_boundaries(raw_map)

% Unvisited pixels must be defined as NaNs

nan_map = isnan(raw_map);
row1 = find(sum(nan_map,2)<size(nan_map,2), 1 );
rown = find(sum(nan_map,2)<size(nan_map,2), 1, 'last' );
col1 = find(sum(nan_map,1)<size(nan_map,1), 1 );
coln = find(sum(nan_map,1)<size(nan_map,1), 1, 'last' );

map = raw_map(row1:rown,col1:coln);

% ----------------------------------------------------------------------------------------------------------------------
