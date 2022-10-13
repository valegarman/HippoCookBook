function [fpositions, fspikes, frate] = smooth_field_plot(positions,spikes,sz,type)

% Returns firing field with smoothing before dividing spikes by positions
% NB doesn't assume zero for unvisited positions: just ignores them.

switch type
    case 'boxcar'
        b = ones(sz);
        % Other filters could be defined in here: should work the same, but might need denom correction
end

c = ones(size(positions));
c(positions==0) = 0;
denom = filter2(b, c);
denom(denom==0) = NaN;
fpositions = filter2(b, positions);
fpositions = fpositions./denom;

if not(isempty(spikes))
    fspikes = filter2(b, spikes);
    fspikes = fspikes./denom;
else
    fspikes = [];
end

if not(isempty(spikes))
    frate = fspikes./fpositions;
else
    frate=[];
end

% Set field = NaN in unoccupied locations
fpositions(positions==0) = NaN;
frate(positions==0) = NaN;
fspikes(positions==0) = NaN;

% ----------------------------------------------------------------------------------------------------------------------
