function plotSTDPcurve(latency_table, varargin)
% plotSTDPcurve(latency_table)
% Plot STDP-like curve using latency (ms) vs z-score.
% Optionally filters for coactivated pairs only.
% 
% OPTIONS:
%   'binSize'     : bin size in ms for latency, default = 5
%   'filterCoact' : if true, only show coactivated pairs (default: false)
%   'plotFit'     : fit a double-exponential curve (default: false)

% Parse inputs
p = inputParser;
p.addRequired('latency_table');
p.addParameter('binSize', 5);
p.addParameter('filterCoact', false);
p.addParameter('plotFit', false);
p.addParameter('plasticitySign', 'all');  % 'all', 'positive', or 'negative'
p.parse(latency_table, varargin{:});

binSize = p.Results.binSize;
filterCoact = p.Results.filterCoact;
plotFit = p.Results.plotFit;

% Filter by coactivation if requested
lat = latency_table.latency_ms;
z = latency_table.z_score;
if filterCoact
    mask = latency_table.coactivated == 1;
    lat = lat(mask);
    z = z(mask);
end

% Remove NaNs and apply plasticity sign filter
valid = ~isnan(lat) & ~isnan(z);

switch lower(p.Results.plasticitySign)
    case 'positive'
        valid = valid & z > 0;
    case 'negative'
        valid = valid & z < 0;
end
lat = lat(valid);
z = z(valid);

% Bin latencies
edges = min(lat):binSize:max(lat);
centers = edges(1:end-1) + binSize/2;

z_mean = zeros(size(centers));
z_sem = zeros(size(centers));

for k = 1:length(centers)
    in_bin = lat >= edges(k) & lat < edges(k+1);
    if sum(in_bin) > 0
        z_mean(k) = mean(z(in_bin));
        z_sem(k) = std(z(in_bin)) / sqrt(sum(in_bin));
    else
        z_mean(k) = NaN;
        z_sem(k) = NaN;
    end
end

% Plot
figure; hold on;
errorbar(centers, z_mean, z_sem, 'o-', 'LineWidth', 1.5);
xline(0, '--k'); yline(0, '--k');
xlabel('Latency (ms)'); ylabel('Plasticity (z-score)');
title('STDP-like relationship');
grid on;

% Optional fit (e.g., double exponential)
if plotFit
    fit_mask = ~isnan(z_mean);
    centers_fit = centers(fit_mask);
    z_fit = z_mean(fit_mask);
    % Fit a double exponential: A*exp(-x/tau1) if x<0, B*exp(x/tau2) if x>0
    % (simplified example, can be improved)
    plot(centers_fit, smooth(z_fit, 5), 'r-', 'LineWidth', 2);
end
end
