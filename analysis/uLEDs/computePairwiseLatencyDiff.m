function latency_table = computePairwiseLatencyDiff(spikes, uledPulses, varargin)
% latency_table = computePairwiseLatencyDiff(spikes, uledPulses, varargin)
% Estimate pairwise latencies based on peak lag in CCG around stimulation
% pulses. Returns a table with one row per neuron pair (i,j).
% 
% OPTIONS:
%   'win'      : time window relative to stim [s], default [0 0.025]
%   'binSize'  : bin size in seconds for CCG, default 0.001
%   'halfBins' : number of bins on each side of 0 in CCG, default 25
%   'coactivationMatrix' : optional NxN logical matrix of coactivated pairs
%   'zScoreMatrix'       : optional NxN matrix of z-scores from plasticity analysis

% Parse inputs
p = inputParser;
p.addRequired('spikes');
p.addRequired('uledPulses');
p.addParameter('win', [0 0.025]);
p.addParameter('binSize', 0.001);
p.addParameter('halfBins', 25);
p.addParameter('coactivationMatrix', []);
p.addParameter('zScoreMatrix', []);
p.addParameter('useFastRestrict', true);
p.parse(spikes, uledPulses, varargin{:});
win = p.Results.win;
binSize = p.Results.binSize;
halfBins = p.Results.halfBins;
coactivationMatrix = p.Results.coactivationMatrix;
zScoreMatrix = p.Results.zScoreMatrix;

nNeurons = length(spikes.times);

% 1. Crear intervalo total de análisis
stim_times = uledPulses.timestamps(:);
intervals = [stim_times + win(1), stim_times + win(2)];
intervals = sortrows(intervals);

% 2. Restricción de spikes
restrictedSpikes = spikes;

if p.Results.useFastRestrict
    fprintf('Using fast Restrict on global interval \n');
    all_bounds = [stim_times + win(1); stim_times + win(2)];
    global_interval = [min(all_bounds), max(all_bounds)];
    restrictedSpikes.times = cellfun(@(x) Restrict(x, global_interval), spikes.times, 'UniformOutput', false);
else
    fprintf('Using full Restrict with all intervals \n');
    restrictedSpikes.times = cellfun(@(x) Restrict(x, intervals), spikes.times, 'UniformOutput', false);
end

% 3. Calcular CCG global
session = loadSession;
Fs = 1 / session.extracellular.sr;
winSize = (halfBins * binSize) * 2;
[allCcg, t_ccg] = CCG(restrictedSpikes.times, [], ...
                     'binSize', binSize, ...
                     'duration', winSize, ...
                     'Fs', Fs, ...
                     'normtype', 'counts');

% 4. Extraer latencias de forma vectorizada
latency_ms = NaN(nNeurons, nNeurons);
[~, maxIdxMatrix] = max(allCcg, [], 1);  % max over time bins

for i = 1:nNeurons
    % fprintf('Analyzing neuron %d of %d \n', i, nNeurons);
    for j = 1:nNeurons
        if i == j
            continue;
        end
        idx = maxIdxMatrix(1, i, j);
        latency_ms(i,j) = t_ccg(idx) * 1000;
    end
end

% 5. Crear tabla a partir de pares válidos
[i_idx, j_idx] = find(~isnan(latency_ms) & ~eye(nNeurons));
lat_vals = latency_ms(sub2ind(size(latency_ms), i_idx, j_idx));

% Coactivación y zscore (si se proporcionan)
if ~isempty(coactivationMatrix)
    coact_vals = coactivationMatrix(sub2ind(size(coactivationMatrix), i_idx, j_idx));
else
    coact_vals = NaN(size(i_idx));
end

if ~isempty(zScoreMatrix)
    zscore_vals = zScoreMatrix(sub2ind(size(zScoreMatrix), i_idx, j_idx));
else
    zscore_vals = NaN(size(i_idx));
end

% 6. Crear tabla
latency_table = table(i_idx, j_idx, lat_vals, coact_vals, zscore_vals, ...
    'VariableNames', {'i','j','latency_ms','coactivated','z_score'});
end
