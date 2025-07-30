function plasticity_output = runPlasticityRule_uLED(spikes, neuronLED_matrix, uledPulses, varargin)
% runPlasticityRule_uLED computes STDP-like plasticity effects using coactivation between neurons and uLED pulses
%
% INPUTS:
%   - spikes: struct from buzcode, with spikes.times{n}
%   - neuronLED_matrix: NxM logical matrix (N neurons x M uLEDs), true if neuron responds to a given uLED
%   - uledPulses: struct with:
%         - uledPulses.timestamps: [start end] timestamps of pulses (in seconds)
%         - uledPulses.code: Mx1 vector with uLED code per pulse
%
% OPTIONALS:
%   - prePostWindow: time window to consider pre- and post-spikes [s], default [0.005 0.025]
%   - binSize: temporal binning [s], default 0.001
%   - doPlot: whether to show plots, default true
%
% OUTPUT:
%   - plasticity_output: struct with prePostHist, t, and coactivePairs
%
% MV-NCL 2025

% Options
p = inputParser;
addParameter(p, 'prePostWindow', [0.005 0.025], @isnumeric);
addParameter(p, 'binSize', 0.001, @isnumeric);
addParameter(p, 'doPlot', true, @islogical);
addParameter(p, 'padding', 0.01, @isnumeric);  % NEW
parse(p, varargin{:});
prePostWindow = p.Results.prePostWindow;
binSize = p.Results.binSize;
doPlot = p.Results.doPlot;
padding = p.Results.padding;

nNeurons = length(spikes.times);
nLEDs = size(neuronLED_matrix, 2);
nBins = round((prePostWindow(2) + prePostWindow(1)) / binSize);
edges = -prePostWindow(1):binSize:prePostWindow(2);
t = edges(1:end-1) + binSize/2;

% Initialize
prePostHist = zeros(nNeurons, nNeurons, length(t));
coactivePairs = false(nNeurons, nNeurons);

% Loop over pulses
for iPulse = 1:length(uledPulses.code)
    led_id = uledPulses.code(iPulse);
    if led_id < 1 || led_id > nLEDs
        warning('LED code %d out of range, skipping...', led_id);
        continue;
    end

    active_neurons = find(neuronLED_matrix(:, led_id));
    if numel(active_neurons) < 2
        continue;  % no pairs to test
    end

    stim_time = mean(uledPulses.timestamps(iPulse,:));
    
    for i = 1:numel(active_neurons)
        for j = 1:numel(active_neurons)
            if i == j, continue; end
            n1 = active_neurons(i);
            n2 = active_neurons(j);

            % Spike times relative to stim
            t1 = spikes.times{n1} - stim_time;
            t2 = spikes.times{n2} - stim_time;

            % Only spikes in (window + padding)
            t1_sel = t1(t1 >= -prePostWindow(1)-padding & t1 <= prePostWindow(2)+padding);
            t2_sel = t2(t2 >= -prePostWindow(1)-padding & t2 <= prePostWindow(2)+padding);
            if isempty(t1_sel) || isempty(t2_sel)
                continue;
            end

            % All pairwise diffs: t2 - t1 (i.e., target - source)
            rel_times = [];
            for spike1 = t1_sel(:)'
                rel_times = [rel_times, t2_sel - spike1];
            end

            % Histogram
            h = histcounts(rel_times, edges);

            prePostHist(n1, n2, :) = squeeze(prePostHist(n1, n2, :)) + h(:);
            coactivePairs(n1, n2) = true;
        end
    end
end

% Output
plasticity_output.prePostHist = prePostHist;
plasticity_output.t = t;
plasticity_output.coactivePairs = coactivePairs;

% Plot if requested
if doPlot
    avgHist = squeeze(mean(prePostHist(coactivePairs), 1));
    figure;
    plot(t, avgHist, 'k', 'LineWidth', 2);
    hold on
    xline(0, '--r');
    xlabel('Time (s)');
    ylabel('Coactivation count');
    title('Average pre-post spike timing');
end
end