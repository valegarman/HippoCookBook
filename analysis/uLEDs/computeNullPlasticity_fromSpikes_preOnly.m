function delta_FC_null = computeNullPlasticity_fromSpikes_preOnly(spikes, pre_interval, window_length, n_windows, n_shuffles, varargin)
% computeNullPlasticity_fromSpikes_preOnly
% Builds a null distribution of delta functional connectivity (∆FC)
% by comparing two sets of random spike intervals inside a pre-stimulation period.
%
% INPUTS:
%   - spikes: structure with spikes.times as {N x 1}
%   - pre_interval: [start_time, end_time] in seconds
%   - window_length: duration of each window (in seconds)
%   - n_windows: number of windows per FC computation (pre1 and pre2)
%   - n_shuffles: how many times to sample new sets of intervals
%   - varargin: passed to getFunctionalConnectivity (e.g. method, binSize, Restrict, etc.)
%
% OUTPUT:
%   - delta_FC_null: [N x N x n_shuffles] array of ∆FC matrices

% Get number of neurons
nNeurons = numel(spikes.times);
delta_FC_null = nan(nNeurons, nNeurons, n_shuffles);

% Unpack time limits
t_start = pre_interval(1);
t_end = pre_interval(2);
interval_range = t_end - t_start;

% Main loop
for s = 1:n_shuffles
    % Sample 2*n_windows random start times (allow overlap)
    random_starts = t_start + (interval_range - window_length) * rand(2 * n_windows, 1);
    all_intervals = [random_starts, random_starts + window_length];

    % Split into pre1 and pre2 sets
    intervals_pre1 = all_intervals(1:n_windows, :);
    intervals_pre2 = all_intervals(n_windows+1:end, :);

    % Compute FC for each set of windows
    FC_pre1 = getFunctionalConnectivity('spikes', spikes, 'Restrict', intervals_pre1, 'doPlot', false, 'saveMat', false, varargin{:});
    FC_pre2 = getFunctionalConnectivity('spikes', spikes, 'Restrict', intervals_pre2, 'doPlot', false, 'saveMat', false, varargin{:});

    % Store delta FC
    delta_FC_null(:,:,s) = FC_pre2.values - FC_pre1.values;
end

end