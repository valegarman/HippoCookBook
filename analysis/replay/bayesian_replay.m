function [replayScores] = bayesian_replay(spikes,ripples,template,include,varargin)
% Estimates the Bayesian method for replay quantification
%
% USAGE
%
% [] = compareReplayMethods(spikes,ripples,template,include)
%
% INPUTS
%
% spikes  - buzcode cellinfo file. A structure that only requires three fields: 
%              -spikes.times, an Nx1 cell array of timestamps in seconds for each neuron
%              -spikes.UID, Nx1 vector of unique identifiers for each neuron in session
%              -spikes.spindices, a sorted vector of [spiketime UID], for all spikes. 
%                               useful input to some functions and plotting rasters
% ripples - buzcode ripples events file. A structure that only requires two fields:
%             -ripples.timestamps, an Nx2 matrix of start/stop times
%             -ripples.peaks, an Nx1 vector of peak timestamps for each event
% template -NxD matrix of N cells and D positions, average firing ratemaps
%           for each cell
% include - indices (1:N) of cells (i.e.place cells) to keep
%
% OUTPUTS 
% 
% bayesRadon - integral under the line of best fit, using the Radon
%              transform (Davidson 2009)
% bayesLinearWeighted - linear weighted correlation of the posterior
%              probability matrix (Grosmark 2016)
%
% HELP
%
% See bz_GetSpikes from the buzcode repo for help with the spikes/ripples
%   data structures
%
% Copyright (C) 2019 Adrien Peyrache & David Tingley.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
%
% Modify by Pablo Abad, NeuCompLab (2026). The original function is
% replay_Bayesian by Tingley and Peyrache.

% Parse options
p = inputParser;

addParameter(p,'nBinsThresh',5); % the minimum number of bins, per event, to analyze event (binSize * nBinsThresh = total event duration)
addParameter(p,'binSize',.01);
addParameter(p,'overlap',1);
addParameter(p,'saveMat',false);
addParameter(p,'n_shuffles',100);
addParameter(p,'plotting',false);

parse(p,varargin{:});

nBinsThresh = p.Results.nBinsThresh;
binSize = p.Results.binSize;
overlap = p.Results.overlap;
saveMat = p.Results.saveMat;
n_shuffles = p.Results.n_shuffles;
plotting = p.Results.plotting;

session = loadSession();

if isstruct(spikes)
    spkmat = bz_SpktToSpkmat(spikes.times,'overlap',overlap,'binSize',binSize * overlap); % converts spike timestamps into matrix format

else
    spkmat = bz_SpktToSpkmat(spikes,'overlap',overlap,'binSize',binSize * overlap); % converts spike timestamps into matrix format

end

if ~isempty(include)
    keep = intersect(include,find(sum(template')>0)); % uncomment to keep all HPC cells that fired
else
    keep = find(sum(template')>0); % just remove cells w/ zeros spikes in template
end

% Preallocating for speed
nEvents = size(ripples.timestamps,1);
bayesLinearWeighted = nan(nEvents,1);
bayesRadon = nan(nEvents,1);
slope_hpc = nan(nEvents,1);
Pr_temp = cell(nEvents,1);
prMax_temp = cell(nEvents,1);
binTimes_temp = cell(nEvents,1);
count_temp = cell(nEvents,1);
nCells = nan(nEvents,1);
nSpks  = nan(nEvents,1);
% shuffles
bayesLinearWeighted_cellID_shuf = cell(nEvents,1);
bayesLinearWeighted_circular_shuf = cell(nEvents,1);
slope_hpc_cellID_shuf = cell(nEvents,1);
slope_hpc_circular_shuf = cell(nEvents,1);
bayesRadon_cellID_shuf = cell(nEvents,1);
bayesRadon_circular_shuf = cell(nEvents,1);
dt = spkmat.dt;
template_keep = template(keep,:);
validCells = nan(nEvents,1);


% parpool
p = gcp('nocreate');
if isempty(p)
    parpool('threads');   % shared memory → ideal para spkmat grande
end

% Progress tracker
D = parallel.pool.DataQueue;
counter = 0;
afterEach(D, @updateProgress);

parfor event = 1:nEvents

    ts = ripples.timestamps(event,:);

    [data, counts, binTimes] = ...
        process_replayData(spkmat, ts, keep, binSize);

    if size(data,1) >= nBinsThresh

        %% REAL REPLAY
        [Pr, prMax] = placeBayes(data(:,keep), template_keep, dt);
        Pr(isnan(Pr)) = 0;

        bayesLinearWeighted(event) = ...
            makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1));

        [slope_hpc(event), bayesRadon(event)] = ...
            Pr2Radon(Pr');

        Pr_temp{event} = Pr;
        prMax_temp{event} = prMax;
        binTimes_temp{event} = binTimes;
        counts_temp{event} = counts;

        %% LOCAL ARRAYS (IMPORTANT FOR PARFOR)
        BLW_cell = zeros(n_shuffles,1);
        BLW_circ = zeros(n_shuffles,1);
        slope_cell = zeros(n_shuffles,1);
        slope_circ = zeros(n_shuffles,1);
        radon_cell = zeros(n_shuffles,1);
        radon_circ = zeros(n_shuffles,1);

        %% SHUFFLES
        for s = 1:n_shuffles

            % ---- cellID shuffle
            tempShuf = bz_shuffleCellID(template_keep);

            PrS = placeBayes(data(:,keep), tempShuf, dt);
            PrS(isnan(PrS)) = 0;

            BLW_cell(s) = makeBayesWeightedCorr1(PrS,ones(size(PrS,1),1));
            [slope_cell(s), radon_cell(s)] = Pr2Radon(PrS');

            % ---- circular shuffle
            tempShuf = bz_shuffleCircular(template_keep);

            PrS = placeBayes(data(:,keep), tempShuf, dt);
            PrS(isnan(PrS)) = 0;

            BLW_circ(s) = makeBayesWeightedCorr1(PrS,ones(size(PrS,1),1));
            [slope_circ(s), radon_circ(s)] = Pr2Radon(PrS');

        end

        %% SAVE LOCAL RESULTS
        bayesLinearWeighted_cellID_shuf{event} = BLW_cell;
        bayesLinearWeighted_circular_shuf{event} = BLW_circ;

        slope_hpc_cellID_shuf{event} = slope_cell;
        slope_hpc_circular_shuf{event} = slope_circ;

        bayesRadon_cellID_shuf{event} = radon_cell;
        bayesRadon_circular_shuf{event} = radon_circ;

    else
        Pr_temp{event} = [];
        prMax_temp{event} = [];
        binTimes_temp{event} = [];
    end

    %% EXTRA INFO
    nCells(event) = sum(sum(counts(:,keep))>0);
    nSpks(event)  = sum(sum(counts(:,keep)));

    % progress update
    send(D,1);

end

% Compute mean and std of shuffling
validCells = cellfun(@(x) ~isempty(x) && numel(x)==n_shuffles, ...
                     bayesLinearWeighted_cellID_shuf);

mu_weighted_linear = zeros(1,nEvents);
sigma_weighted_linear = zeros(1,nEvents);
mu_radon_linear = zeros(1,nEvents);
sigma_radon_linear = zeros(1,nEvents);
mu_weighted_circular = zeros(1,nEvents);
sigma_weighted_circular = zeros(1,nEvents);
mu_radon_circular = zeros(1,nEvents);
sigma_radon_circular = zeros(1,nEvents);

mu_weighted_linear(validCells) = cellfun(@mean, ...
                         bayesLinearWeighted_cellID_shuf(validCells));
sigma_weighted_linear(validCells) = cellfun(@std, ...
                            bayesLinearWeighted_cellID_shuf(validCells));
mu_radon_linear(validCells) = cellfun(@mean, ...
                         bayesRadon_cellID_shuf(validCells));
sigma_radon_linear(validCells) = cellfun(@std, ...
                            bayesRadon_cellID_shuf(validCells));

mu_weighted_circular(validCells) = cellfun(@mean, ...
                         bayesLinearWeighted_circular_shuf(validCells));
sigma_weighted_circular(validCells) = cellfun(@std, ...
                            bayesLinearWeighted_circular_shuf(validCells));

mu_radon_circular(validCells) = cellfun(@mean, ...
                         bayesRadon_circular_shuf(validCells));
sigma_radon_circular(validCells) = cellfun(@std, ...
                            bayesRadon_circular_shuf(validCells));

% stats

% create data struct to return
replayScores.bayesLinearWeighted = bayesLinearWeighted;
replayScores.bayesRadon = bayesRadon;

replayScores.bayesLinearWeighted_cellID_shuf = bayesLinearWeighted_cellID_shuf;
replayScores.bayesRadon_cellID_shuf = bayesRadon_cellID_shuf;
replayScores.bayesLinearWeighted_circular_shuf = bayesLinearWeighted_circular_shuf;
replayScores.bayesRadon_circular_shuf = bayesRadon_circular_shuf;

replayScores.bayesLinearWeighted_cellID_mean = mu_weighted_linear';
replayScores.bayesLinearWeighted_cellID_std = sigma_weighted_linear';
replayScores.bayesRadon_cellID_mean = mu_radon_linear';
replayScores.bayesRadon_cellID_std = sigma_radon_linear';
replayScores.bayesLinearWeighted_circular_mean = mu_weighted_circular';
replayScores.bayesLinearWeighted_circular_std = sigma_weighted_circular';
replayScores.bayesRadon_circular_mean = mu_radon_circular';
replayScores.bayesRadon_circular_std = sigma_radon_circular';

replayScores.Pr = Pr_temp;
replayScores.prMax = prMax_temp;
replayScores.binTimes = binTimes_temp;
replayScores.counts = counts_temp;

replayScores.nCells = nCells;
replayScores.nSpks = nSpks;
replayScores.mode = 'bayes';

% stats
replayScores.is_replay_linearWeighted_cellID = abs(((bayesLinearWeighted - mu_weighted_linear')./sigma_weighted_linear')) > 1.96;
replayScores.is_replay_linearRadon_cellID = abs(((bayesRadon - mu_radon_linear')./sigma_radon_linear')) > 1.96;
replayScores.is_replay_circularWeighted_cellID = abs(((bayesLinearWeighted - mu_weighted_circular')./sigma_weighted_circular')) > 1.96;
replayScores.is_replay_circularRadon_cellID = abs(((bayesRadon - mu_radon_circular')./sigma_radon_circular')) > 1.96;

% save variable
if saveMat
    save([session.general.name,'.bayesReplay.events.mat'],'replayScores');
end

if plotting
    figure;
    histogram(bayesLinearWeighted,50,'FaceColor',[.5 .6 .5], 'EdgeColor','k');
    hold on;
    histogram(bayesLinearWeighted(replayScores.is_replay_linearWeighted_cellID),50,'FaceColor',[.1 .1 .1], 'EdgeColor','k');
    legend('All ripples','Significant ripples');
    xlim([-1 1]);
    ylabel('Number of ripples');
    xlabel('Bayes Linear Weighted');
    saveas(gca,['SummaryFigures\replay_bayes_linearWeighted.png']);

    figure;
    histogram(bayesLinearWeighted,50,'FaceColor',[.5 .6 .5], 'EdgeColor','k');
    hold on;
    histogram(bayesLinearWeighted(replayScores.is_replay_circularWeighted_cellID),50,'FaceColor',[.1 .1 .1], 'EdgeColor','k');
    legend('All ripples','Significant ripples');
    xlim([-1 1]);
    ylabel('Number of ripples');
    xlabel('Bayes Circular Weighted');
    saveas(gca,['SummaryFigures\replay_bayes_circularWeighted.png']);

    figure;
    histogram(bayesRadon,50,'FaceColor',[.5 .6 .5], 'EdgeColor','k');
    hold on;
    histogram(bayesRadon(replayScores.is_replay_linearRadon_cellID),50,'FaceColor',[.1 .1 .1], 'EdgeColor','k');
    legend('All ripples','Significant ripples');
    ylabel('Number of ripples');
    xlabel('Bayes circular Radon');
    saveas(gca,['SummaryFigures\replay_bayes_linearRadon.png']);

    figure;
    histogram(bayesRadon,50,'FaceColor',[.5 .6 .5], 'EdgeColor','k');
    hold on;
    histogram(bayesRadon(replayScores.is_replay_circularRadon_cellID),50,'FaceColor',[.1 .1 .1], 'EdgeColor','k');
    legend('All ripples','Significant ripples');
    ylabel('Number of ripples');
    xlabel('Bayes circular Radon');
    saveas(gca,['SummaryFigures\replay_bayes_circularRadon.png']);


end

function updateProgress(~)
    
    counter = counter + 1;

    elapsed = toc;
    rate = counter / elapsed;
    remaining = (nEvents-counter)/rate;

    fprintf(['Ripple %d/%d | %.1f%% | ' ...
        'Elapsed %.1fs | ETA %.1fs\n'],...
        counter,nEvents,...
        100*counter/nEvents,...
        elapsed,remaining);
    end
end