function [replayScores] = rankOrder_replay(spikes,ripples,template,include,varargin)
% Estimates the rank order correlation method for replay quantification
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
% replayScore.rankOrd - rank order correlations w/ template for each event
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
% Created by Pablo Abad (NCL; 2026) based on replay_rankOrder by Peyrache
% and Tingley.

% Parse options
p = inputParser;

addParameter(p,'nCellsPerEvt',4); % the minimum number of bins, per event, to analyze event (binSize * nBinsThresh = total event duration)
addParameter(p,'binSize',.01);
addParameter(p,'overlap',1);
addParameter(p,'saveMat',false);
addParameter(p,'n_shuffles',100);
addParameter(p,'plotting',false);

parse(p,varargin{:});

nCellsPerEvt = p.Results.nCellsPerEvt;
binSize = p.Results.binSize;
overlap = p.Results.overlap;
saveMat = p.Results.saveMat;
n_shuffles = p.Results.n_shuffles;
plotting = p.Results.plotting;

spkmat = bz_SpktToSpkmat(spikes.times,'overlap',overlap,'binSize',binSize * overlap); % converts spike timestamps into matrix format

if ~isempty(include)
    keep = intersect(include,find(sum(template')>0)); % uncomment to keep all HPC cells that fired
else
    keep = find(sum(template')>0); % just remove cells w/ zeros spikes in template
end

% Preallocating for speed
nEvents = size(ripples.timestamps,1);
data_temp = cell(nEvents,1);
counts_temp = cell(nEvents,1);
binTimes_temp = cell(nEvents,1);
data_original_temp = cell(nEvents,1);

ord_template_temp = cell(nEvents,1);
ord_firstSpike_temp = cell(nEvents,1);
rankOrd_temp = cell(nEvents,1);
pvals_temp = cell(nEvents,1);


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
    parpool('threads');
end

% Progress tracker
D = parallel.pool.DataQueue;
counter = 0;
afterEach(D,@updateProgress);

for event = 1:nEvents

    ts = ripples.timestamps(event,:); % gets start/stop timestamp of ripple event
    [data counts binTimes data_original] = process_replayData(spkmat, ts, keep, binSize); % processing func to get out the 'event' using common FR heuristics
    
    idx = intersect(find(sum(data)>0),keep); % find cells that fired and are in the 'include' array (i.e. active place cells)
    if length(idx) >= nCellsPerEvt
        [~,~,ord_template] = sort_cells(template(idx,:));
        [~,ord_firstSpk] = sortrows(data(:,idx)','descend');
        [rankOrd(event),pvals(event)] = corr(ord_template,ord_firstSpk,'rows','complete');

        data_temp{event} = data;
        counts_temp{event} = counts;
        binTimes_temp{event} = binTimes;
        data_original_temp{event} = data_original;
        cell_idx_temp{event} = idx;
        % Shuffling
        rankOrd_shuf_cell = zeros(n_shuffles,1);
        pvals_shuff_cell = zeros(n_shuffles,1);
        for s = 1:n_shuffles
            [rankOrd_shuff_cell(s),pvals_shuff_cell(s)] = corr(ord_template(randperm(length(idx))),ord_firstSpk,'rows','complete');
        end
        
        ord_template_temp{event} = ord_template;
        ord_firstSpk_temp{event} = ord_firstSpk;
        rankOrd_shuff{event} = rankOrd_shuff_cell;
        pvals_shuff{event} = pvals_shuff_cell;
    else
        rankOrd(event) = nan;
        pvals(event) = nan;
        rankOrd_shuff{event} = nan;
        pvals_shuff{event} = nan;
    end       
    
    % extra data to send up
    nCells(event) = sum(sum(counts(:,keep))>0);
    nSpks(event) = sum(sum(counts(:,keep)));

    % Progress update
    send(D,1);
end

% Compute mean and std of shuffling
validCells = cellfun(@(x) ~isempty(x) && numel(x)==n_shuffles, ...
                     rankOrd_shuff);

rankOrd_mu = zeros(1,nEvents);
rankOrd_sigma = zeros(1,nEvents);

rankOrd_mu(validCells) = cellfun(@mean, ...
                         rankOrd_shuff(validCells));

rankOrd_sigma(validCells) = cellfun(@std, ...
                         rankOrd_shuff(validCells));

% create data struct to return
replayScores.rankOrd = rankOrd;
replayScores.pvals = pvals;
replayScores.rankOrd_shuf = rankOrd_shuff;
replayScores.pvals_shuf = pvals_shuff;
replayScores.ord_template = ord_template_temp;
replayScores.ord_firstSpk = ord_firstSpk_temp;

replayScores.data = data_temp;
replayScores.counts = counts_temp;
replayScores.binTimes = binTimes_temp;
replayScores.data_original = data_original_temp;
replayScores.cell_idx = cell_idx_temp;

replayScores.nCells = nCells;
replayScores.nSpks = nSpks;
replayScores.mode = 'rankOrder';

% stats
replayScores.is_replay = abs(((rankOrd - rankOrd_mu)./rankOrd_sigma)) > 1.96;

% save
if saveMat
    save([session.general.name,'rankOrdReplay.events.mat'],'replayScores');
end

if plotting
    figure;
    histogram(rankOrd,50,'FaceColor',[.5 .6 .5], 'EdgeColor','k');
    hold on;
    histogram(rankOrd(replayScores.is_replay),50,'FaceColor',[.1 .1 .1], 'EdgeColor','k');
    legend('All ripples','Significant ripples');
    xlim([-1 1]);
    ylabel('Number of ripples');
    xlabel('Rank Order Correlation');
    saveas(gca,['SummaryFigures\replay_rankOrder.png']);
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






