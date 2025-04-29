
function [spatialModulation] = getSpatialModulation(varargin)
% [spatialModulation] = getSpatialModulation(varargin)
%
% Get spatial modulation structure
%
% <OPTIONALS>
% firingMaps        Buzcode firingMap structure. By default, looks in basepath.
% placeFieldStats   Buzcode placeFieldStats structure. By default, looks in basepath.
% spikes            Buzcode spikes structure. By default, looks in basepath.
% behaviour         Buzcode behaviour structure. By default, looks in basepath.
% firingTrialsMap   Buzcode firingTrialsMap structure. By default, looks in basepath.
% basepath          Default, pwd
% saveMat   	    Saves file, logical (default: true) 
%
% OUTPUTS
% spatialModulation
%
% MV-BuzsakiLab 2022
%
% to do: For now it only uses the biggest field if a neuron has more than
% one field. 
% to do: summary plot
% Parse options
p = inputParser;
addParameter(p,'firingMaps',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'behaviour',[],@isstruct);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'placeFieldStats',[],@isstruct);
addParameter(p,'firingTrialsMap',[],@isstruct);
addParameter(p,'saveMat', true, @islogical);
addParameter(p,'plotOpt', true, @islogical);
addParameter(p,'force', false, @islogical);

parse(p, varargin{:});
firingMaps = p.Results.firingMaps;
basepath = p.Results.basepath;
behaviour = p.Results.behaviour;
spikes = p.Results.spikes;
placeFieldStats = p.Results.placeFieldStats;
saveMat = p.Results.saveMat;
plotOpt = p.Results.plotOpt;
firingTrialsMap = p.Results.firingTrialsMap;
force = p.Results.force;

% Deal with inputs
prevPath = pwd;
cd(basepath);

filename = basenameFromBasepath(pwd);
if ~isempty(dir([basenameFromBasepath(pwd) '.spatialModulation.cellinfo.mat'])) && ~force
    disp('Spatial modulation already computed! Loading file');
    file =dir([basenameFromBasepath(pwd) '.spatialModulation.cellinfo.mat']);
    load(file.name);
    return
end

if isempty(behaviour)
    behaviour = getSessionLinearize('forceReload',false);  
end

if isempty(spikes)
    spikes = loadSpikes;
end

if isempty(firingMaps)
    if exist([basenameFromBasepath(pwd) '.firingMapsAvg.cellinfo.mat']) == 2
        load([basenameFromBasepath(pwd) '.firingMapsAvg.cellinfo.mat']);
    else
        try
            firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',false);
        catch
            error('Firing maps could not be obtained.')
        end
    end
end

if isempty(placeFieldStats)
    if exist([basenameFromBasepath(pwd) '.placeFields.cellinfo.mat']) == 2
        load([basenameFromBasepath(pwd) '.placeFields.cellinfo.mat']);
    else
        try
            placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps,'maxSize',.75,'sepEdge',0.03); %
        catch
            error('Place fields could not be obtained.')
        end
    end
end

if isempty(firingTrialsMap)
    firingTrialsMap = firingMapPerTrial;
end

% rearrange maps
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
for jj = 1:length(firingMaps.rateMaps{1})
    for ii = 1:length(firingMaps.UID)
        spatialModulation.(['map_' num2str(jj) '_rateMaps'])(ii,:) = firingMaps.rateMaps{ii}{jj};
        spatialModulation.(['map_' num2str(jj) '_countMaps'])(ii,:) = firingMaps.countMaps{ii}{jj};
        spatialModulation.(['map_' num2str(jj) '_occupancy'])(ii,:) = firingMaps.occupancy{ii}{jj};
        spatialModulation.(['map_' num2str(jj) '_rateMapsZ'])(ii,:) = zscor_xnan(firingMaps.rateMaps{ii}{jj});
    end
    spatialModulation.(['map_' num2str(jj) '_timestamps'])...
        = 0:firingMaps.cmBin:firingMaps.cmBin*(length(firingMaps.rateMaps{1}{jj})-1);
end

%% compute statistics

% 1 Spatial corr
for ii = 1:length(firingMaps.UID)
    maps_pairs = nchoosek(1:length(firingMaps.rateMaps{ii}), 2);
    for jj = 1:size(maps_pairs,1) % for each pairs of maps
        % 1 Spatial corr
        [rho, pval]= corr(firingMaps.rateMaps{ii}{maps_pairs(jj,1)}',firingMaps.rateMaps{ii}{maps_pairs(jj,2)}',...
            'Type','Spearman','rows','complete');
        spatialModulation.(['corr_maps_' num2str(maps_pairs(jj,1)) '_' num2str(maps_pairs(jj,2))])(ii,1) = rho;
        spatialModulation.(['corr_pval_maps_' num2str(maps_pairs(jj,1)) '_' num2str(maps_pairs(jj,2))])(ii,1) = pval;
        clear rho pval
    end
end

% 2 Spatial stats
for jj = 1:length(firingMaps.rateMaps{1})
    for ii = 1:length(firingMaps.UID)
        sTimeSpent = firingMaps.occupancy{ii}{jj};
        snSpikes = firingMaps.countMaps{ii}{jj};
        Rate_Map = firingMaps.rateMaps{ii}{jj};
        
        T = sum(sum(sTimeSpent));
        occupancy = sTimeSpent./T;
        meanFiringRate = sum(sum(Rate_Map.*sTimeSpent))./T;
        logArg = Rate_Map./meanFiringRate;
        logArg(logArg == 0) = 1;
        
        bits_spike = sum(sum(occupancy.*logArg.*log2(logArg))); % bits per spike.
        bits_second = sum(sum(occupancy.*Rate_Map.*log2(logArg))); % bits per second.
        sparsity = ((sum(sum(occupancy.*Rate_Map))).^2)/sum(sum(occupancy.*(Rate_Map.^2)));
        selectivity = max(max(Rate_Map))./meanFiringRate;
        
        pf_position = placeFieldStats.mapStats{ii}{jj}.x(1);
        pf_peak = placeFieldStats.mapStats{ii}{jj}.peak(1);
        if pf_peak == 0
            pf_peak =  NaN;
        end
        pf_size = placeFieldStats.mapStats{ii}{jj}.size(1); 
        if pf_size == 0
            pf_size = NaN;
        end
        fieldX = placeFieldStats.mapStats{ii}{jj}.fieldX(1,:);
        is_placeField = double(~isnan(placeFieldStats.mapStats{ii}{jj}.x(1)));

        spatialModulation.(['meanRate_map_' num2str(jj)])(ii,1) = meanFiringRate;
        spatialModulation.(['bits_spike_map_' num2str(jj)])(ii,1) = bits_spike;
        spatialModulation.(['bits_second_map_' num2str(jj)])(ii,1) = bits_second;
        spatialModulation.(['sparsity_map_' num2str(jj)])(ii,1) = sparsity;
        spatialModulation.(['selectivity_map_' num2str(jj)])(ii,1) = selectivity;

        spatialModulation.(['PF_position_map_' num2str(jj)])(ii,1) = pf_position;
        spatialModulation.(['PF_peak_map_' num2str(jj)])(ii,1) = pf_peak;
        spatialModulation.(['PF_size_map_' num2str(jj)])(ii,1) = pf_size;
        spatialModulation.(['PF_boundaries_map_' num2str(jj)])(ii,:) = fieldX;
        spatialModulation.(['is_placeField_map_' num2str(jj)])(ii,:) = double(is_placeField);
        
        clear pf_position pf_peak pf_size fieldX is_placeField
        clear meanFiringRate bits_spike bits_second sparsity selectivity T logArg occupancy
    end
end

% 3 Mutual information between space and rate
params.nsmooth = 1;%  # smoothing windows (x)
for jj = 1:length(firingTrialsMap.raster_rate{1})
    for ii = 1:length(firingTrialsMap.raster_rate)
        params.npos = size(firingTrialsMap.raster_rate{ii}{jj},2);   % - # position bins

        raster_rate = firingTrialsMap.raster_rate{ii}{jj}; 
        raster_rate(isnan(raster_rate)) = 0; raster_rate(isinf(raster_rate)) = -1;
        raster_rate(raster_rate==-1) = max(raster_rate(:));
        raster_rate_smooth = smoothdata(raster_rate','gaussian',8)';
        rr{1} = raster_rate_smooth;
        [~, mut_inf] = spatialInfo(rr, params);
        spatialModulation.(['rate_space_MI_map' num2str(jj)])(ii,:) = mut_inf;
        clear raster_rate rr mut_inf
    end
end

% 4 General metrics
% 4.1 Stability (it does not separate between recordings... )
if size(behaviour.masks.recording,1) < size(behaviour.masks.recording,2)
    behaviour.masks.recording = behaviour.masks.recording';
end
list_of_recordings = unique(behaviour.masks.recording);
count2 = 1;
for jj = 1:length(list_of_recordings)
    masks_trials = behaviour.masks.trials(behaviour.masks.recording == list_of_recordings(jj));
    half_trial = round(max(masks_trials)/2);
    
    if ~isfield(behaviour.masks,'direction')
        behaviour.masks.direction = behaviour.masks.arm; 
    end
    if size(behaviour.masks.direction,1) < size(behaviour.masks.direction,2)
        behaviour.masks.direction = behaviour.masks.direction';
    end
    masks_directions = behaviour.masks.direction(behaviour.masks.recording == list_of_recordings(jj));
    list_of_directions = unique(masks_directions);
    count = 1;
    clear idx_second_half idx_first_half
    if size(behaviour.masks.trials,1) < size(behaviour.masks.trials,2)
        behaviour.masks.trials = behaviour.masks.trials';
    end
    for ii = 1:length(list_of_directions)
        idx_first_half{count} = behaviour.masks.trials<half_trial & behaviour.masks.direction == list_of_directions(ii) & behaviour.masks.recording == list_of_recordings(jj);
        idx_second_half{count} = behaviour.masks.trials>half_trial & behaviour.masks.direction == list_of_directions(ii) & behaviour.masks.recording == list_of_recordings(jj);
        count = count + 1;
    end
    % get maps
    clear maps
    count = 1;
    for ii = 1:length(idx_first_half)
        maps{count} = [behaviour.timestamps(idx_first_half{ii}) behaviour.position.lin(idx_first_half{ii})];
        count = count + 1;
        maps{count} = [behaviour.timestamps(idx_second_half{ii}) behaviour.position.lin(idx_second_half{ii})];
        count = count + 1;
    end    
    % get fields
    firingMaps_stability = bz_firingMapAvg(maps, spikes,'saveMat',false,'speedThresh',0.1);
    
    % computing stability
    count = 0;
    for ii = 1:length(list_of_directions)
        stability_temp = ones(size(firingMaps_stability.UID));
        for jj = 1:length(stability_temp)
            stability_temp(jj) = corr(firingMaps_stability.rateMaps{jj}{1+count}', firingMaps_stability.rateMaps{jj}{2+count}','Type','Spearman');
        end
        count = count + 2;
        spatialModulation.(['stability_map' num2str(count2)]) = stability_temp';
        count2 = count2 + 1;
    end
end

if saveMat
    save([basenameFromBasepath(pwd) '.spatialModulation.cellinfo.mat'],'spatialModulation'); 
end

cd(prevPath);
end